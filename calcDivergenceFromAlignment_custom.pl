#!/usr/bin/perl
use strict;
use Getopt::Long;
use FileHandle;
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::Bin/..";
use RepeatMaskerConfig;
use SearchResult;
use CrossmatchSearchEngine;

my $Version = $RepeatMaskerConfig::VERSION;
my $DEBUG   = 0;
my @getopt_args = ('-version', '-noCpGMod', '-a=s', '-s=s', '-tbl=s');
my %options = ();
Getopt::Long::config("noignorecase", "bundling_override");
unless ( GetOptions(\%options, @getopt_args) ) {
    die "Usage: calcDivergenceFromAlign_custom.pl -s <divsum> [-a <new.align>] [-tbl <RepeatMasker.tbl>] file.align[.gz]\n";
}

my %classMap;
if ($options{'tbl'}) {
    my $tbl = $options{'tbl'};
    open my $tfh, "<", $tbl or die "Cannot open $tbl: $!";
    while (<$tfh>) {
        next if /^#/ || /^\s*$/;
        my @f = split;
        next if @f < 2;
        # from .out
        if ($f[-2] !~ /bp/ && $f[-1] =~ m{[A-Za-z]+/[A-Za-z]+}) {
            my $repname  = $f[-2];
            my $repclass = $f[-1];
            $classMap{$repname} = $repclass;
        }
        # from .tbl summary
        elsif ($f[0] !~ /bases|length|number/i && $f[0] ne "File" && $f[0] ne "===") {
            my $repclass = $f[0];
            $repclass =~ s/://;
            $classMap{$repclass} = $repclass;
        }
    }
    close $tfh;
    print STDERR "Loaded ", scalar(keys %classMap), " TE class mappings from $options{tbl}\n";
}

my $alignFile = $ARGV[0] or die "Missing alignment file!\n";
my $maxDiv    = 70;
my %repeatMuts;
my %classDivWCLen;
my $cntAlign = 0;

my $searchResultsFH = new FileHandle;
if ( $alignFile =~ /\.gz$/ ) {
    open $searchResultsFH, "gunzip -c $alignFile |" or die "Cannot gunzip $alignFile: $!\n";
} else {
    open $searchResultsFH, "<", $alignFile or die "Cannot open $alignFile: $!\n";
}

if ($options{'s'}) {
    open SOUT, ">", $options{'s'} or die "Cannot open $options{'s'} for writing\n";
    print SOUT "Weighted average Kimura divergence for each repeat family\n";
    print SOUT "Class\tRepeat\tabsLen\twellCharLen\tKimura%\n";
    print SOUT "-----\t------\t------\t-----------\t-------\n";
}

my $outAlign = 0;
if ($options{'a'}) {
    open COUT, ">", $options{'a'} or die "Cannot open $options{'a'} for writing\n";
    $outAlign = 1;
}

CrossmatchSearchEngine::parseOutput(
    searchOutput => $searchResultsFH,
    callback     => sub {
        my $result = shift;
        return unless $result;

        my $subjName = $result->getSubjName();
        my ($hitname, $class);

        # Extract class from align label
        if ($subjName =~ /(\S+)[#\|](\S+)/) {
            ($hitname, $class) = ($1, $2);
        } else {
            $hitname = $subjName;
            $class   = "Unknown";
        }

        # Reassign from mapping if possible
        if (exists $classMap{$hitname}) {
            $class = $classMap{$hitname};
        } elsif (exists $classMap{$class}) {
            $class = $classMap{$class};
        }

        $class = "Unclassified" if $class eq "" || $class =~ /Unknown/i;

        my $queryStart = $result->getQueryStart();
        my $queryEnd   = $result->getQueryEnd();
        my $alen       = $queryEnd - $queryStart + 1;
        my $wellCharBases = $alen - int($alen * ($result->getPctInsert() / 100));

        my $div = $result->getPctKimuraDiverge();
        if ($div eq "") {
            my ($d2, undef, undef, undef, undef) = $result->calcKimuraDivergence(divCpGMod => 1);
            $div = sprintf("%4.2f", $d2);
            $result->setPctKimuraDiverge($div);
        }

        $repeatMuts{$class}->{$hitname}->{'sumdiv'}      += $div * $wellCharBases;
        $repeatMuts{$class}->{$hitname}->{'wellCharLen'} += $wellCharBases;
        $repeatMuts{$class}->{$hitname}->{'absLen'}      += $alen;

        my $key = "$class " . int($div);
        $classDivWCLen{$key} += $wellCharBases;

        print COUT $result->toStringFormatted(SearchResult::AlignWithQuerySeq)."\n" if $outAlign;
        print STDERR "." if (++$cntAlign % 5000 == 0);
    }
);

# Write divsum
if ($options{'s'}) {

    # Weighted table
    foreach my $class (sort keys %repeatMuts) {
        foreach my $id (sort keys %{ $repeatMuts{$class} }) {
            my $wlen = $repeatMuts{$class}{$id}{'wellCharLen'} || 0;
            my $kimura = ($wlen > 0)
                ? sprintf("%4.2f", $repeatMuts{$class}{$id}{'sumdiv'} / $wlen)
                : "----";
            $kimura = $maxDiv if ($kimura ne "----" && $kimura > $maxDiv);
            print SOUT join("\t", $class, $id,
                            $repeatMuts{$class}{$id}{'absLen'} || 0,
                            $wlen, $kimura), "\n";
        }
    }

    print SOUT "\nCoverage for each repeat class and divergence (Kimura)\n";
    print SOUT "Div ";
    foreach my $class ( sort keys %repeatMuts ) {
        print SOUT "$class ";
    }
    print SOUT "\n";

    for (my $j = 0; $j <= $maxDiv; ++$j) {
        print SOUT "$j ";
        foreach my $class ( sort keys %repeatMuts ) {
            my $label = "$class $j";
            $classDivWCLen{$label} = 0 unless $classDivWCLen{$label};
            print SOUT "$classDivWCLen{$label} ";
        }
        print SOUT "\n";
    }

    close SOUT;
}

print STDERR "\nDone. Wrote divergence summary to $options{s}\n" if $options{'s'};
exit 0;
