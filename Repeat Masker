###will use the output files created from repeatmodeler to calculate the divergence (age) of TEs in the genome

#!/bin/bash

#SBATCH -p 128x24   # Partition name
#SBATCH -J repeatmasker        # Job name
#SBATCH --mail-user=aymays@ucsc.edu
#SBATCH --mail-type=ALL
#SBATCH -o repeatmasker_aurelia.out    # Name of stdout output file
#SBATCH -N 1       # Total number of nodes requested (128x24/Instructional only)
#SBATCH -n 16        # Total number of mpi tasks requested per node
#SBATCH -t 64:00:00  # Run Time (hh:mm:ss) - 1.5 hours (optional)
#SBATCH --mem=24G # Memory to be allocated PER NODE

#sets language and character settings for the output files, put before code
export LANGUAGE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
export LANG=en_US.UTF-8

#load the repeatmasker module (program)
module load repeatmasker

#run repeatmasker, -lib specifies the library created by repeatmodeler (families.fna file), -pa is the number of threads (power) you want to use to run it, then put the genome fasta file, -a says you want it to do alignment analysis
RepeatMasker -lib aurelia_aurita_db-families.fa -pa 16 aurelia_aurita_genome.fna -a

