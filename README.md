# Overview

Manuscript: 
## Platforms used
  - UCSC Hummingbird Computing cluster in command line
  - RStudio

## Packages and dependencies
  - repeatmodeler v2.0.5
  - DeepTE

### Use conda to install (per https://github.com/LiLabAtVT/DeepTE)
  - conda create -n py36 python=3.6
  - conda activate py36
  - conda install tensorflow-gpu=1.14.0
  - conda install biopython
  - conda install keras=2.2.4
  - conda install numpy=1.16.0
  - pip install sklearn

## Methodology
We downloaded genomes from NCBI genome bank and processed them through multiple functions of the repeatmodeler module (repeatmodeler, repeatmasker, calcdiv, etc). 
 1. Use repeatmodeler to build database for each genome
 2. merge libraries together using cat (families file)
 3. remove redundancies using CD-hit
 4. extract out Unknown TEs with Seqkt
 5. Annotate unknowns with  DeepTE
 7. Merge all libraries and use as input for repeatmasker and downstream analysis
 8. Run each genome through repeatmasker with custom Medusozoa library
 9. Format .divsum file and run through calcdiv using custom python script
 10. Create repeatlandscape in Rstudio (per https://github.com/cejuliano/brown_hydra_genomes/blob/main/02_repeatMasking/04_visRep/kimuraPlot.R)

# Data availability

