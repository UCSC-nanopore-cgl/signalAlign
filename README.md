## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.

### Cheat sheet/Quick Start

### Pre-installation on `toil-box` (if you're using `cgcloud`)
1. `sudo apt-get update && sudo apt-get install zlib1g-dev g++ git`

### Installation:
1. Recursively clone this repo `git clone --recursive -b update_manual https://github.com/UCSC-nanopore-cgl/signalAlign.git`
2. Set required environment variables. Use path to local anaconda or miniconda. `$HOME/anaconda3` is just a common install location for anaconda.  
   `export PATH="$HOME/anaconda3/bin:$PATH"`  
   `export C_INCLUDE_PATH="$HOME/anaconda3/envs/signalalign/include:$C_INCLUDE_PATH"`  
   `export LD_LIBRARY_PATH="$HOME/anaconda3/envs/signalalign/lib:$LD_LIBRARY_PATH"`  
   `export LIBRARY_PATH="$HOME/anaconda3/envs/signalalign/lib:$LIBRARY_PATH"`  
3. Create conda environment   
`conda env create -f requirements.yml python=3.6`
4. Activate conda environment  
`source activate signalalign`
5. Make project  
`make`
6. Add bin to path  
`export PATH=$PATH:$PWD/bin`
7. Test install   
`make test`

### Introduction
Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, then applying a voltage across the membrane and monitoring the ionic current through the nanopore. The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are enzymatically guided through the nanopore. **SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a hierarchical Dirichlet process (HDP) mixture of normal distributions. The HDP models enable mapping of methylated bases to your reference sequence. Instructions for usage including building/training HDP models can be found in the [manual](https://github.com/UCSC-nanopore-cgl/signalAlign/blob/master/Manual.md).

### Requirements
* GCC 4.4.7 or newer (tested on 4.4.7 and 5.0)
* All requirements can be found in `requirements.yml`

