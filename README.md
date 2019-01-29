## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.

### Cheat sheet/Quick Start

### Docker
You can take a look at the config files which point to data already stored within the Docker image. If you want to run SignalAlign using your own data, you just need to place your data and config into a 
directory and set the config up so that it points to the data correctly once you mount the local directory to the docker image. 


1. `docker run -v /complete/path/to/signalAlign/tests/test_Docker/:/data adbailey4/train_sa_models run --config /data/trainModels-config.json`
2. `docker run -v /complete/path/to/signalAlign/tests/test_Docker/:/data adbailey4/signalalign run --config /data/runSignalAlign-config.json`


### Pre-installation
1. `sudo apt-get update`
2. `sudo apt-get install wget git make g++ zlib1g-dev libbz2-dev liblzma-dev python3-dev python3-setuptools` 

### Installation:
1. Recursively clone this repo `git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git`
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

