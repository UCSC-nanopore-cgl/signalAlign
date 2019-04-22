## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.

### Cheat sheet/Quick Start

### Docker
You can take a look at the config files which point to data already stored within the Docker image. If you want to run SignalAlign using your own data, you just need to place your data and config into a 
directory and set the config up so that it points to the data correctly once you mount the local directory to the docker image. 


1. `docker run -v /complete/path/to/signalAlign/tests/test_Docker/:/data adbailey4/train_sa_models run --config /data/trainModels-config.json`
2. `docker run -v /complete/path/to/signalAlign/tests/test_Docker/:/data adbailey4/signalalign run --config /data/runSignalAlign-config.json`


### Requirements
* git, gcc, g++, make, zlib, libbz2, liblzma, numpy, setuptools, python 3.5, wget
* Needed in the path
    * [samtools](https://www.biostars.org/p/328831/)
    * [bwa](https://github.com/lh3/bwa)

### Installation:
1. Recursively clone this repo `git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git`
2. Make project  
`make`
3. Add bin to path  
`export PATH=$PATH:$PWD/bin`
4. Test install   
`make test`

### Introduction
Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, then applying a voltage across the membrane and monitoring the ionic current through the nanopore. The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are enzymatically guided through the nanopore. **SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a hierarchical Dirichlet process (HDP) mixture of normal distributions. The HDP models enable mapping of methylated bases to your reference sequence. Instructions for usage including building/training HDP models can be found in the [manual](https://github.com/UCSC-nanopore-cgl/signalAlign/blob/master/Manual.md).

