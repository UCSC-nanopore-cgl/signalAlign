## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.

### Cheat sheet/Quick Start

### Docker
You can take a look at the config files which point to data already stored within the Docker image. If you want to run SignalAlign using your own data, you just need to place your data and config into a 
directory and set the config up so that it points to the data correctly once you mount the local directory to the docker image. 


1. `docker run -v /complete/path/to/signalAlign/tests/test_Docker/:/data adbailey4/train_sa_models run --config /data/trainModels-config.json`
2. `docker run -v /complete/path/to/signalAlign/tests/test_Docker/:/data adbailey4/signalalign run --config /data/runSignalAlign-config.json`


### Requirements
* git-2.17.1, gcc-5, g++-5, make-4.1, zlib-1:1.2.11.dfsg-0ubuntu2, libbz2-1.0.6-8.1ubuntu0.2, liblzma-5.2.2-1.3, numpy-1.16.4, setuptools-41.0.1, python-3.5, wget-1.19.4-1ubuntu2.2
* Note: If you have newer versions of any of these installed and it does not work please file an issue. 
* Needed in the path
    * [samtools](https://www.biostars.org/p/328831/)
    * [bwa](https://github.com/lh3/bwa)
* If for some reason your server does not have pip3 associated with python3.5 you can make your own local pip
    * `wget https://bootstrap.pypa.io/get-pip.py`
    * `python3 get-pip.py --user`

### Installation:
0. Sometimes it is necessary to install pysam and cython before running the pip install command so before installing  
    * `pip install pysam`
    * `pip install cython`
1. Recursively clone this repo `git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git`
2. Make project  
`make`
3. Install all required python packages.  
`pip install -e . --user` 
4. Add bin to path  
`export PATH=$PATH:$PWD/bin`
5. Test install   
`make test`

### Introduction
Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, then applying a voltage across the membrane and monitoring the ionic current through the nanopore. The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are enzymatically guided through the nanopore. **SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a hierarchical Dirichlet process (HDP) mixture of normal distributions. The HDP models enable mapping of methylated bases to your reference sequence. Instructions for usage including building/training HDP models can be found in the [manual](https://github.com/UCSC-nanopore-cgl/signalAlign/blob/master/Manual.md).

