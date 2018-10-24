## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.
Documentation is still being worked on, apologies for this.  

### Cheat sheet/Quick Start

### Pre-installation
1. `sudo apt-get update`
2. `sudo apt-get install wget git make g++ zlib1g-dev libbz2-dev liblzma-dev python3-dev python3-setuptools` 

### Getting BWA
1. `git clone https://github.com/lh3/bwa.git`
2. `cd bwa`
3. `make`
4. `export $PATH=$(pwd):$PATH`

### Installation:
1. Clone this repo `git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git && cd signalAlign`
3. Create a python3 virtual environment `virtualenv -p python3 venv && . venv/bin/activate`
4. Compile the executables `make`
5. Test the program `make test`
5. All of the programs can be found in the `bin/` directory

### Introduction
Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, then applying a voltage across the membrane and monitoring the ionic current through the nanopore. The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are enzymatically guided through the nanopore. **SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a hierarchical Dirichlet process (HDP) mixture of normal distributions. The HDP models enable mapping of methylated bases to your reference sequence. Instructions for usage including building/training HDP models can be found in the [manual](https://github.com/ArtRand/signalAlign/blob/master/Manual.md).

### Requirements
* Python 3
    1. H5Py
    2. Numpy
    3. Pandas
    3. Scipy
    4. Pysam
* BWA-MEM (Li H. (2013), instructions can be found (https://github.com/lh3/bwa))
    * Needs to be in path
* GCC 4.4.7 or newer (tested on 4.4.7 and 5.0)


*Code in this repo is based on cPecan (https://github.com/benedictpaten/cPecan)*
