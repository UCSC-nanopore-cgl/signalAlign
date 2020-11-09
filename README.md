## SignalAlign

#### MinION signal-level alignment and methylation detection using hidden Markov Models with hierarchical Dirichlet process kmer learning.

### Cheat sheet/Quick Start

### Docker image
ucscbailey/signalalign
### Installation & Requirements
See [dockerfile](Docker/Dockerfile)

### Introduction
Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, then applying a voltage across the membrane and monitoring the ionic current through the nanopore. The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are enzymatically guided through the nanopore. **SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a hierarchical Dirichlet process (HDP) mixture of normal distributions. The HDP models enable mapping of methylated bases to your reference sequence. Instructions for usage including building/training HDP models can be found in the [manual](https://github.com/UCSC-nanopore-cgl/signalAlign/blob/master/Manual.md).

