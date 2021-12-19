[![Build Status](https://travis-ci.org/UCSC-nanopore-cgl/signalAlign.svg?branch=master)](https://github.com/UCSC-nanopore-cgl/signalAlign)

# SignalAlign

## Introduction

Nanopore sequencing is based on the principal of isolating a nanopore in a membrane separating buffered salt solutions, 
then applying a voltage across the membrane and monitoring the ionic current through the nanopore. 
The Oxford Nanopore Technologies (ONT) MinION sequences DNA by recording the ionic current as DNA strands are 
enzymatically guided through the nanopore. 
**SignalAlign** will align the ionic current from the MinION to a reference sequence using a trainable hidden Markov model (HMM). 
The emissions model for the HMM can either be the table of parametric normal distributions provided by ONT or a 
hierarchical Dirichlet process (HDP) mixture of normal distributions. 
The HDP models enable mapping of methylated bases to your reference sequence. 

###### Description of programs

* runSignalAlign
  * Aligns ionic current events from a directory of basecalled MinION reads (.fast5) to a reference sequence. With appropriate inputs, characters in the reference sequence can be flagged as _ambiguous_, meaning multiple bases will be probabilistically aligned to a position. Right now, the program the program supports aligning cytosine variants (5-mC and 5-hmC) and adenine variants (6-mA) to a position.
* trainModels
  * Trains the transitions and/or emissions of the HMM. Uses a directory of basecalled reads and a reference sequence to learn transition parameters to the model. Once enough assignments have been accumulated it will (optionally) rebuild a hierarchical Dirichlet process (HDP) from these assignments.


## Installation:

Given the installation is often long, tedious and somewhat fragile. So, we recommend using Docker to run signalAlign.

#### Docker image

`docker run -it -v "$(pwd)":/data ucscbailey/signalalign:latest run --config tests/docker_runSignalAlign-config.json`

### Installation:

There is an installation script which works and has been tested on a clean ubuntu18.08 server from aws. 

#### Requirements

  * Cmake v3.17.0
  * htslib v1.9
  * boost v1.69.0
  * hdf5 v1.10.4
  * python 3.7
  * bwa
  * samtools
  * [vbz_compression](https://github.com/nanoporetech/vbz_compression.git) 
(This became a requirement recently for accessing new fast5 files)


1. install embed_fast5 
```
git clone --recursive https://github.com/adbailey4/embed_fast5.git
cd embed_fast5
python3.7 -m pip install .
python3.7 -m pytest
```

2. Install signalAlign
```
git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git
cd signalAlign && mkdir -p build && cd build
cmake .. -DCMAKE_INSTALL_PREFIX=. -DCMAKE_VERBOSE_MAKEFILE=ON -DBUILD_SHARED_LIBS=ON -DCMAKE_BUILD_TYPE=RELEASE
make -j 4
cd ..
python3.7 -m pip install .
python3.7 -m pytest
```

### Recommened Workflow

1. Install
2. Generate a BAM file, preferably with MD field. We have a filtering steps which removes secondary and supplementary reads as well as reads with a mean phred quality score of less than 7. NOTE: Use minimap2 for RNA mapping so BAM format is correct. If you are interested in specific regions of the genome or some other subsets of reads, I would recommend filtering your BAM to specified regions.
3. Split Multi-fast5's into single fast5s. Yes I know this is tedious and disk intensive but we have no plans on supporting multi-fast5s. [ONT's](https://github.com/nanoporetech/ont_fast5_api) software is recommended.
4. Generate an index file with mapping between read_id and fast5 file. We use nanopolish's index and can be accessed via 'embed_main index'.
5. Generate positions file (see [Positions File Specification](#positions-file-specification))
6. Create a config file. An example config file is located in [tests](tests/runSignalAlign-config.json).
```
runSignalAlign.py run --config tests/runSignalAlign-config.json
```
7. Run signalAlign 

## Signalalign Config Specifications

`{`  
`"signal_alignment_args": {`  
`"target_regions": null, `deprecated  
`"track_memory_usage": false,`deprecated   
`"threshold": 0.01,` probablity threshold for outputing an alignment between kmer and event  
`"event_table": null,` specify location of event table in fast5 (deprecated)  
`"embed": false,` embed signalalign data into fast5 (not recommened)  
`"delete_tmp": true,` delete temorary folders (highly recommened)  
`"output_format": "full"` specify output format "full", "assignments", "variantCaller" (see [Output Formats](#output-formats))  
`},`  
`"samples": [`  
`{`  
`"positions_file": "path/to/file.positions",` locations to call mods (see [Positions File Specification](#positions-file-specification))  
`"fast5_dirs": ["/path/to/fast5/"],`path to fast5 files  
`"bwa_reference": "/path/to/reference.fa",` path to reference   
`"fofns": [],` deprecated  
`"readdb": "/path/to/index.readdb",` path to readdb/index file  
`"fw_reference": null,` deprecated    
`"bw_reference": null,` deprecated   
`"kmers_from_reference": false,` deprecated  
`"motifs": null,` deprecated  
`"name": "wild_type",` name of sample  
`"alignment_file": "/path/to/alignment.bam",` path to alignment file  
`"recursive": false,` deprecated  
`},`
`],`  
`"path_to_bin": "/path/to/signalalign/bin",` path to signalAlign bin where executables are stored  
`"complement_hdp_model": null,` path to complement_hdp_model   
`"template_hdp_model": null,` path to template_hdp_model   
`"complement_hmm_model": null,` path to complement_hmm_model  
`"template_hmm_model": "path/to/model.model",` path to template_hmm_model  
`"job_count": 96,` number of threads / jobs to use  
`"debug": false,` if things are failing set this to true  
`"two_d": false,` if 2d reads set to true  
`"output_dir": "/path/to/output",` path to output directory   
`"perform_kmer_event_alignment": true,` force pre alignment between kmers and events   
`"overwrite": true,` overwrite previous runs  
`"rna": true,` set to true if rna  
`"ambig_model": path/to/ambig.model",` ambiguous model (see [Ambiguous Model Specification](#ambiguous-model-specification)   
`"built_alignments": null,` deprecated  
`}`  

#### Output Formats

There are three output formats. `full`, `variantCaller`, and `assignments`. Each read that aligns to the target regions (if you specified that option) will come as a separate file.

`full` has the following tab-separated-format:

| Contig | Reference Index | Reference k-mer | Read File | Strand | Event Index | Event Mean | Event Noise | Event Duration | Aligned k-mer | Scaled Mean Current | Scaled Noise | Posterior Probability | Descaled Event Mean | Model (ONT) Mean | Path k-mer |
|--- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |

`variantCaller` has the following tab-separated format:

| Event Index | Reference Position | Base | Posterior Probability | Strand | Forward Mapped | Read File |
|--- | --- | --- | --- | --- | --- | --- |

finally, `assignments` has the following tab-separated format:

| k-mer | Read File | Descaled Event Mean | Posterior Probability | 
|--- | --- | --- | --- |



####Positions File Specification

Tsv file with the columns seen below. Used to change the `bwa_reference` characters to ambiguous characters defined
in both the [ambig model](#ambiguous-model-specification) and hmm/hdp model. see [Example](tests/test_position_files/CCWGG_ecoli_k12_mg1655.positions)

| contig (str) | 0_based_ref_index (int)| strand (str)| ref_char (str)| replacement_char (str)|
|--- | --- | --- | --- | --- |


####Ambiguous Model Specification

Tsv file with the columns seen below. If a character in the `find` column is found in the reference, the model
will create branch points for each of the `replace` characters. This model allows you to specify a set of arbitrary characters
to represent a set of arbitrary modifications or variants. Note: Must be consistent with a HMM model.

| find (str) | replace (str) |
|--- | --- |


## trainModels Config Specifications

#### Input

_Required_
* A config file. The default/test model can be found at `tests/trainModels-config.json`

#### Options

`{  `  
`"signal_alignment_args": {` SignalAlign arguments to be used for each sample  
`"target_regions": null,` If you want to specify regions of the genome to process
`"track_memory_usage": false,`  Option to print memory usage statistics   
`"threshold": 0.01,` Minimum output threshold for SignalAlign     
`"event_table": false,` Select specific path to Basecalled Event table    
`"embed": false,` If set, will embed output into fast5 and run the maximum expected accuracy algorithm    
`"delete_tmp": false` Will delete temporary folders after each run. (False is recommended so the temp files do not need to be generated each round of training)    
`},`  
`"samples": [` List of samples to process. If you have a control and an experiment you can create 2 separate samples  
`{`  
`"positions_file": null,` Way of selecting find and replace characters for reference genome    
`"fast5_dirs": ["./tests/minion_test_reads/one_R9_canonical_ecoli"],` List of directories to search for Fast5s  
`"bwa_reference": "./tests/test_sequences/E.coli_K12.fasta",` Path to unedited genome  
`"readdb": null,`  read_id and fast5 path mapping   
`"fw_reference": null,` Forward reference with all characters replaced correctly  
`"bw_reference": null,` Backward reference with all characters replaced correctly  
`"kmers_from_reference": false,` Use only kmers from reference (Recommended when not using Human Genome)  
`"motifs": null,` Motif find and replace must be in specific list within a list format. eg: [['CCAGG', 'CEAGG'], ['CCTGG', 'CETGG']]   
`"name": "canonical",` Name of your sample  
`"probability_threshold": 0.8,` Minimum threshold for including event-kmer mapping in training data for HDP
`"number_of_kmer_assignments": 10,` Number of assignments to create per kmer  
`"alignment_file": null,` BAM file of aligned reads  
`"recursive": false,` Search fast5_dirs recursively  
`"assignments_dir": null` Directory where assignments have been created. (Used for traingin HDP)  
`}`  
`],`  
`"hdp_args": {`  
`"grid_start": 30.0,` Minimum pA coverage by HDP distribution  
`"grid_end": 120.0,`  Maximum pA coverage by HDP distribution  
`"grid_length": 1200,` Number of points that cover the distribution  
`"base_alpha": 1.0,` HDP priors  
`"base_beta": 1.0,`  
`"base_gamma": 1.0,`  
`"middle_alpha": 1.0,`  
`"middle_beta": 1.0,`  
`"middle_gamma": 1.0,`   
`"leaf_alpha": 1.0,`    
`"leaf_beta": 1.0,`   
`"leaf_gamma": 1.0,`   
`"thinning": 100,`    
`"gibbs_samples": 1000,` Number of draws to sample to create distribution  
`"burnin_multiplier": 32,` Helps calculate number of burnin samples to make  
`"hdp_type": "singleLevelFixedCanonical",`   
`"threshold": 0.0,`   
`"number_of_assignments": 100,` Number of assignments for each kmer  
`"built_alignments": null` If Alignments were already built   
`},`  
`"transitions_args": {`  
`"training_bases": 10000,` Number of bases to use for getting transition expectations   
`"iterations": 2,` number of iterations  
`"test": false` If set to true, will check if total probablity increases with each iteration     
`},`  
`"training": {`  
`"transitions": false,`   
`"normal_emissions": true,`  
`"hdp_emissions": true,`  
`"expectation_maximization": false,`  
`"em_iterations": 3`  
`},`  
`"path_to_bin": "./bin",`  
`"complement_hdp_model": null,`  
`"template_hdp_model": null,`  
`"complement_hmm_model": "./models/testModelR9_5mer_acgt_complement.model",`  
`"template_hmm_model": "./models/testModelR9_5mer_acgt_template.model",`  
`"job_count": 2,` number of processes to run  
`"debug": false,` Will fail if any reads fail  
`"two_d": false,`   
`"output_dir": "./",`   
`"constraint_trim": null,` Number of anchor pairs to trim on either side of an anchor block  
`"diagonal_expansion": null,` Size of diagonal to calculate posterior probabilities  
`"traceBackDiagonals": 100,` Number of diagonals to calculate before using an estimate for total probability   
`"filter_reads": 7` Minimum average phred quality score for a read  
}

##### HDP_TYPES_ACEGOT

    ("singleLevelFixed", 0),
    ("singleLevelPrior", 1),
    ("multisetFixed", 2),
    ("multisetPrior", 3),
    ("compFixed", 4),
    ("compPrior", 5),
    ("middleNtsFixed", 6),
    ("middleNtsPrior", 7),
    ("groupMultisetFixed", 8),
    ("groupMultisetPrior", 9),


##### HDP_TYPES_1D

    ("singleLevelPrior2", 10),
    ("multisetPrior2", 11),
    ("singleLevelFixedCanonical", 14),
    ("singleLevelFixedM6A", 15),


##### HDP_TYPES_ACEGT

    ("singleLevelPrior2", 10),
    ("multisetPrior2", 11),


##### HDP_TYPES_ACGT

    ("singleLevelFixedCanonical", 14)

##### HDP_TYPES_ACEGIT

    ("multisetPriorEcoli", 12),
    ("singleLevelPriorEcoli", 13),

##### HDP_TYPES_ACFGT

    F = m6A
    ("singleLevelFixedM6A", 15),

#### Output

* `template_trained.hmm` template HMM with trained parameters
* `complement_trained.hmm` complement HMM with trained parameters if 2D
  If HDPs were used, a copy of the input HDP will also be here, they have `.nhdp` suffix.



