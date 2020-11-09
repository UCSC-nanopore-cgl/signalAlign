[![Build Status](https://travis-ci.org/UCSC-nanopore-cgl/signalAlign.svg?branch=master)](https://github.com/UCSC-nanopore-cgl/signalAlign)


## SignalAlign Manual

### Introduction
SignalAlign is a hidden Markov model (HMM) software package for aligning the ionic current signal from the Oxford Nanopore Technologies (ONT) MinION to a reference sequence and inferring properties of the reference sequence.

### Installation:
1. Recursively clone this repo `git clone --recursive https://github.com/UCSC-nanopore-cgl/signalAlign.git`
2. Make project  
`make`
3. Add bin to path  
`export PATH=$PATH:$PWD/bin`
4. Test install   
`make test`


### Recommened Workflow
1. Install
2. Extract fastq's and generate a readdb file with mapping between read_id and fast5 file. Note: A sequencing summary file can be passed in as a `readdb` file. 
3. Generate BAM file, preferably with MD field. We have a filtering steps which removes secondary and supplementary reads as well as reads with a mean phred quality score of less than 7. NOTE: Use minimap2 for RNA mapping so BAM format is correct. 
4. If you are interested in specific regions of the genome or some other read mapping based information I would recommend filtering your BAM to specified regions.
5. To speed up computation, I would also recommend using `filterReads` to move fast5 files into another directory so the program does not have to loop through all of the failed reads. 
6. A config file is setup in tests/ to run from the signalAlign home directory
```bash
runSignalAlign run --config tests/runSignalAlign-config.json &> log_file.txt
```

### Description of programs
* filterReads
    * Filter out reads based on mapping information and base quality score if MD flag is present in the BAM file. 
* runSignalAlign
    * Aligns ionic current events from a directory of basecalled MinION reads (.fast5) to a reference sequence. With appropriate inputs, characters in the reference sequence can be flagged as _ambiguous_, meaning multiple bases will be probabilistically aligned to a position. Right now, the program the program supports aligning cytosine variants (5-mC and 5-hmC) and adenine variants (6-mA) to a position.
* trainModels
    * Trains the transitions and/or emissions of the HMM. Uses a directory of basecalled reads and a reference sequence to learn transition parameters to the model. Once enough assignments have been accumulated it will (optionally) rebuild a hierarchical Dirichlet process (HDP) from these assignments.

### Nucleotide Encodings
| Character | Nucleotide         |
|:---------:|--------------------|
| A         | Adenine            |
| T         | Thymine/Uracil     |
| G         | Guanine            |
| C         | Cytosine           |
|---------|--------------------|
| F         | 6-Methyladenine  |
| E         | 5-Methylcytosine   |
| O         | 5-Hydroxymethylcytosine   |
| J         | Bromodeoxyuridine   |
| p         | Pseudouridine   |
| b         | 7-Methylguanosine   |
| d         | 2-methylguanosine   |
| e         | n4-methylcytidine   |
| h         | 2'-O-methyluridine   |
| i         | n6,n6-dimethyladenosine   |

| Ambig Character | Other Characters |
|:---------:|--------------------|
| R         |AG |
| Y         |CT |
| S         |CG |
| W         |AT |
| K         |GT |
| M         |AC |
| B         |CGT |
| D         |AGT |
| H         |ACT |
| V         |ACG |
| X         |ACGT |
| L         |CEO |
| P         |CE |
| Q         |AI |
| f         |AF |
| U         |ACEGOT |
| Z         |JT |
| j         |Tp |
| k         |Gb |
| l         |Gd |
| m         |Ce |
| n         |Th |
| o         |Ai |


### extract
Will search subdirectories if recursive flag is set
#### Input
* `-d`: fast5 dir
* `-o`: name of output file 
* `-r`: recursive

### filterReads
If you have several sub_directories of fast5 data, then use the recursive flag and I would recommend setting copy_dir_structure so you know the original locations of the moved fast5s
#### Input
*  `--alignment_file` Bam file with all alignment data  
* `--fast5_dir` Directory of all fast5 files  
* `--readdb` Path to readdb file  
* `--pass_output_dir` Location where all pass reads will be moved  
* `--quality_threshold` Minimum average base quality threshold. Default = 7  
* `--recursive` Search directory recursively to find fast5 files  
* `--copy_dir_structure`  Step into directory and copy directory structure for output files  
* `--trim`    Only move as many files which contain a total of some number of bases set by trim.  
* `--jobs`    Number of jobs to start if copy_dir_structure is set  
* `--debug`   Will run copy_dir_structure with only one job and fail if errors arise


### runSignalAlign-update
Both command line and config file options are available. The config file is almost identical to the trainModels config.  
*  `runSignalAlign run --config config.json` if you want to use a config file  
* `runSignalAlign -d /path/to/minionReads/ -r /path/to/reference.fasta -o /path/to/output/ --in_tempalte_hmm /path/to/model &> log_file.txt` For command line usage

#### Input
_Required_
* A directory of MinION reads, `-d`, (*.fast5) that have been basecalled. Currently, the following versions of Metrichor basecalled files are supported: 1.15.0, 1.19.0, 1.20.0, 1.22.2, 1.22.4. If you have a more recent or unsupported version open an issue and I'll modify the program (or feel free to implement it yourself and issue a pull request!). 
* A reference sequence, `-r`, in FASTA format.
* Output location, `-o`, a path to use as working directory. A new directory will be made here, so the program won't pollute this directory.  
* A file containing trained HMM transitions parameters `-T` (template HMM) `-C` (complement HMM).

_HIGHLY Recommended_ 
* Pass BAM file into `--alignment_file`
* Set the `--path_to_bin` so the program can find all of the executables  
* Pass in a readdb or sequencing summary file `--readdb`
* A file containing a HDP model or other emissions (normal distributions) parameters `-tH` (template HDP) `-cH` (complement HDP).
* Set `--filter_reads` to avoid low quality reads 

_Optional_
* Target regions file. Only reads that map to these regions will follow on to event-alignment `-q`
* Ambiguity positions file, `-p` flags positions to be aligned to multiple bases (variant calling).

#### Options and flags
* `--file_directory`, `-d` directory with MinION fast5 reads to align  
* `--in_tempalte_hmm`, `-T` template HMM parameters file
* `--in_complement_hmm`, `-C` complement HMM parameters file 
* `--template_hdp`, `-tH` template HDP model file
* `--complement_hdp`, `-cH` complement HDP model file 
* `--degenerate`, `-x` nucleotide options for degenerate or _ambiguous_ positions. `m6a` = {AF}, `variant` = {A,C,G,T} `cytosine2` = {CE} `cytosine3` = {CEO} `adenosine` = {AI}. **n.b.** E = 5-methylcytosine, O = 5-hydroxymethylcytosine, I = 6-methyladenine
* `--stateMachineType`, `-smt` HMM to use. Options: `threeState` and `threeStateHdp`. Default: `threeState`.
* `--file_of_files`, `-fofn` a file containing the absolute path to files to align with, one file path per line
* `--threshold`, `-t`. Minimum posterior match probability threshold (matches below this threshold will not be tabulated). Default: 0.01.
* `--diagonalExpansion`, `-e` Mumber of diagonals to expand around each anchor, Default: 50.
* `--constraintTrim`, `-m` Amount to remove from an anchor constraint. Default: 14.
* `--target_regions`, `-q` File containing target regions to align to, if the read doesn't get mapped to a region in this file the signal-level alignment will not precede. The format is `start \t end \t kmer` where `start` and `end` are the genomic coordinates and the `kmer` is the sequence at those coordinates. The `kmer` is mostly historical, so it doesn't have to be perfect. See below for details.
* `--ambiguity_positions`, `p` file containing positions to flag as ambiguous see below for details.
* `--jobs, -j` Number of jobs to run concurrently, Default: 4.
* `--nb_files`, `-n` Maximum number of files to align, will randomly select if you point to a directory/fofn with more files than this so if you want to make sure you use all of the files set a large number (or count the lines in your fofn). Default: 500.
* `--ambig_char`, `-X` in development, will be for specifying specific subsets of ambiguous positions **leave as default**
* `--output_format`, `-f` format of output (this is different than the summary that goes to `stdout`). Options are: `full`, `variantValler`, and `assignments`, see below for description.
* `--output_location`, `-o` place to put results, it will make a new working directory in this location
* `--forward_ref` forward reference sequence for SignalAlignment align to, in FASTA  
* `--backward_ref` backward reference sequence for SignalAlignment align to, in FASTA
* `--alignment_file` BAM file of alignments if FAST5s do not have fasta info or events    
* `--bwa_reference`, `-r` Reference sequence required for generating guide alignment  
* `--motifs` Motif find and replace must be in specific list within a list format. eg: [['CCAGG', 'CEAGG'], ['CCTGG', 'CETGG']]  
* `--embed`   Embed full output into fast5 file  
* `--event_table` Specify path withing Fast5 event table  
* `--force_kmer_event_alignment` If passed, force SignalAlign to infer kmer to event alignment. Must include alignment_file  
* `--allow_unsupported_nanopore_read_versions` Will attempt to complete execution with unsupported nanopore read versions  
* `--filter_reads`        Will filter reads out if average fastq quality scores are below 7.  
* `--path_to_bin` Path to bin to find signalMachine  
* `--readdb`        Path to readdb file or sequencing summary file for easy filtering 
* `--keep_tmp_folder`     Keep the temporary folder with files fed into SignalMachine
* `--recursive`           Recursively search top directory for fast5 files

#### Output

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

### trainModels
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


#### Mapping methylation by using substitution and target files
**signalAlign** uses two files to map methylation in a reference sequence. First it uses a *target file* that specifies positions in the reference, only reads that map to regions specified in the file will be aligned on the signal-level. The *target file* has format: `start \t end \t sequence`. An example can be seen here:`/signalAlign/tests/test_regions/test_sites_bal_1.tgt`. The second file is a *label file* that tells signalAlign which bases to flag as ambiguous. This file has the format `X \t position \t ... \n` one line for both the template and complement strand. An example is at `/signalAlign/tests/test_regions/test_labels_bal_1.tsv`.

An example command that would produce methylation call probabilities in _E. coli_
(**n.b.** to run this example, download the HDP models from the [HDP_models](https://github.com/ArtRand/HDP_models) repo)
```bash
./runSignalAlign \
> -d ../tests/minion_test_reads/ecoli/ \
> -r ../tests/test_sequences/E.coli_K12.fasta \
> -T ../../HDP_models/ecoli_r7.3_models/template_trained.hmm \
> -C ../../HDP_models/ecoli_r7.3_models/complement_trained.hmm \
> -tH ../../HDP_models/ecoli_r7.3_models/template.multisetPrior2.nhdp \
> -cH ../../HDP_models/ecoli_r7.3_models/complement.multisetPrior2.nhdp \
> -x cytosine2 \
> -smt=threeStateHdp \
> -q ../tests/test_regions/test_sites_bal_1.tgt \
> -p ../tests/test_regions/test_labels_bal_1.tsv \
> -f variantCaller \
> -o ../../ \
> 2> ../../a.err
```

