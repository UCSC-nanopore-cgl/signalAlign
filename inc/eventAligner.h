//
// Created by Andrew Bailey on 6/8/18.
//  c++ converted to C functions from nanopolish (https://github.com/jts/nanopolish)

#ifndef NANOPORE_RNN_EVENTAILIGNER_H
#define NANOPORE_RNN_EVENTAILIGNER_H
#define MAX_KMER_SIZE 6

#include <stdlib.h>
#include <hdf5.h>
#include "event_detection.h"
#include "stateMachine.h"
#include "scrappie_common.h"
#include "nanopore.h"

typedef struct {
    double start;
    double length;
    double mean;
    double stdv;
    char model_state[MAX_KMER_SIZE + 1];
    int move;
    uint64_t raw_start;
    uint64_t raw_length;
    double p_model_state;
} basecalled_event;

typedef struct {
    size_t n;
    size_t aln_n;
    size_t start;
    size_t end;
    basecalled_event *event;
} basecalled_event_table;

//join directory and file
char* path_join_two_strings(char* directory, char* file_name);

// check if file has correct extension
bool check_file_ext(char* file_path, char* ext);

// write fastqs from directory to outpath
int write_fastqs_to_file(char* fast5_dir, char* output_path);

// get a string dataset
char *fast5_get_string(hid_t hdf5_file, char* path);

// get fastq from fast5
char *fast5_get_fastq(hid_t hdf5_file);

// check if group exists
bool hdf5_group_exists(hid_t hdf5_file, char* path);

// open the file and return the hdf ID
hid_t fast5_open(char* filename);

// write readdb file: sequence_name /t fast5_name /n
int write_readdb_file1(char* fast5_dir, char* output_path);

//write both fastq and readdb file
int write_fastq_and_readdb_file1(char* fast5_dir, char* fastq_output_path, char* readdb_output_path);

typedef struct {
    char* name;
    char* comment;
    char* seq;
    char* qual;
} fastq_entry;

// parse fastq string from fast5 file
fastq_entry* parse_fastq_string(char* fastq_string);

void fastq_entry_destruct(fastq_entry* fastq);

// From scrappie
typedef struct {
    //  Information for scaling raw data from ADC values to pA
    float digitisation;
    float offset;
    float range;
    float sample_rate;
} fast5_raw_scaling;


// close the file
herr_t fast5_close(hid_t hdf5_file);

// get the raw samples from this file
raw_table fast5_get_raw_samples(hid_t hdf5_file, fast5_raw_scaling scaling);
//
// get the name of the raw read in the file (eg Read_1234)
char* fast5_get_raw_read_name(hid_t hdf5_file);

// get the name of the raw read group (eg /Raw/Read/Read_1234)
char* fast5_get_raw_read_group(hid_t hdf5_file);

// Get the identifier of a read from the hdf5 file
char* fast5_get_read_id(hid_t hdf5_file);

// Get the experiment type attribute
char* fast5_get_experiment_type(hid_t hdf5_file);

// Get sample rate, and ADC-to-pA scalings
fast5_raw_scaling fast5_get_channel_params(hid_t hdf5_file);

// Look for a group and if it does not exist, create it
// (group_location = /Old/New) -> creates New group
herr_t fast5_create_group(hid_t hdf5_file, char* group_location);

// Look for a group and if it does not exist, create it
// (group_location = /Old/New/New2/New3) -> creates New/New2/New3 groups
herr_t fast5_create_all_groups(hid_t hdf5_file, char* group_location);

// set basecalled events table
herr_t fast5_set_basecall_event_table(hid_t hdf5_file, char* table_location, basecalled_event_table *et);

// get basecalled events table
herr_t fast5_get_basecall_events(hid_t hdf5_file, char* table_location, basecalled_event *dst_buf);

//get start time of a read
float fast5_get_start_time(hid_t hdf5_file);

//char* fast5_get_fixed_string_attribute(hid_t hdf5_file, const std::string& group_name, const std::string& attribute_name);
char* fast5_get_fixed_string_attribute(hid_t hdf5_file, char* group_name, char* attribute_name);

// get a fast5 float attribute
float fast5_read_float_attribute(hid_t group, const char *attribute);

// use nanopolish's algorithm to estimate the shift and scale
NanoporeReadAdjustmentParameters estimate_scalings_using_mom(stList* kmer_list, StateMachine pore_model, event_table et);

// update the SignalMachine with new NanoporeReadAdjustementParameters
void* update_SignalMachineWithNanoporeParameters(NanoporeReadAdjustmentParameters npp, StateMachine *sM);

// create adaptive banded alignment in C using our model
stList* adaptive_banded_simple_event_align(event_table et, StateMachine *pore_model, stList* kmer_list);
stList* adaptive_banded_simple_event_align2(event_table et, StateMachine *pore_model, stList* kmer_list, bool writeFailedAlignment);

// convert event table into basecalled event table. Assume's start and length are "raw_start" and "raw_length"
basecalled_event_table* event_table_to_basecalled_table(event_table *et, fast5_raw_scaling scaling, float start_time);

// fill basecalled_event_table from alignment
void alignment_to_base_event_map(stList *event_alignment, basecalled_event_table* b_et,
                                 stList *kmer_list, StateMachine *pore_model);

void rna_alignment_to_base_event_map(stList *event_alignment, basecalled_event_table* b_et,
                                     stList *kmer_list, StateMachine *pore_model);

herr_t load_from_raw(char* fast5_file_path, char* templateModelFile, char* sequence, char* path_to_embed, bool rna);
herr_t load_from_raw2(char* fast5_file_path, char* templateModelFile, char* sequence, char* path_to_embed, bool writeFailedAlignment, bool rna);

event_table reverse_events(event_table et);

basecalled_event_table* reverse_basecalled_events(basecalled_event_table *bet);

// struct from nanopolish
struct AlignedPair
{
    int ref_pos;
    int read_pos;
};

// helper functions which deal with the stList of AlignedPairs
void alignedPair_destruct(struct AlignedPair *aligned_pair);

struct AlignedPair *alignedPair_construct(int ref_pos, int read_pos);

// generate list of kmers to align to events.
// If rna, we replace the U's with T's and create 3'-5' kmers in reverse order because we align events in reverse order to
// limit effect of the beginning long stalls during polyA sequencing
//
// DNA - ATGCATGC -> ATGCA, TGCAT, GCATG, CATGC
// RNA - AUGCAUGC -> ACGTA, TACGT, GTACG, CGTAC
stList* build_kmer_list(const char* sequence, int64_t kmer_len, bool rna);



#endif //NANOPORE_RNN_EVENTAILIGNER_H
