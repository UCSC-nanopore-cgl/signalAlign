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


// open the file and return the hdf ID
hid_t fast5_open(char* filename);


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

// set an events table
void* fast5_set_event_table(hid_t hdf5_file, char* table_name, event_table *et);

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
NanoporeReadAdjustmentParameters estimate_scalings_using_mom(const char* sequence, StateMachine pore_model, event_table et);

// update the SignalMachine with new NanoporeReadAdjustementParameters
void* update_SignalMachineWithNanoporeParameters(NanoporeReadAdjustmentParameters npp, StateMachine *sM);

// create adaptive banded alignment in C using our model
stList* adaptive_banded_simple_event_align(event_table et, StateMachine *pore_model, char* sequence);

// convert event table into basecalled event table. Assume's start and length are "raw_start" and "raw_length"
basecalled_event_table* event_table_to_basecalled_table(event_table *et, fast5_raw_scaling scaling, float start_time);

// fill basecalled_event_table from alignment
void alignment_to_base_event_map(stList *event_alignment, basecalled_event_table* b_et,
                                                    char *sequence, StateMachine *pore_model);

herr_t load_from_raw(char* fast5_file_path, char* templateModelFile, char* sequence, char* path_to_embed);


// struct from nanopolish
struct AlignedPair
{
    int ref_pos;
    int read_pos;
};

// helper functions which deal with the stList of AlignedPairs
void alignedPair_destruct(struct AlignedPair *aligned_pair);

struct AlignedPair *alignedPair_construct(int ref_pos, int read_pos);

#endif //NANOPORE_RNN_EVENTAILIGNER_H
