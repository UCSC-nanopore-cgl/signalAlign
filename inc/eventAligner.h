//
// Created by Andrew Bailey on 6/8/18.
//  c++ converted to C functions from nanopolish (https://github.com/jts/nanopolish)

#ifndef NANOPORE_RNN_EVENTAILIGNER_H
#define NANOPORE_RNN_EVENTAILIGNER_H

#include <stdlib.h>
#include "eventAligner.h"
#include <hdf5.h>
#include "event_detection.h"
#include "stateMachine.h"
#include "scrappie_common.h"
#include "nanopore.h"

//typedef struct {
//    double scale;
//    double shift;
//    double drift;
//    double var;
//    double scale_sd;
//    double var_sd;
//
//    // derived parameters that are cached for efficiency
//    double log_var;
//    double scaled_var;
//    double log_scaled_var;
//} SquiggleScalings;
//
//
//static SquiggleScalings const SquiggleScalings_default = {
//        .scale = 1.0,
//        .shift = 0.0,
//        .drift = 0.0,
//        .var = 1.0,
//        .scale_sd = 1.0,
//        .var_sd = 1.0,
//        .log_var = 0.0,
//        .scaled_var = 1.0,
//        .log_scaled_var = 0.0
//};
//
//
//
//SquiggleScalings set6_SquiggleScalings(double _shift, double _scale, double _drift, double _var, double _scale_sd,
//                                       double _var_sd);
//
//SquiggleScalings set4_SquiggleScalings(double _shift, double _scale, double _drift, double _var);

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
void fast5_close(hid_t hdf5_file);

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
//
//
// Internal utility functions
//
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

//void load_from_raw(hid_t hdf5_file, StateMachine sM, char* sequence);

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
