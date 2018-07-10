//
// Created by Andrew Bailey on 6/8/18.
//



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include "hdf5_hl.h"
#include <scrappie_structures.h>
#include <eventAligner.h>

#include "signalMachineUtils.h"

#define RAW_ROOT "/Raw/Reads/"
//#define DEBUG_FAST5_IO 1
//#define DEBUG_PRINT_STATS 1

hid_t fast5_open(char* filename)
    {
    hid_t hdf5file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    return hdf5file;
    }


void* fast5_close(hid_t hdf5_file)
{
    H5Fclose(hdf5_file);
}



char* fast5_get_raw_read_name(hid_t hdf5_file)
{
    // This code is From scrappie's fast5_interface

    // retrieve the size of the read name
    ssize_t size =
            H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, NULL,
                               0, H5P_DEFAULT);

    if (size < 0) {
    return "";
    }

    // copy the read name out of the fast5
    char* name = (char*)calloc(1 + size, sizeof(char));
    H5Lget_name_by_idx(hdf5_file, RAW_ROOT, H5_INDEX_NAME, H5_ITER_INC, 0, name, 1 + size, H5P_DEFAULT);

    return name;
}

char* fast5_get_raw_read_group(hid_t hdf5_file)
{
    char* read_name = fast5_get_raw_read_name(hdf5_file);

    if (read_name != "") {
        char* read_group = stString_concat(RAW_ROOT, read_name);
        return read_group;
    } else {
        return "";
    }
}


char* fast5_get_read_id(hid_t hdf5_file)
{
    int ret;
    hid_t read_name_attribute, raw_group, attribute_type;
    size_t storage_size = 0;
    char* read_name_str = NULL;

    char* out = "";

    // Get the path to the raw read group
    char* raw_read_group = fast5_get_raw_read_group(hdf5_file);
    if(raw_read_group == "") {
        return out;
    }

    return fast5_get_fixed_string_attribute(hdf5_file, raw_read_group, "read_id");
}

raw_table fast5_get_raw_samples(hid_t hdf5_file, fast5_raw_scaling scaling)
{
    float* rawptr = NULL;
    hid_t space;
    hsize_t nsample;
    herr_t status;
    float raw_unit;
    raw_table rawtbl = { 0, 0, 0, NULL };

    // mostly from scrappie
    char* raw_read_group = fast5_get_raw_read_group(hdf5_file);

    // Create data set name
    char* signal_path = stString_concat(raw_read_group, "/Signal");
    free(raw_read_group);
    hid_t dset = H5Dopen(hdf5_file, signal_path, H5P_DEFAULT);
    if (dset < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open dataset '%s' to read raw signal from.\n", signal_path);
#endif
        goto cleanup2;
    }

    space = H5Dget_space(dset);
    if (space < 0) {
        fprintf(stderr, "Failed to create copy of dataspace for raw signal %s.\n", signal_path);
        goto cleanup3;
    }

    H5Sget_simple_extent_dims(space, &nsample, NULL);
    rawptr = (float*)calloc(nsample, sizeof(float));
    status = H5Dread(dset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, rawptr);

    if (status < 0) {
        free(rawptr);
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to read raw data from dataset %s.\n", signal_path.c_str());
#endif
        goto cleanup4;
    }

    // convert to pA
    rawtbl = (raw_table) { nsample, 0, nsample, rawptr };
    raw_unit = scaling.range / scaling.digitisation;
    for (size_t i = 0; i < nsample; i++) {
        rawptr[i] = (rawptr[i] + scaling.offset) * raw_unit;
    }

    free(signal_path);
    cleanup4:
    H5Sclose(space);
    cleanup3:
    H5Dclose(dset);
    cleanup2:
    return rawtbl;
}


char* fast5_get_experiment_type(hid_t hdf5_file)
{
    return fast5_get_fixed_string_attribute(hdf5_file, "/UniqueGlobalKey/context_tags", "experiment_type");
}
//
//// from scrappie
float fast5_read_float_attribute(hid_t group, const char *attribute) {
    float val = NAN;
    if (group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Invalid group passed to %s:%d.", __FILE__, __LINE__);
#endif
        return val;
    }

    hid_t attr = H5Aopen(group, attribute, H5P_DEFAULT);
    if (attr < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to open attribute '%s' for reading.", attribute);
#endif
        return val;
    }

    H5Aread(attr, H5T_NATIVE_FLOAT, &val);
    H5Aclose(attr);

    return val;
}


fast5_raw_scaling fast5_get_channel_params(hid_t hdf5_file)
{
    // from scrappie
    fast5_raw_scaling scaling = { NAN, NAN, NAN, NAN };
    const char *scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(hdf5_file, scaling_path, H5P_DEFAULT);
    if (scaling_group < 0) {
#ifdef DEBUG_FAST5_IO
        fprintf(stderr, "Failed to group %s.", scaling_path);
#endif
        return scaling;
    }

    scaling.digitisation = fast5_read_float_attribute(scaling_group, "digitisation");
    scaling.offset = fast5_read_float_attribute(scaling_group, "offset");
    scaling.range = fast5_read_float_attribute(scaling_group, "range");
    scaling.sample_rate = fast5_read_float_attribute(scaling_group, "sampling_rate");

    H5Gclose(scaling_group);

    return scaling;
}


char* fast5_get_fixed_string_attribute(hid_t hdf5_file, char* group_name, char* attribute_name)
    {
    size_t storage_size;
    char* buffer;
    hid_t group, attribute, attribute_type;
    int ret;
    char* out;

    // according to http://hdf-forum.184993.n3.nabble.com/check-if-dataset-exists-td194725.html
    // we should use H5Lexists to check for the existence of a group/dataset using an arbitrary path
    ret = H5Lexists(hdf5_file, group_name, H5P_DEFAULT);
    if(ret <= 0) {
        return "";
    }

    // Open the group /Raw/Reads/Read_nnn
    group = H5Gopen(hdf5_file, group_name, H5P_DEFAULT);
    if(group < 0) {
    #ifdef DEBUG_FAST5_IO
    fprintf(stderr, "could not open group %s\n", group_name);
    #endif
    goto close_group;
    }

    // Ensure attribute exists
    ret = H5Aexists(group, attribute_name);
    if(ret <= 0) {
    goto close_group;
    }

    // Open the attribute
    attribute = H5Aopen(group, attribute_name, H5P_DEFAULT);
    if(attribute < 0) {
    #ifdef DEBUG_FAST5_IO
    fprintf(stderr, "could not open attribute: %s\n", attribute_name.c_str());
    #endif
    goto close_attr;
    }

    // Get data type and check it is a fixed-length string
    attribute_type = H5Aget_type(attribute);
    if(H5Tis_variable_str(attribute_type)) {
    #ifdef DEBUG_FAST5_IO
    fprintf(stderr, "variable length string detected -- ignoring attribute\n");
    #endif
    goto close_type;
    }

    // Get the storage size and allocate
    storage_size = H5Aget_storage_size(attribute);
    buffer = (char*)calloc(storage_size + 1, sizeof(char));

    // finally read the attribute
    ret = H5Aread(attribute, attribute_type, buffer);
    if(ret >= 0) {
        out = buffer;
    }

    // clean up
    free(buffer);
    close_type:
            H5Tclose(attribute_type);
    close_attr:
            H5Aclose(attribute);
    close_group:
            H5Gclose(group);

    return out;
    }


void* fast5_set_event_table(hid_t hdf5_file, char* table_name, event_table *et) {
//    create empty event table with correct number of elements
//    int NFIELDS = 6;

    size_t dst_size = sizeof(event_t);
    size_t n_events = et->n;
//    event_t dst_buf[n_events];
    event_t dst_buf[n_events];
    for (int i = 0; i<n_events; i++){
        dst_buf[i].stdv = et->event[i].stdv;
        dst_buf[i].state = et->event[i].state;
        dst_buf[i].mean = et->event[i].mean;
        dst_buf[i].start = et->event[i].start;
        dst_buf[i].pos = et->event[i].pos;
        dst_buf[i].length = et->event[i].length;
    }


    size_t dst_offset[6]=  {HOFFSET(event_t, start),
                                  HOFFSET(event_t, length),
                                  HOFFSET(event_t, mean),
                                  HOFFSET(event_t, stdv),
                                  HOFFSET(event_t, pos),
                                  HOFFSET(event_t, state)};

//    size_t dst_sizes[NFIELDS] = {sizeof(dst_buf[0].start),
//                                 sizeof(dst_buf[0].length),
//                                 sizeof(dst_buf[0].mean),
//                                 sizeof(dst_buf[0].stdv),
//                                 sizeof(dst_buf[0].pos),
//                                 sizeof(dst_buf[0].state)};

    const char *field_names[6] = {"start", "length", "mean", "stdv", "pos", "state"};

    hid_t field_type[6];
    hsize_t chunk_size = 10;
    int *fill_data = NULL;
    int compress = 0;

    field_type[0] = H5T_NATIVE_UINT64;
    field_type[1] = H5T_NATIVE_FLOAT;
    field_type[2] = H5T_NATIVE_FLOAT;
    field_type[3] = H5T_NATIVE_FLOAT;
    field_type[4] = H5T_NATIVE_INT;
    field_type[5] = H5T_NATIVE_INT;


    H5TBmake_table("Table Title", hdf5_file, table_name, 6, n_events,
                   dst_size, field_names, dst_offset, field_type,
                   chunk_size, fill_data, compress, dst_buf);

}

void* fast5_set_basecall_event_table(hid_t hdf5_file, char* table_name, basecalled_event_table *et) {

    size_t n_events = et->n;
    size_t k = strlen(et->event[1].model_state);
//    printf("hey this is it %i \n", k+1);
    typedef struct {
        float start;
        float length;
        float mean;
        float stdv;
        char model_state[k+1];
        int move;
        uint64_t raw_start;
        uint64_t raw_length;
        double p_model_state;
    } _basecalled_event;

    size_t dst_size = sizeof(_basecalled_event);

    _basecalled_event dst_buf[n_events];

    for (int i = 0; i<n_events; i++){
        dst_buf[i].stdv = et->event[i].stdv;
        strcpy(dst_buf[i].model_state, et->event[i].model_state);
//        dst_buf[i].model_state =  et->event[i].model_state;
        dst_buf[i].mean = et->event[i].mean;
        dst_buf[i].start = et->event[i].start;
        dst_buf[i].move = et->event[i].move;
        dst_buf[i].length = et->event[i].length;
        dst_buf[i].raw_start = et->event[i].raw_start;
        dst_buf[i].raw_length = et->event[i].raw_length;
        dst_buf[i].p_model_state = et->event[i].p_model_state;

    }

    size_t dst_offset[9]=  {HOFFSET(_basecalled_event, start),
                            HOFFSET(_basecalled_event, length),
                            HOFFSET(_basecalled_event, mean),
                            HOFFSET(_basecalled_event, stdv),
                            HOFFSET(_basecalled_event, raw_start),
                            HOFFSET(_basecalled_event, raw_length),
                            HOFFSET(_basecalled_event, model_state),
                            HOFFSET(_basecalled_event, move),
                            HOFFSET(_basecalled_event, p_model_state)};

//    size_t dst_sizes[NFIELDS] = {sizeof(dst_buf[0].start),
//                                 sizeof(dst_buf[0].length),
//                                 sizeof(dst_buf[0].mean),
//                                 sizeof(dst_buf[0].stdv),
//                                 sizeof(dst_buf[0].pos),
//                                 sizeof(dst_buf[0].state)};

    const char *field_names[9] = {"start", "length", "mean", "stdv", "raw_start", "raw_length", "model_state",
                                  "move", "p_model_state"};

    hid_t field_type[9];
    hsize_t chunk_size = 10;
    int *fill_data = NULL;
    int compress = 0;
    hid_t string_type = H5Tcopy( H5T_C_S1 );
    H5Tset_size( string_type, k+1 );
    field_type[0] = H5T_NATIVE_FLOAT;
    field_type[1] = H5T_NATIVE_FLOAT;
    field_type[2] = H5T_NATIVE_FLOAT;
    field_type[3] = H5T_NATIVE_FLOAT;
    field_type[4] = H5T_NATIVE_UINT64;
    field_type[5] = H5T_NATIVE_UINT64;
    field_type[6] = string_type;
    field_type[7] = H5T_NATIVE_INT;
    field_type[8] = H5T_NATIVE_DOUBLE;

    H5TBmake_table("Table Title", hdf5_file, table_name, 9, n_events,
                   dst_size, field_names, dst_offset, field_type,
                   chunk_size, fill_data, compress, dst_buf);

    H5Tclose(string_type);

}


float fast5_get_start_time(hid_t hdf5_file){
    char* read_name = fast5_get_raw_read_group(hdf5_file);

    hid_t scaling_group = H5Gopen(hdf5_file, read_name, H5P_DEFAULT);
    return fast5_read_float_attribute(scaling_group, "start_time");
}

basecalled_event_table event_table_to_basecalled_table(event_table *et, fast5_raw_scaling scaling, float start_time){

    basecalled_event_table basecalled_et = { 0 };
    basecalled_et.event = calloc(et->n, sizeof(basecalled_event));
    basecalled_et.start = et->start;
    basecalled_et.end = et->end;
    basecalled_et.n = et->n;

    for (int i=0; i < et->n; i++){
        basecalled_et.event[i].raw_start = et->event[i].start;
        basecalled_et.event[i].raw_length = (uint64_t) et->event[i].length;
        basecalled_et.event[i].mean = et->event[i].mean;
        basecalled_et.event[i].stdv = et->event[i].stdv;
        basecalled_et.event[i].start = (((float) et->event[i].start) / scaling.sample_rate) + (start_time / scaling.sample_rate);
        basecalled_et.event[i].length = et->event[i].length / scaling.sample_rate;
        basecalled_et.event[i].p_model_state = 0.0;
    }
    return basecalled_et;

}


NanoporeReadAdjustmentParameters estimate_scalings_using_mom(const char* sequence, StateMachine pore_model, event_table et) {

    NanoporeReadAdjustmentParameters out;
    int64_t k = pore_model.kmerLength;
    size_t n_kmers = (strlen(sequence) - k + 1);
    char* alphabet = pore_model.alphabet;
    int64_t alphabet_size = pore_model.alphabetSize;
    double *eventModel = pore_model.EMISSION_MATCH_MATRIX;

//    size_t n_kmers = sequence.size() - k + 1;
//    const Alphabet* alphabet = pore_model.alphabet;

    // Calculate summary statistics over the events and
    // the model implied by the read
    double event_level_sum = 0.0f;
    for(int64_t i = 0; i < et.n; ++i) {
        event_level_sum += et.event[i].mean;
    }

    double kmer_level_sum = 0.0f;
    double kmer_level_sq_sum = 0.0f;
    for(int64_t i = 0; i < n_kmers; ++i) {
//        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
//        double l = pore_model.get_parameters(kmer_rank).level_mean;
        char* kmer = stString_getSubString(sequence, i, k);
        int64_t kmerIndex = kmer_id(kmer, alphabet, alphabet_size, k);
        // get the µ for level and noise for the model
        double levelMean = emissions_signal_getModelLevelMean(eventModel, kmerIndex);
        kmer_level_sum += levelMean;
        kmer_level_sq_sum += pow(levelMean, 2.0f);
        free(kmer);
    }

    double shift = event_level_sum / et.n - kmer_level_sum / n_kmers;

    // estimate scale
    double event_level_sq_sum = 0.0f;
    for(size_t i = 0; i < et.n; ++i) {
        event_level_sq_sum += pow(et.event[i].mean - shift, 2.0);
    }

    double scale = (event_level_sq_sum / et.n) / (kmer_level_sq_sum / n_kmers);

    out = set4_NanoporeReadAdjustmentParameters(shift, scale, 0.0, 1.0);

    #if DEBUG_PRINT_STATS
    fprintf(stderr, "event mean: %.2lf kmer mean: %.2lf shift: %.2lf\n", event_level_sum / et.n, kmer_level_sum / n_kmers, out.shift);
            fprintf(stderr, "event sq-mean: %.2lf kmer sq-mean: %.2lf scale: %.2lf\n", event_level_sq_sum / et.n, kmer_level_sq_sum / n_kmers, out.scale);
            fprintf(stderr, "truth shift: %.2lf scale: %.2lf\n", pore_model.shift, pore_model.scale);
    #endif
    return out;
}

void* update_SignalMachineWithNanoporeParameters(NanoporeReadAdjustmentParameters npp, StateMachine *sM){
    sM->scale = npp.scale;
    sM->shift = npp.shift;
    sM->var   = npp.var;
}



#define event_kmer_to_band(ei, ki) (((ei) + 1) + ((ki) + 1))
#define band_event_to_offset(bi, ei) band_lower_left[bi].event_idx - (ei)
#define band_kmer_to_offset(bi, ki) ((ki) - band_lower_left[bi].kmer_idx)
#define is_offset_valid(offset) (offset) >= 0 && (offset) < bandwidth
#define event_at_offset(bi, offset) band_lower_left[(bi)].event_idx - (offset)
#define kmer_at_offset(bi, offset) band_lower_left[(bi)].kmer_idx + (offset)

//#define move_down(curr_band) { (curr_band).event_idx + 1, (curr_band).kmer_idx }
//#define move_right(curr_band) { (curr_band).event_idx, (curr_band).kmer_idx + 1 }
#define new_max(x,y) ((x) >= (y)) ? (x) : (y)
#define new_min(x,y) ((x) <= (y)) ? (x) : (y)

//#define len(arr) ((int) (sizeof (arr) / sizeof (arr)[0]))
struct EventKmerPair {
    int event_idx;
    int kmer_idx;
};


struct EventKmerPair move_down(struct EventKmerPair curr_band) {
    struct EventKmerPair moved_down = {curr_band.event_idx + 1, curr_band.kmer_idx};
    return moved_down;
}

struct EventKmerPair move_right(struct EventKmerPair curr_band) {
    struct EventKmerPair moved_right = {curr_band.event_idx, curr_band.kmer_idx + 1};
    return moved_right;
}


void alignedPair_destruct(struct AlignedPair *aligned_pair) {
    free(aligned_pair);
}

struct AlignedPair *alignedPair_construct(int ref_pos, int read_pos) {
    struct AlignedPair *alignedPair = st_malloc(sizeof(struct AlignedPair));
//    kmer indx
    alignedPair->ref_pos = ref_pos;
//    event index
    alignedPair->read_pos = read_pos;

    return alignedPair;
}


stList* adaptive_banded_simple_event_align(event_table et, StateMachine *pore_model, char* sequence) {

    StateMachine3 *sM3 = (StateMachine3 *) pore_model;

//    size_t strand_idx = 0;
    int64_t k = sM3->model.kmerLength;
//    char *alphabet = sM3->model.alphabet;
    size_t n_kmers = (strlen(sequence) - k + 1);
    size_t n_events = et.n;
//    int64_t alphabet_size = sM3->model.alphabetSize;

//    const Alphabet* alphabet = pore_model.pmalphabet;
//    size_t n_events = read.events[strand_idx].size();
//    size_t n_kmers = sequence.size() - k + 1;

    // backtrack markers
    const uint8_t FROM_D = 0;
    const uint8_t FROM_U = 1;
    const uint8_t FROM_L = 2;

    // qc
    double min_average_log_emission = -5.0;
    int max_gap_threshold = 50;

    // banding
    int bandwidth = 100;
    int half_bandwidth = bandwidth / 2;

    // transition penalties
    double events_per_kmer = (double) n_events / n_kmers;
    double p_stay = 1 - (1 / (events_per_kmer + 1));

    // setting a tiny skip penalty helps keep the true alignment within the adaptive band
    // this was empirically determined
    double epsilon = 1e-10;
    double lp_skip = log(epsilon);
    double lp_stay = log(p_stay);
    double lp_step = log(1.0 - exp(lp_skip) - exp(lp_stay));
    double lp_trim = log(0.01);

    // dp matrix
    size_t n_rows = n_events + 1;
    size_t n_cols = n_kmers + 1;
    size_t n_bands = n_rows + n_cols;

    // Initialize

    // Precompute k-mer ranks to avoid doing this in the inner loop
//    std::vector<size_t> kmer_ranks(n_kmers);
    char *kmer_list[n_kmers];

    for (size_t i = 0; i < n_kmers; ++i) {
        char *kmer = stString_getSubString(sequence, i, k);
        kmer_list[i] = kmer;
//        kmer_ranks[i] = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
    }

//    typedef std::vector<float> bandscore;
//    typedef std::vector<uint8_t> bandtrace;
    double bands[n_bands][bandwidth];
    uint8_t trace[n_bands][bandwidth];
    for (size_t j = 0; j < n_bands; j++) {
        for (size_t k = 0; k < bandwidth; k++) {
            trace[j][k] = 0;
            bands[j][k] = -INFINITY;
        }
    }

//    std::vector<bandscore> bands(n_bands);
//    std::vector<bandtrace> trace(n_bands);

//    for(size_t i = 0; i < n_bands; ++i) {
//        bands[i].resize(bandwidth, -INFINITY);
//        trace[i].resize(bandwidth, 0);
//    }
//
    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two bands have their coordinates initialized, the rest are computed adaptively
//    std::vector<EventKmerPair> band_lower_left(n_bands);
//
    struct EventKmerPair band_lower_left[n_bands];

    // initialize range of first two bands
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    // band 0: score zero in the central cell
    int start_cell_offset = band_kmer_to_offset(0, -1);
    assert(is_offset_valid(start_cell_offset));
    assert(band_event_to_offset(0, -1) == start_cell_offset);
    bands[0][start_cell_offset] = 0.0f;

    // band 1: first event is trimmed
    int first_trim_offset = band_event_to_offset(1, 0);
    assert(kmer_at_offset(1, first_trim_offset) == -1);
    assert(is_offset_valid(first_trim_offset));
    bands[1][first_trim_offset] = lp_trim;
    trace[1][first_trim_offset] = FROM_U;

    int fills = 0;
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1, first_trim_offset, 0, -1, bands[1][first_trim_offset]);
#endif

    // fill in remaining bands
    for (int band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores
        double ll = bands[band_idx - 1][0];
        double ur = bands[band_idx - 1][bandwidth - 1];
        bool ll_ob = ll == -INFINITY;
        bool ur_ob = ur == -INFINITY;

        bool right = false;
        if (ll_ob && ur_ob) {
            right = band_idx % 2 == 1;
        } else {
            right = ll < ur; // Suzuki's rule
        }

        if (right) {
            band_lower_left[band_idx] = move_right(band_lower_left[band_idx - 1]);
        } else {
            band_lower_left[band_idx] = move_down(band_lower_left[band_idx - 1]);
        }

        /*
                float max_score = -INFINITY;
                int tmp_max_offset = 0;
                for(int tmp = 0; tmp < bandwidth; ++tmp) {
                    float s = bands[band_idx - 1][tmp];
                    if(s > max_score) {
                        max_score = s;
                        tmp_max_offset = tmp;
                    }
                }
                fprintf(stderr, "bi: %d ll: %.2f up: %.2f [%d %d] [%d %d] max: %.2f [%d %d] move: %s\n",
                    band_idx, bands[band_idx - 1][0], bands[band_idx - 1][bandwidth - 1],
                    band_lower_left[band_idx - 1].event_idx, band_lower_left[band_idx - 1].kmer_idx,
                    event_at_offset(band_idx - 1, bandwidth - 1), kmer_at_offset(band_idx - 1, bandwidth - 1),
                    max_score, event_at_offset(band_idx - 1, tmp_max_offset), kmer_at_offset(band_idx - 1, tmp_max_offset),
                    (right ? "RIGHT" : "DOWN"));
        */

        // If the trim state is within the band, fill it in here
        int trim_offset = band_kmer_to_offset(band_idx, -1);
        if (is_offset_valid(trim_offset)) {
            int event_idx = event_at_offset(band_idx, trim_offset);
            if (event_idx >= 0 && event_idx < n_events) {
                bands[band_idx][trim_offset] = (lp_trim * (event_idx + 1));
                trace[band_idx][trim_offset] = FROM_U;
            } else {
                bands[band_idx][trim_offset] = -INFINITY;
            }
        }

        // Get the offsets for the first and last event and kmer
        // We restrict the inner loop to only these values
        int kmer_min_offset = band_kmer_to_offset(band_idx, 0);
        int kmer_max_offset = band_kmer_to_offset(band_idx, n_kmers);
        int event_min_offset = band_event_to_offset(band_idx, n_events - 1);
        int event_max_offset = band_event_to_offset(band_idx, -1);

        int min_offset = new_max(kmer_min_offset, event_min_offset);
        min_offset = new_max(min_offset, 0);

        int max_offset = new_min(kmer_max_offset, event_max_offset);
        max_offset = new_min(max_offset, bandwidth);

        for (int offset = min_offset; offset < max_offset; ++offset) {
            int event_idx = event_at_offset(band_idx, offset);
            int kmer_idx = kmer_at_offset(band_idx, offset);

            char *kmer = kmer_list[kmer_idx];
            double y[2] = {(double) et.event[event_idx].mean, (double) et.event[event_idx].stdv};

            int offset_up = band_event_to_offset(band_idx - 1, event_idx - 1);
            int offset_left = band_kmer_to_offset(band_idx - 1, kmer_idx - 1);
            int offset_diag = band_kmer_to_offset(band_idx - 2, kmer_idx - 1);

#ifdef DEBUG_ADAPTIVE
            // verify loop conditions
                        assert(kmer_idx >= 0 && kmer_idx < n_kmers);
                        assert(event_idx >= 0 && event_idx < n_events);
                        assert(offset_diag == band_event_to_offset(band_idx - 2, event_idx - 1));
                        assert(offset_up - offset_left == 1);
                        assert(offset >= 0 && offset < bandwidth);
#endif

            float up = is_offset_valid(offset_up) ? bands[band_idx - 1][offset_up] : -INFINITY;
            float left = is_offset_valid(offset_left) ? bands[band_idx - 1][offset_left] : -INFINITY;
            float diag = is_offset_valid(offset_diag) ? bands[band_idx - 2][offset_diag] : -INFINITY;

            double lp_emission = sM3->getMatchProbFcn(pore_model, kmer, y, TRUE);

//            float lp_emission = log_probability_match_r9(read, pore_model, kmer_rank, event_idx, strand_idx);
            float score_d = (float) (diag + lp_step + lp_emission);
            float score_u = (float) (up + lp_stay + lp_emission);
            float score_l = (float) (left + lp_skip);

            float max_score = score_d;
            uint8_t from = FROM_D;

            max_score = score_u > max_score ? score_u : max_score;
            from = max_score == score_u ? FROM_U : from;
            max_score = score_l > max_score ? score_l : max_score;
            from = max_score == score_l ? FROM_L : from;

#ifdef DEBUG_ADAPTIVE
            fprintf(stderr, "[adafill] offset-up: %d offset-diag: %d offset-left: %d\n", offset_up, offset_diag, offset_left);
                        fprintf(stderr, "[adafill] up: %.2lf diag: %.2lf left: %.2lf\n", up, diag, left);
                        fprintf(stderr, "[adafill] bi: %d o: %d e: %d k: %d s: %.2lf f: %d emit: %.2lf\n", band_idx, offset, event_idx, kmer_idx, max_score, from, lp_emission);
#endif
            bands[band_idx][offset] = max_score;
            trace[band_idx][offset] = from;
            fills += 1;
        }
    }

    /*
// Debug, print some of the score matrix
for(int col = 0; col <= 10; ++col) {
    for(int row = 0; row < 100; ++row) {
        int kmer_idx = col - 1;
        int event_idx = row - 1;
        int band_idx = event_kmer_to_band(event_idx, kmer_idx);
        int offset = band_kmer_to_offset(band_idx, kmer_idx);
        assert(offset == band_event_to_offset(band_idx, event_idx));
        assert(event_idx == event_at_offset(band_idx, offset));
        fprintf(stdout, "ei: %d ki: %d bi: %d o: %d s: %.2f\n", event_idx, kmer_idx, band_idx, offset, bands[band_idx][offset]);
    }
}
*/

//     Backtrack to compute alignment

    double sum_emission = 0;
    double n_aligned_events = 0;

    float max_score = -INFINITY;
    int curr_event_idx = 0;
    int curr_kmer_idx = (int) (n_kmers - 1);

    // Find best score between an event and the last k-mer. after trimming the remaining evnets
    for (int event_idx = 0; event_idx < n_events; ++event_idx) {
        int band_idx = event_kmer_to_band(event_idx, curr_kmer_idx);
//        assert(band_idx < bands.size());
        assert(band_idx < n_bands);

        int offset = band_event_to_offset(band_idx, event_idx);
        if (is_offset_valid(offset)) {
            float s = (float) (bands[band_idx][offset] + (n_events - event_idx) * lp_trim);
            if (s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, max_score);
#endif
//    std::vector<AlignedPair> out;
    stList *out = stList_construct3(0, (void (*)(void *)) alignedPair_destruct);
//    stList_append(AlignedPair_list, alignedPair_construct(1, 2));
//    stList_append(AlignedPair_list, alignedPair_construct(3, 4));

    int out_index = 0;
    int curr_gap = 0;
    int max_gap = 0;
    while (curr_kmer_idx >= 0 && curr_event_idx >= 0) {
        // emit alignment
        stList_append(out, alignedPair_construct(curr_kmer_idx, curr_event_idx));
#ifdef DEBUG_ADAPTIVE
        fprintf(stderr, "[adaback] ei: %d ki: %d\n", curr_event_idx, curr_kmer_idx);
#endif
        out_index += 1;
        // qc stats
//        int event_idx = event_at_offset(band_idx, offset);
//        int kmer_idx = kmer_at_offset(band_idx, offset);

        char *kmer = kmer_list[curr_kmer_idx];
        double y[2] = {(double) et.event[curr_event_idx].mean, (double) et.event[curr_event_idx].stdv};
//        size_t kmer_rank = alphabet->kmer_rank(sequence.substr(curr_kmer_idx, k).c_str(), k);
//        sum_emission += log_probability_match_r9(read, pore_model, kmer_rank, curr_event_idx, strand_idx);
        sum_emission += sM3->getMatchProbFcn(pore_model, kmer, y, TRUE);

        n_aligned_events += 1;

        int band_idx = event_kmer_to_band(curr_event_idx, curr_kmer_idx);
        int offset = band_event_to_offset(band_idx, curr_event_idx);
        assert(band_kmer_to_offset(band_idx, curr_kmer_idx) == offset);

        uint8_t from = trace[band_idx][offset];
        if (from == FROM_D) {
            curr_kmer_idx -= 1;
            curr_event_idx -= 1;
            curr_gap = 0;
        } else if (from == FROM_U) {
            curr_event_idx -= 1;
            curr_gap = 0;
        } else {
            curr_kmer_idx -= 1;
            curr_gap += 1;
            max_gap = new_max(curr_gap, max_gap);
//            max_gap = std::max(curr_gap, max_gap);
        }
    }

//    std::reverse(out.begin(), out.end());
    stList_reverse(out);
    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
//    bool spanned = out.front().ref_pos == 0 && out.back().ref_pos == n_kmers - 1;
    struct AlignedPair *front = (struct AlignedPair *) stList_get(out, 0);
    struct AlignedPair *back = (struct AlignedPair *) stList_get(out, stList_length(out)-1);

    bool spanned = front->ref_pos == 0 && back->ref_pos == n_kmers - 1;

    bool failed = false;
    if (avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
        failed = true;
        stList_destruct(out);
        stList *out = stList_construct3(0, (void (*)(void *)) alignedPair_destruct);
    }

    fprintf(stderr, "ada\t%s\t%.2lf\t%zu\t%.2lf\t%d\t%d\t%d\n", failed ? "FAILED" : "OK", events_per_kmer, strlen(sequence) , avg_log_emission, curr_event_idx, max_gap, fills);
    return out;
}

// embed event table from fast5 path, template model path and nucleotide sequence
void* load_from_raw(char* fast5_file_path, char* templateModelFile, char* sequence, char* path_to_embed) {
//    load model
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling_MeanOnly,
                                                  stateMachine3_loadTransitionsFromFile, NULL);


//    load fast5
    hid_t hdf5_file = fast5_open(fast5_file_path);
    // Get start time and event detection defaults
    float start_time = fast5_get_start_time(hdf5_file);
    const detector_param *ed_params = &event_detection_defaults;
    char *experiment_type = fast5_get_experiment_type(hdf5_file);
    char *RNA = "rna";
// if experiment is RNA, change defaults and nucleotide sequence
    if (strcmp(experiment_type, RNA) == 0) {
        ed_params = &event_detection_rna;
        char *new_sequence = stString_ReverseString(stString_replace(sequence, "U", "T"));
        sequence = new_sequence;
    }
//    get channel parameters and scale raw ADC counts to get pA raw current
    fast5_raw_scaling channel_params = fast5_get_channel_params(hdf5_file);
    raw_table rt = fast5_get_raw_samples(hdf5_file, channel_params);
    // trim using scrappie's internal method
    // parameters taken directly from scrappie defaults
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
//  get events
    event_table et = detect_events(rt, *ed_params);
    assert(rt.n > 0);
    assert(et.n > 0);
// estimate scalings via method of moments using scrappie methods
    NanoporeReadAdjustmentParameters scalings_template = estimate_scalings_using_mom(sequence, *sM, et);
//    update signal machine with new parameters
    update_SignalMachineWithNanoporeParameters(scalings_template, sM);
// banded kmer to event alignment
    stList *event_alignment = adaptive_banded_simple_event_align(et, sM, sequence);
// create new event table with our own data structure
    basecalled_event_table b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    basecalled_event_table *event_table = alignment_to_base_event_map(event_alignment, &b_et, sequence, sM);
    fast5_set_basecall_event_table(hdf5_file, path_to_embed, event_table);
    fast5_close(hdf5_file);
    stateMachine_destruct(sM);
    return 0;
}

basecalled_event_table* alignment_to_base_event_map(stList *event_alignment, basecalled_event_table* b_et,
                                                   char *sequence, StateMachine *pore_model) {

    StateMachine3 *sM3 = (StateMachine3 *) pore_model;
    int k = (int) sM3->model.kmerLength;
    // transform alignment into the base-to-event map
    int64_t alignment_length = stList_length(event_alignment);
    assert(alignment_length > 0) ;
    int k_idx;
    int event_idx;
    char *kmer;
    double lp_emission;
    int prev_event_idx = -1;
    int prev_kmer_indx = 0;

//      loop through the alignment and create assignments

    for (int64_t i = 0; i < alignment_length; ++i) {

        struct AlignedPair *aligned_pair = stList_get(event_alignment, i);
        k_idx = aligned_pair->ref_pos;
        event_idx = aligned_pair->read_pos;

        kmer = stString_getSubString(sequence, k_idx, k);
        double y[2] = {b_et->event[event_idx].mean, b_et->event[event_idx].stdv};
        lp_emission = sM3->getMatchProbFcn(pore_model, kmer, y, TRUE);
        if (event_idx == prev_event_idx) {
            if (k_idx == prev_kmer_indx){
                fprintf(stderr, "[Event Aligner] There was an error in the event map.");
            } else {
                b_et->event[event_idx].p_model_state = exp(lp_emission);
                b_et->event[event_idx].model_state = kmer;
                b_et->event[event_idx].move += k_idx-prev_kmer_indx;
//                    printf("We assigned to %i\n", event_idx);
                prev_kmer_indx = k_idx;
                prev_event_idx = event_idx;
            }
        } else {
            b_et->event[event_idx].p_model_state = exp(lp_emission);
            b_et->event[event_idx].model_state = kmer;

            if (k_idx == prev_kmer_indx){
                b_et->event[event_idx].move = 0;
            } else {
                b_et->event[event_idx].move = k_idx-prev_kmer_indx;
            }
            prev_kmer_indx = k_idx;
            prev_event_idx = event_idx;
        }
    }
    return b_et;
}




