//
// Created by Andrew Bailey on 6/8/18.
//



#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <hdf5.h>
#include <scrappie_structures.h>
#include <eventAligner.h>

#include "signalMachineUtils.h"

#define RAW_ROOT "/Raw/Reads/"
//#define DEBUG_ADAPTIVE 1
//#define DEBUG_FAST5_IO 1
//#define DEBUG_PRINT_STATS 1

hid_t fast5_open(char* filename)
    {
    hid_t hdf5file = H5Fopen(filename, H5F_ACC_RDWR, H5P_DEFAULT);
    return hdf5file;
    }


herr_t fast5_close(hid_t hdf5_file)
{
    return H5Fclose(hdf5_file);
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
    char* out = "";

    // according to http://hdf-forum.184993.n3.nabble.com/check-if-dataset-exists-td194725.html
    // we should use H5Lexists to check for the existence of a group/dataset using an arbitrary path
    ret = H5Lexists(hdf5_file, group_name, H5P_DEFAULT);
    if(ret <= 0) {
        return "";
    }

    // Open thoute group /Raw/Reads/Read_nnn
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
        out = stString_copy(buffer);
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


void fast5_basecall_event_type(hid_t* types) {

//    hid_t string_type = H5Tcopy( H5T_C_S1 );
//    H5Tset_size( string_type, MAX_KMER_SIZE + 1 );


    hid_t string_type = H5Tcopy( H5T_C_S1 );
    H5Tset_size( string_type, MAX_KMER_SIZE + 1 );
    hid_t bce_tid = H5Tcreate (H5T_COMPOUND, sizeof(basecalled_event));
    H5Tinsert(bce_tid, "start", HOFFSET(basecalled_event, start), H5T_NATIVE_DOUBLE);
    H5Tinsert(bce_tid, "length", HOFFSET(basecalled_event, length), H5T_NATIVE_DOUBLE);
    H5Tinsert(bce_tid, "mean", HOFFSET(basecalled_event, mean), H5T_NATIVE_DOUBLE);
    H5Tinsert(bce_tid, "stdv", HOFFSET(basecalled_event, stdv), H5T_NATIVE_DOUBLE);
    H5Tinsert(bce_tid, "raw_start", HOFFSET(basecalled_event, raw_start), H5T_NATIVE_UINT64);
    H5Tinsert(bce_tid, "raw_length", HOFFSET(basecalled_event, raw_length), H5T_NATIVE_UINT64);
    H5Tinsert(bce_tid, "model_state", HOFFSET(basecalled_event, model_state), string_type);
    H5Tinsert(bce_tid, "move", HOFFSET(basecalled_event, move), H5T_NATIVE_INT);
    H5Tinsert(bce_tid, "p_model_state", HOFFSET(basecalled_event, p_model_state), H5T_NATIVE_DOUBLE);

    types[0] = bce_tid;
    types[1] = string_type;
}


herr_t fast5_set_basecall_event_table(hid_t hdf5_file, char* table_location, basecalled_event_table *et) {

    /* prep */
    size_t n_events = 0;
    size_t k = 0;
    for (int i = 0; i < et->n; i++) {
        if (strlen(et->event[i].model_state) != 0) {
            if (k == 0) k = strlen(et->event[i].model_state);
            n_events++;
        }
    }
    assert(k != 0);
    assert(n_events != 0);
    basecalled_event *dst_buf = calloc(n_events, sizeof(basecalled_event));
    hid_t      dataset, space; /* Handles */
    herr_t     status;
    hsize_t    dim = n_events;   /* Dataspace dimensions */
    hid_t      dtypes[2];
    fast5_basecall_event_type(dtypes);
    hid_t      bce_tid = dtypes[0];
    hid_t      str_tid = dtypes[1];

    /* Initialize the data */
    int j = 0;
    for (int i = 0; i<et->n; i++){
        // skip unaligned events (these were trimmed)
        if (strlen(et->event[i].model_state) == 0) continue;

        dst_buf[j].stdv = et->event[i].stdv;
        strcpy(dst_buf[j].model_state, et->event[i].model_state);
        dst_buf[j].mean = et->event[i].mean;
        dst_buf[j].start = et->event[i].start;
        dst_buf[j].move = et->event[i].move;
        dst_buf[j].length = et->event[i].length;
        dst_buf[j].raw_start = et->event[i].raw_start;
        dst_buf[j].raw_length = et->event[i].raw_length;
        dst_buf[j].p_model_state = et->event[i].p_model_state;
        j++;
    }
    assert(j == n_events);

    /* Create and write the dataset. */
    space = H5Screate_simple(1, &dim, NULL);
    dataset = H5Dcreate2(hdf5_file, table_location, bce_tid, space, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset, bce_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst_buf);

    /* Release resources */
    H5Tclose(str_tid);
    H5Tclose(bce_tid);
    H5Sclose(space);
    H5Dclose(dataset);
    free(dst_buf);

    return status;
}

// There is no size checking in H5Dread, so dst_buf MUST BE of the right size.
// I recommend that this only be used in unit tests
herr_t fast5_get_basecall_events(hid_t hdf5_file, char* table_location, basecalled_event *dst_buf) {
    hid_t      dataset;
    herr_t     status;
    hid_t      dtypes[2];
    fast5_basecall_event_type(dtypes);
    hid_t      bce_tid = dtypes[0];
    hid_t      str_tid = dtypes[1];

    dataset = H5Dopen2(hdf5_file, table_location, H5P_DEFAULT);

    status = H5Dread(dataset, bce_tid, H5S_ALL, H5S_ALL, H5P_DEFAULT, dst_buf);

    H5Tclose(str_tid);
    H5Tclose(bce_tid);
    H5Dclose(dataset);

    return status;
}


float fast5_get_start_time(hid_t hdf5_file){
    char* read_name = fast5_get_raw_read_group(hdf5_file);

    hid_t scaling_group = H5Gopen(hdf5_file, read_name, H5P_DEFAULT);
    return fast5_read_float_attribute(scaling_group, "start_time");
}

basecalled_event_table* event_table_to_basecalled_table(event_table *et, fast5_raw_scaling scaling, float start_time){

    basecalled_event_table* basecalled_et = malloc(sizeof(basecalled_event_table));
    basecalled_et->event = calloc(et->n, sizeof(basecalled_event));
    basecalled_et->start = et->start;
    basecalled_et->end = et->end;
    basecalled_et->n = et->n;
    basecalled_et->aln_n = 0;

    for (int i=0; i < et->n; i++){
        basecalled_et->event[i].raw_start = et->event[i].start;
        basecalled_et->event[i].raw_length = (uint64_t) et->event[i].length;
        basecalled_et->event[i].mean = et->event[i].mean;
        basecalled_et->event[i].stdv = et->event[i].stdv;
        basecalled_et->event[i].start = (((double) et->event[i].start) / scaling.sample_rate) + (start_time / scaling.sample_rate);
        basecalled_et->event[i].length = et->event[i].length / scaling.sample_rate;
        basecalled_et->event[i].p_model_state = 0.0;
        basecalled_et->event[i].move = 0;
        basecalled_et->event[i].model_state[0] = '\0';

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
#define band_array_access_2d(bi, offset) (bi) * bandwidth + (offset)


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
    int64_t kmer_length = sM3->model.kmerLength;
//    char *alphabet = sM3->model.alphabet;
    size_t n_kmers = (strlen(sequence) - kmer_length + 1);
    size_t n_events = et.n;
//    int64_t alphabet_size = sM3->model.alphabetSize;

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
        char *kmer = stString_getSubString(sequence, i, kmer_length);
        kmer_list[i] = kmer;
//        kmer_ranks[i] = alphabet->kmer_rank(sequence.substr(i, k).c_str(), k);
    }

//    typedef std::vector<float> bandscore;
//    typedef std::vector<uint8_t> bandtrace;
    //double bands[n_bands][bandwidth];
    //uint8_t trace[n_bands][bandwidth];
    double* bands = (double*) malloc(sizeof(double) * n_bands * bandwidth);
    uint8_t* trace = (uint8_t*) malloc(sizeof(uint8_t) * n_bands * bandwidth);
    for (size_t j = 0; j < n_bands; j++) {
        for (size_t k = 0; k < bandwidth; k++) {
            trace[band_array_access_2d(j,k)] = 0;
            bands[band_array_access_2d(j,k)] = -INFINITY;
        }
    }

    // Keep track of the event/kmer index for the lower left corner of the band
    // these indices are updated at every iteration to perform the adaptive banding
    // Only the first two bands have their coordinates initialized, the rest are computed adaptively
    struct EventKmerPair band_lower_left[n_bands];

    // initialize range of first two bands
    band_lower_left[0].event_idx = half_bandwidth - 1;
    band_lower_left[0].kmer_idx = -1 - half_bandwidth;
    band_lower_left[1] = move_down(band_lower_left[0]);

    // band 0: score zero in the central cell
    int start_cell_offset = band_kmer_to_offset(0, -1);
    assert(is_offset_valid(start_cell_offset));
    assert(band_event_to_offset(0, -1) == start_cell_offset);
    bands[band_array_access_2d(0, start_cell_offset)] = 0.0f;

    // band 1: first event is trimmed
    int first_trim_offset = band_event_to_offset(1, 0);
    assert(kmer_at_offset(1, first_trim_offset) == -1);
    assert(is_offset_valid(first_trim_offset));
    bands[band_array_access_2d(1, first_trim_offset)] = lp_trim;
    trace[band_array_access_2d(1, first_trim_offset)] = FROM_U;

    int fills = 0;
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[trim] bi: %d o: %d e: %d k: %d s: %.2lf\n", 1, first_trim_offset, 0, -1, bands[band_array_access_2d(1,first_trim_offset)]);
#endif

    // fill in remaining bands
    for (int band_idx = 2; band_idx < n_bands; ++band_idx) {
        // Determine placement of this band according to Suzuki's adaptive algorithm
        // When both ll and ur are out-of-band (ob) we alternate movements
        // otherwise we decide based on scores

        double ll = bands[band_array_access_2d(band_idx - 1, 0)];
        double ur = bands[band_array_access_2d(band_idx - 1, bandwidth - 1)];
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
                bands[band_array_access_2d(band_idx, trim_offset)] = (lp_trim * (event_idx + 1));
                trace[band_array_access_2d(band_idx, trim_offset)] = FROM_U;
            } else {
                bands[band_array_access_2d(band_idx, trim_offset)] = -INFINITY;
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

            float up = is_offset_valid(offset_up) ? bands[band_array_access_2d(band_idx - 1, offset_up)] : -INFINITY;
            float left = is_offset_valid(offset_left) ? bands[band_array_access_2d(band_idx - 1, offset_left)] : -INFINITY;
            float diag = is_offset_valid(offset_diag) ? bands[band_array_access_2d(band_idx - 2, offset_diag)] : -INFINITY;

            double lp_emission = sM3->getMatchProbFcn(pore_model, kmer, y, TRUE);

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
            bands[band_array_access_2d(band_idx, offset)] = max_score;
            trace[band_array_access_2d(band_idx, offset)] = from;
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
        assert(band_idx < n_bands);

        int offset = band_event_to_offset(band_idx, event_idx);
        if (is_offset_valid(offset)) {
            float s = (float) (bands[band_array_access_2d(band_idx, offset)] + (n_events - event_idx) * lp_trim);
            if (s > max_score) {
                max_score = s;
                curr_event_idx = event_idx;
            }
        }
    }
#ifdef DEBUG_ADAPTIVE
    fprintf(stderr, "[adaback] ei: %d ki: %d s: %.2f\n", curr_event_idx, curr_kmer_idx, max_score);
#endif
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

        uint8_t from = trace[band_array_access_2d(band_idx, offset)];
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
        }
    }

    stList_reverse(out);
    // QC results
    double avg_log_emission = sum_emission / n_aligned_events;
    struct AlignedPair *front = (struct AlignedPair *) stList_get(out, 0);
    struct AlignedPair *back = (struct AlignedPair *) stList_get(out, stList_length(out)-1);

    bool spanned = front->ref_pos == 0 && back->ref_pos == n_kmers - 1;

    bool failed = false;
    if (avg_log_emission < min_average_log_emission || !spanned || max_gap > max_gap_threshold) {
        failed = true;
        stList_destruct(out);
        stList *out = stList_construct3(0, (void (*)(void *)) alignedPair_destruct);
    }

    fprintf(stderr, "%s\tevents_per_kmer:%.2lf\tsequence_len:%zu\tavg_log_emission:%.2lf\tcurr_event_idx:%d\tmax_gap:%d\tfills:%d\n", failed ? "FAILED" : "OK", events_per_kmer, strlen(sequence) , avg_log_emission, curr_event_idx, max_gap, fills);
    return out;
}

// embed event table from fast5 path, template model path and nucleotide sequence
herr_t load_from_raw(char* fast5_file_path, char* templateModelFile, char* sequence, char* path_to_embed) {
    // prep
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling_MeanOnly,
                                                  stateMachine3_loadTransitionsFromFile, NULL);
    hid_t hdf5_file = fast5_open(fast5_file_path);
    float start_time = fast5_get_start_time(hdf5_file);
    const detector_param *ed_params = &event_detection_defaults;

    // if experiment is RNA, change defaults and nucleotide sequence
    char *experiment_type = fast5_get_experiment_type(hdf5_file);
    if (strcmp(experiment_type, "rna") == 0) {
        ed_params = &event_detection_rna;
        char *new_sequence = stString_ReverseString(stString_replace(sequence, "U", "T"));
        sequence = new_sequence;
    }

    // get channel parameters and scale raw ADC counts to get pA raw current
    fast5_raw_scaling channel_params = fast5_get_channel_params(hdf5_file);
    raw_table rt = fast5_get_raw_samples(hdf5_file, channel_params);

    // trim using scrappie's internal method
    // parameters taken directly from scrappie defaults
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);

    // get events
    event_table et = detect_events(rt, *ed_params);
    assert(rt.n > 0);
    assert(et.n > 0);

    // estimate scalings via method of moments using scrappie methods
    NanoporeReadAdjustmentParameters scalings_template = estimate_scalings_using_mom(sequence, *sM, et);
    // update signal machine with new parameters
    update_SignalMachineWithNanoporeParameters(scalings_template, sM);

    // banded kmer to event alignment
    stList *event_alignment = adaptive_banded_simple_event_align(et, sM, sequence);

    // create new event table with our own data structure
    basecalled_event_table* b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    alignment_to_base_event_map(event_alignment, b_et, sequence, sM);
    herr_t write_success = fast5_set_basecall_event_table(hdf5_file, path_to_embed, b_et);

    // cleanup
    free(b_et->event);
    free(b_et);
    herr_t close_success = fast5_close(hdf5_file);
    stateMachine_destruct(sM);
    return write_success | close_success;
}

void alignment_to_base_event_map(stList *event_alignment, basecalled_event_table* b_et,
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
    b_et->aln_n = 0;

    // loop through the alignment and create assignments
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
                strcpy(b_et->event[event_idx].model_state, kmer);
                b_et->event[event_idx].move += k_idx-prev_kmer_indx;
                prev_kmer_indx = k_idx;
                prev_event_idx = event_idx;
            }
        } else {
            b_et->event[event_idx].p_model_state = exp(lp_emission);
            strcpy(b_et->event[event_idx].model_state, kmer);

            if (k_idx == prev_kmer_indx){
                b_et->event[event_idx].move = 0;
            } else {
                b_et->event[event_idx].move = k_idx-prev_kmer_indx;
            }
            prev_kmer_indx = k_idx;
            prev_event_idx = event_idx;
        }
        b_et->aln_n += 1;
        free(kmer);
    }
}
