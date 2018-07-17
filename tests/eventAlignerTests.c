//
// Created by Andrew Bailey on 6/25/18.
//

#include <stdlib.h>
#include "hdf5.h"

#include <math.h>
#include <inttypes.h>
#include <unistd.h>
#include <eventAligner.h>
#include <scrappie_structures.h>
#include "CuTest.h"
#include "eventAligner.h"
#include "sonLibString.h"
#include "sonLibList.h"
#include "signalMachine.h"
#include "stateMachine.h"

//#define HOME "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/"
#define HOME "../" //this is based on where travis runs tests from

static void test_fast5_get_raw_read_name(CuTest *testCase) {
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    char* answer = "Read_61";
    hid_t fast5_handle = fast5_open(path);
    char* read_name = fast5_get_raw_read_name(fast5_handle);
    fast5_close(fast5_handle);
    CuAssertStrEquals(testCase, answer, read_name);
    free(path);
}

static void test_fast5_get_raw_read_group(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);

    char* group = fast5_get_raw_read_group(fast5_handle);
    char* answer = "/Raw/Reads/Read_61";
    CuAssertStrEquals(testCase, answer, group);
    free(group);
    fast5_close(fast5_handle);

}

static void test_fast5_get_fixed_string_attribute(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    char* raw_read_group = "/Raw/Reads/Read_61";

    char* attribute = fast5_get_fixed_string_attribute(fast5_handle, raw_read_group, "read_id");
    char* answer = "8898d755-e46d-4cbe-842c-545a97718b9d";

    CuAssertStrEquals(testCase, answer, attribute);
    fast5_close(fast5_handle);

}

static void test_fast5_get_read_id(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    char* attribute = fast5_get_read_id(fast5_handle);
    char* answer = "8898d755-e46d-4cbe-842c-545a97718b9d";
    CuAssertStrEquals(testCase, answer, attribute);
    fast5_close(fast5_handle);

}

static void test_fast5_get_experiment_type(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    char* attribute = fast5_get_experiment_type(fast5_handle);
    char* answer = "rna";
    CuAssertStrEquals(testCase, answer, attribute);
    fast5_close(fast5_handle);

}

static void test_fast5_read_float_attribute(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    const char *scaling_path = "/UniqueGlobalKey/channel_id";

    hid_t scaling_group = H5Gopen(fast5_handle, scaling_path, H5P_DEFAULT);
    float digitisation = fast5_read_float_attribute(scaling_group, "digitisation");

    float d_answer = 8192.0;
    CuAssertDblEquals(testCase, d_answer, digitisation, .001);
    fast5_close(fast5_handle);

}

static void test_fast5_get_channel_params(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    fast5_raw_scaling rna_scaling_params = fast5_get_channel_params(fast5_handle);
    double digitisation_answer = 8192.0;
    double range_answer = 1182.140015;
    double sample_rate_answer = 3012.0;
    double offset_answer = 3.0;

    CuAssertDblEquals(testCase, digitisation_answer, rna_scaling_params.digitisation, 0.001);
    CuAssertDblEquals(testCase, range_answer, rna_scaling_params.range, 0.001);
    CuAssertDblEquals(testCase, sample_rate_answer, rna_scaling_params.sample_rate, 0.001);
    CuAssertDblEquals(testCase, offset_answer, rna_scaling_params.offset, 0.001);
    fast5_close(fast5_handle);


}

static void test_fast5_get_raw_samples(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    fast5_raw_scaling rna_scaling_params = fast5_get_channel_params(fast5_handle);
    raw_table rna_table = fast5_get_raw_samples(fast5_handle, rna_scaling_params);
    size_t start_answer = 0;
    size_t end_answer = 25794;
    size_t n_answer = 25794;

    CuAssertDblEquals(testCase, (double) 93.797729, (double) rna_table.raw[0], 0.0001);
    CuAssertDblEquals(testCase, (double) 94.230644, (double) rna_table.raw[1], 0.0001);
    CuAssertIntEquals(testCase, (unsigned int) start_answer, (unsigned int) rna_table.start);
    CuAssertIntEquals(testCase, (unsigned int) end_answer, (unsigned int) rna_table.end);
    CuAssertIntEquals(testCase, (unsigned int) n_answer, (unsigned int) rna_table.n);
    fast5_close(fast5_handle);


}

//TODO do we need "set_event_table" without "setbasecalled_event_table"?
//static void test_fast5_set_event_table(CuTest *testCase){
//    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
//    hid_t fast5_handle = fast5_open(path);
//    fast5_raw_scaling channel_params = fast5_get_channel_params(fast5_handle);
//    raw_table rt = fast5_get_raw_samples(fast5_handle, channel_params);
//    const detector_param* ed_params = &event_detection_rna;
//
//    event_table et = detect_events(rt, *ed_params);
//    fast5_set_event_table(fast5_handle, "Analyses/UnitTest_Events", &et);
//    size_t n_events = et.n;
//    event_t dst_buf[n_events];
//    size_t dst_offset[6]=  {HOFFSET(event_t, start),
//                            HOFFSET(event_t, length),
//                            HOFFSET(event_t, mean),
//                            HOFFSET(event_t, stdv),
//                            HOFFSET(event_t, pos),
//                            HOFFSET(event_t, state)};
//
//    size_t dst_sizes[6] = {sizeof(dst_buf[0].start),
//                           sizeof(dst_buf[0].length),
//                           sizeof(dst_buf[0].mean),
//                           sizeof(dst_buf[0].stdv),
//                           sizeof(dst_buf[0].pos),
//                           sizeof(dst_buf[0].state)};
//    size_t dst_size = sizeof(event_t);
//
//    //TODO
//    //H5TBread_table( fast5_handle, "Analyses/UnitTest_Events", dst_size, dst_offset, dst_sizes, dst_buf );
//
//    CuAssertIntEquals(testCase, 0, (int) dst_buf[0].start);
//    CuAssertDblEquals(testCase, 7.000000, dst_buf[0].length, 0.001);
//    CuAssertDblEquals(testCase, 92.086693, dst_buf[0].mean, 0.001);
//    CuAssertDblEquals(testCase, 3.655048, dst_buf[0].stdv, 0.001);
//    CuAssertIntEquals(testCase, -1, dst_buf[0].pos);
//    CuAssertIntEquals(testCase, -1, dst_buf[0].state);
//
//    CuAssertIntEquals(testCase, 7, (int) dst_buf[1].start);
//    CuAssertDblEquals(testCase, 15.000000, dst_buf[1].length, 0.001);
//    CuAssertDblEquals(testCase, 87.082771, dst_buf[1].mean, 0.001);
//    CuAssertDblEquals(testCase, 1.637721, dst_buf[1].stdv, 0.001);
//    CuAssertIntEquals(testCase, -1, dst_buf[1].pos);
//    CuAssertIntEquals(testCase, -1, dst_buf[1].state);
//    //TODO
//    //H5TBdelete_record (fast5_handle, "Analyses/SignalAlign_Basecall_1D_000/Events", 0, n_events);
//    H5Ldelete( fast5_handle, "Analyses/UnitTest_Events", H5P_DEFAULT );
//
//    fast5_close(fast5_handle);
//
//}

static void test_fast5_get_start_time(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    float start_time = fast5_get_start_time(fast5_handle);
    CuAssertIntEquals(testCase, 232505, (int) start_time);
}

static void test_event_table_to_basecalled_table(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    float start_time = fast5_get_start_time(fast5_handle);

    fast5_raw_scaling channel_params = fast5_get_channel_params(fast5_handle);
    raw_table rt = fast5_get_raw_samples(fast5_handle, channel_params);
    const detector_param* ed_params = &event_detection_rna;

    event_table et = detect_events(rt, *ed_params);
    basecalled_event_table b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    CuAssertIntEquals(testCase, 7, (int) b_et.event[1].raw_start);
    CuAssertIntEquals(testCase, 15, (int)  b_et.event[1].raw_length);
    CuAssertDblEquals(testCase, 87.082771, b_et.event[1].mean, 0.001);
    CuAssertDblEquals(testCase, 1.637721, b_et.event[1].stdv, 0.001);
    CuAssertDblEquals(testCase, 77.195221, b_et.event[1].start, 0.0001);
    CuAssertDblEquals(testCase, 0.004980, b_et.event[1].length, 0.0001);
    CuAssertDblEquals(testCase, 0.0, b_et.event[1].p_model_state, 0.0001);

}

static void test_fast5_set_basecall_event_table(CuTest *testCase){

    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);

    basecalled_event_table basecalled_et = { 0 };
    basecalled_et.event = calloc(2, sizeof(basecalled_event));
    basecalled_et.start = 0;
    basecalled_et.end = 1;
    basecalled_et.n = 2;

    basecalled_et.event[0].raw_start = 0;
    basecalled_et.event[0].raw_length = 1;
    basecalled_et.event[0].mean = 100.0;
    basecalled_et.event[0].stdv = 10.0;
    basecalled_et.event[0].start = 0.1;
    basecalled_et.event[0].length = 0.1;
    basecalled_et.event[0].p_model_state = 0.99;
    basecalled_et.event[0].model_state = stString_copy("GATTA");
    basecalled_et.event[0].move = 0;

    basecalled_et.event[1].raw_start = 1;
    basecalled_et.event[1].raw_length = 2;
    basecalled_et.event[1].mean = 200.0;
    basecalled_et.event[1].stdv = 20.0;
    basecalled_et.event[1].start = 0.2;
    basecalled_et.event[1].length = 0.2;
    basecalled_et.event[1].p_model_state = 0.5;
    basecalled_et.event[1].model_state = stString_copy("TTACA");
    basecalled_et.event[1].move = 1;

    basecalled_event_table* bet_ptr = &basecalled_et;

    herr_t status = fast5_set_basecall_event_table(fast5_handle, "Analyses/UnitTest_Events", bet_ptr);
    fast5_close(fast5_handle);
    CuAssertIntEquals(testCase, 0, (int) status);

    //open and read
    fast5_handle = fast5_open(path);
    basecalled_event dst_buf[2];
    status = fast5_get_basecall_events(fast5_handle, "Analyses/UnitTest_Events", dst_buf);
    CuAssertIntEquals(testCase, 0, (int) status);

    CuAssertIntEquals(testCase, 0, (int) dst_buf[0].raw_start);
    CuAssertIntEquals(testCase, 1, (int) dst_buf[0].raw_length);
    CuAssertDblEquals(testCase, 100.0, dst_buf[0].mean, 0.00001);
    CuAssertDblEquals(testCase, 10.0, dst_buf[0].stdv, 0.00001);
    CuAssertDblEquals(testCase, 0.1, dst_buf[0].start, 0.00001);
    CuAssertDblEquals(testCase, 0.1, dst_buf[0].length, 0.00001);
    CuAssertDblEquals(testCase, 0.99, dst_buf[0].p_model_state, 0.00001);
    CuAssertStrEquals(testCase, "GATTA", dst_buf[0].model_state);
    CuAssertIntEquals(testCase, 0, dst_buf[0].move);

    CuAssertIntEquals(testCase, 1, (int) dst_buf[1].raw_start);
    CuAssertIntEquals(testCase, 2, (int) dst_buf[1].raw_length);
    CuAssertDblEquals(testCase, 200.0, dst_buf[1].mean, 0.00001);
    CuAssertDblEquals(testCase, 20.0, dst_buf[1].stdv, 0.00001);
    CuAssertDblEquals(testCase, 0.2, dst_buf[1].start, 0.00001);
    CuAssertDblEquals(testCase, 0.2, dst_buf[1].length, 0.00001);
    CuAssertDblEquals(testCase, 0.5, dst_buf[1].p_model_state, 0.00001);
    CuAssertStrEquals(testCase, "TTACA", dst_buf[1].model_state);
    CuAssertIntEquals(testCase, 1, dst_buf[1].move);

    H5Ldelete( fast5_handle, "Analyses/UnitTest_Events", H5P_DEFAULT );

    fast5_close(fast5_handle);

}



static void test_alignment_to_base_event_map(CuTest *testCase){
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling_MeanOnly,
                                                  stateMachine3_loadTransitionsFromFile, NULL);

//    fprintf(stderr, "%s\n",sM->alphabet);
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    fast5_raw_scaling channel_params = fast5_get_channel_params(fast5_handle);
    raw_table rt = fast5_get_raw_samples(fast5_handle, channel_params);
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    const detector_param* ed_params = &event_detection_rna;

    const char* sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA";
//    fprintf(stderr, "%s\n", sequence);
    char* new_sequence = stString_ReverseString(stString_replace(sequence, "U", "T"));
//    fprintf(stderr, "%s\n", new_sequence);
    event_table et = detect_events(rt, *ed_params);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(new_sequence, *sM, et);
    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    stList* event_alignment = adaptive_banded_simple_event_align(et, sM, new_sequence);
    float start_time = fast5_get_start_time(fast5_handle);
    basecalled_event_table b_et = event_table_to_basecalled_table(&et, channel_params, start_time);

    basecalled_event_table* bet_ptr = alignment_to_base_event_map(event_alignment, &b_et, new_sequence, sM);
    char* prev_kmer = bet_ptr->event[0].model_state;
    int move;
    int k = (int) sM->kmerLength;
    for (int i = 1; i < bet_ptr->n; i++){
//        printf("%s: %f\n", bet_ptr->event[i].model_state, bet_ptr->event[i].p_model_state);
//        printf("%i\n", bet_ptr->event[i].move);
        move = bet_ptr->event[i].move;
        CuAssertStrEquals(testCase, stString_getSubString(prev_kmer, move, k-move),
                          stString_getSubString(bet_ptr->event[i].model_state, 0, k-move));
        prev_kmer = bet_ptr->event[i].model_state;
    }
    fast5_close(fast5_handle);

}

static void test_estimate_scalings_using_mom(CuTest *testCase){
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling,
                                                  stateMachine3_loadTransitionsFromFile, NULL);
//    fprintf(stderr, "%s\n",sM->alphabet);
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    fast5_raw_scaling channel_params = fast5_get_channel_params(fast5_handle);
    raw_table rt = fast5_get_raw_samples(fast5_handle, channel_params);
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    const detector_param* ed_params = &event_detection_rna;

    const char* sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA";
//    fprintf(stderr, "%s\n", sequence);
    char* new_sequence = stString_ReverseString(stString_replace(sequence, "U", "T"));
//    fprintf(stderr, "%s\n", new_sequence);
    event_table et = detect_events(rt, *ed_params);
    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(new_sequence, *sM, et);
    CuAssertDblEquals(testCase, 1.016111, test_squiggle_scalings.scale, 0.0001);
    CuAssertDblEquals(testCase, 20.720264, test_squiggle_scalings.shift, 0.0001);
    free(new_sequence);
    fast5_close(fast5_handle);
}

static void test_set_NanoporeReadAdjustmentParameters(CuTest *testCase){
    double scale = 2.0;
    double shift = 2.0;
    double drift = 2.0;
    double var = 2.0;
    double scale_sd = 2.0;
    double var_sd = 2.0;
    double shift_sd = 2.0;

    NanoporeReadAdjustmentParameters set_6_test = set7_NanoporeReadAdjustmentParameters(scale, shift, drift, var, scale_sd, var_sd, shift_sd);
    CuAssertDblEquals(testCase, scale, set_6_test.scale, 0.0001);
    CuAssertDblEquals(testCase, shift, set_6_test.shift, 0.0001);
    CuAssertDblEquals(testCase, drift, set_6_test.drift, 0.0001);
    CuAssertDblEquals(testCase, var, set_6_test.var, 0.0001);
    CuAssertDblEquals(testCase, scale_sd, set_6_test.scale_sd, 0.0001);
    CuAssertDblEquals(testCase, var_sd, set_6_test.var_sd, 0.0001);
    CuAssertDblEquals(testCase, shift_sd, set_6_test.shift_sd, 0.0001);

    double scale_sd_4 = 1.0;
    double var_sd_4 = 1.0;
    double shift_sd_4 = 0.0;

    NanoporeReadAdjustmentParameters set_4_test = set4_NanoporeReadAdjustmentParameters(scale, shift, drift, var);
    CuAssertDblEquals(testCase, scale, set_4_test.scale, 0.0001);
    CuAssertDblEquals(testCase, shift, set_4_test.shift, 0.0001);
    CuAssertDblEquals(testCase, drift, set_4_test.drift, 0.0001);
    CuAssertDblEquals(testCase, var, set_4_test.var, 0.0001);
    CuAssertDblEquals(testCase, scale_sd_4, set_4_test.scale_sd, 0.0001);
    CuAssertDblEquals(testCase, var_sd_4, set_4_test.var_sd, 0.0001);
    CuAssertDblEquals(testCase, shift_sd_4, set_4_test.shift_sd, 0.0001);


}

static void test_sonlib_lists(CuTest *testCase){
    stList *AlignedPair_list = stList_construct3(0, (void (*)(void *)) alignedPair_destruct);
    stList_append(AlignedPair_list, alignedPair_construct(1, 2));
    stList_append(AlignedPair_list, alignedPair_construct(3, 4));
    struct AlignedPair *something = stList_get(AlignedPair_list, 1);
    CuAssertIntEquals(testCase, 3, something->ref_pos);
    CuAssertIntEquals(testCase, 4, something->read_pos);
    stList_destruct(AlignedPair_list);
}

static void test_update_SignalMachineWithNanoporeParameters(CuTest *testCase){
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling,
                                                  stateMachine3_loadTransitionsFromFile, NULL);
    NanoporeReadAdjustmentParameters test_squiggle_scalings = set4_NanoporeReadAdjustmentParameters(3.0, 4.0, 5.0, 6.0);

    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    CuAssertDblEquals(testCase, 4.0, sM->scale, 0.0001);
    CuAssertDblEquals(testCase, 6.0, sM->var, 0.0001);
    CuAssertDblEquals(testCase, 3.0, sM->shift, 0.0001);

}

static void test_adaptive_banded_simple_event_align(CuTest *testCase){
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling_MeanOnly,
                                                  stateMachine3_loadTransitionsFromFile, NULL);

    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    fast5_raw_scaling channel_params = fast5_get_channel_params(fast5_handle);
    raw_table rt = fast5_get_raw_samples(fast5_handle, channel_params);
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    const detector_param* ed_params = &event_detection_rna;

    const char* sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA";
    char* new_sequence = stString_ReverseString(stString_replace(sequence, "U", "T"));
    event_table et = detect_events(rt, *ed_params);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(new_sequence, *sM, et);
    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    stList* something = adaptive_banded_simple_event_align(et, sM, new_sequence);
    struct AlignedPair *test_end = (struct AlignedPair*) stList_pop(something);
    CuAssertIntEquals(testCase, 453, test_end->ref_pos);
    CuAssertIntEquals(testCase, 1219, test_end->read_pos);
    fast5_close(fast5_handle);
}


static void test_load_from_raw_rna(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    char* sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA";
    herr_t status = load_from_raw(path, templateModelFile, sequence, "Analyses/SignalAlign_Basecall_1D_000/Events");
    CuAssertIntEquals(testCase, 0, (int) status);

    hid_t fast5_handle = fast5_open(path);

    size_t n_events = 1220;

    basecalled_event dst_buf[n_events];
    size_t dst_offset[9]=  {HOFFSET(basecalled_event, start),
                            HOFFSET(basecalled_event, length),
                            HOFFSET(basecalled_event, mean),
                            HOFFSET(basecalled_event, stdv),
                            HOFFSET(basecalled_event, raw_start),
                            HOFFSET(basecalled_event, raw_length),
                            HOFFSET(basecalled_event, model_state),
                            HOFFSET(basecalled_event, move),
                            HOFFSET(basecalled_event, p_model_state)};

    size_t dst_sizes[9] = {sizeof(dst_buf[0].start),
                           sizeof(dst_buf[0].length),
                           sizeof(dst_buf[0].mean),
                           sizeof(dst_buf[0].stdv),
                           sizeof(dst_buf[0].raw_start),
                           sizeof(dst_buf[0].raw_length),
                           sizeof(dst_buf[0].model_state),
                           sizeof(dst_buf[0].move),
                           sizeof(dst_buf[0].p_model_state)};

    size_t dst_size = sizeof(basecalled_event);

    //TODO
    //H5TBread_table( fast5_handle, "Analyses/SignalAlign_Basecall_1D_000/Events", dst_size, dst_offset, dst_sizes, dst_buf );

    CuAssertIntEquals(testCase, 0, (int) dst_buf[0].raw_start);
    CuAssertDblEquals(testCase, 7.000000, dst_buf[0].raw_length, 0.001);
    CuAssertDblEquals(testCase, 92.086693, dst_buf[0].mean, 0.001);
    CuAssertDblEquals(testCase, 3.655048, dst_buf[0].stdv, 0.001);
    CuAssertStrEquals(testCase, "AACCT", dst_buf[0].model_state);
    CuAssertIntEquals(testCase, 0, dst_buf[0].move);

    CuAssertIntEquals(testCase, 7, (int) dst_buf[1].raw_start);
    CuAssertDblEquals(testCase, 15.000000, dst_buf[1].raw_length, 0.001);
    CuAssertDblEquals(testCase, 87.082771, dst_buf[1].mean, 0.001);
    CuAssertDblEquals(testCase, 1.637721, dst_buf[1].stdv, 0.001);
    CuAssertStrEquals(testCase, "ACCTA", dst_buf[1].model_state);
    CuAssertIntEquals(testCase, 1, dst_buf[1].move);

    //TODO
    H5Ldelete( fast5_handle, "Analyses/SignalAlign_Basecall_1D_000/Events", H5P_DEFAULT );
}


static void test_load_from_raw_dna(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5");
    char* templateModelFile = stString_concat(HOME, "models/testModelR9p4_5mer_acegt_template.model");
    char* sequence = "TGCATGCCGTTTCCGTTACGTATTGCTAATCACGGGTGACGCCGTTTTCGCGCATTGAGCGAATCAGCAAAACCATCGCTAAACCACGGCTAACCCGGCGATGTGTGCTCCGTTCTCATCGACATCCCAAACCGTCAAACCATCCGGCGACAATCAGATCAGCGCAAAGATAATTAACCCACGTTGCAGGTAAATGCCACTTTGCGGATCGCGTTCGCCGTAGCCAGACGTAGCCCATCAGCACGCCACCACGCCAAGCCCGCCAAACCGGCCGCTGAATTTTGCTGCACATAGCCGCTTAACAGGGCGCTGATGGCGTAATGACGAGGTGGCTTACCGCTACCGAGGCGTTTTCCACCGCACCGCCGAGATACCACCACCAAGAGCAGGTTAAAGGGATCAGCGAAATTGCATTAACGCGTGGGTGGTAACGCCAGAACTCAAATTTCGTGTTGGATCCGAATGGCGGGGCCAGCCATAACGTCACTTCCTGATCGCCGGGGTACATGGCAATAAACACCACCACGCAGGCGATCATCATCACCCGGTTACCGGACCTGCGCGTTCACACCAAGGCGGCAAAGGATAACGGCGATAATGCAGGCCCTGCCGACCATGGCCTGCCTGCCAGCTCGCCGCCAGATAACACGCGCGGATCTGCCGGTTTCGAAAACGCGCCAGCTCCGCCCGTACGCTATCGGCACAGGACTCATCCGCCAGCCAGACATCGCTTTGGTTATGTTGTTGGTCGTGGGGATAATACCCTGCGTCGCCATGTAATCAACAAACGCCTGCGCCACGCGGGGATTGTAAAAGTGTCATCAGCATCGTTGCTGTCATATTCCACAAGGGACAGTATAAAGCGTTACGCGCCGTACGCCACCTCTGCGGAAACTGACGTTGCCGGGCTTCAAGCCGCCGTCAATGCTATGAACCACATCGTAGCCCTGTTGCAGCAGATGCTGCGCCGCGCCTTTGCTGCTGTGCCGTGATAACACATCACCATCACCGGGTGTCAAAGTCGTTATCACGCATAAAAGCGCCCAGCGTGTCGTTGGTTAAATGGAAAGCCTGCACCACTGTCCATTGCGAAACTCTGTGGATCCGCCTTGAATATCAGAGCCAGCACCGCCTCTTTCCTGCAACTTCTGGTGCGCGTCGGCAGCGTTAATGCATTGAACTGATCCATGCGTCTCTCTTTCTTTGACAAGTGGGCAGAATTACCGCACAGTTTACGTCGAAGCGGCAGATAAACGCCATAATGTTGCCATATCATAAAATGTTTTCAATGTTACCCCGCGATTCTTTGCTAATATGTTCGATAACGAACATTTATGAGCTTTAACGAAAGTGAATGAGGGCGGCATGACCAAAGATCTGATTGTGATGGGGGCGGCATCAATGGTGCTGGTATCGCGGCAGACGCCGCTGGACGCGGACCTCCGTGCTGATGCTGGGCGCAGGATCTCGCTTGCGGCCTCTTCCGCCAGTTCAAACTCATTCGGTGGCCTGCATACCTTGAGCACTATAATTCCGCTTTGGTCAGCGAGGCGCTGAACGTGAAGTGCTGCTGAAAATGGCCCCGCATATCGCGCCTTCCCGATGCGTTTCGCCTGCCATCGTCCGCATCTGCGCCCGGCGTGGATGATTCGCATTGGTCTGTTTATGTACGATCATCTGGGTAAACGCACCAGCTTGCCGGGATCAACTGGTTTGCGTTTTGGCGCAAATTCAGTGTTGAAATTAGCGCGGATTCCAGATATTCTGACTGTTGGTGGGCGACGCCCGTCTGGTACTCGCCAACGCCCAGATGGTGTGGTGCGTAAAGGCGGCGAAGTGCTACTCGGACTCGCGCCACCTCTGCTCGCCGCGAAACGGCCTGTGGATTGTGGAAGCGGAAGATCGATACCGGCAAAAATATAGCTGGCAAGCGCGCGGCAGTTAACGCGCCACCGGCCCGTGGGGTGAAACAGTTCTTCGACGACGGGATGCATCTGCCTTCGCCTTATGGCATTCGCCTGATCAAGGCAGCCATATTGTGGTGCCGCGCGTGCATGCAAGCAAGCCTACATTCTGCAAAACGAGATAAACGTATTGTGTTCGTGATCCCGTGGATGGACGAGTTTTCCATCATCGGCACTACCGATGTCGAGTACAAAGGCGAATCGAAAGCGGTGAAGATTGAAGAGTGAAATCAATTACCTGCTGAATAACACGCACTTTAAAAAGCCAGTTAAGCCATTGACGATATCGTCTGGACCTACTCCGGTGTGCGTCCGCTGTGTGATGATAGGTCGACTCGCCGCAGGCTATTACCCGTGATTACACCCTTGATATTCGCGGTGAAATGGCAAGCACCGCTGCTGTCGGTATTCGGCGGTAAGCTGACCACCTACCGAAAACTGGCGGAACATGCGCTGGAAAACTAACGCCGTATTATCAGGGTATTGGCCCGGCATGGACGAAAGAGAGTGTGCTACCGGGTGGCGCCATTGAAGCGACCGCGACGTTAATACCGCTCGCCTGCGCCGCCGCTATCCGTTCCTGACTGAATCGCTGGCGCGTCGTACGCTCGCACTTACGGCAGCAACAGCGAGCTGCTGCTCGGCAATGCAGGAACGGTAAGCGATCTCAGGGAAGATTTCGGTCATGGGTTCTACAAGCGGAGCTGAAATACGGTGGATGGGTCCACCGCCGACGACGCCCTGTAGTCGCACAAGTAAGGCATGTGGCTAAATGCGGATCAACAATCTCGTGTGAGTCGGTGGCTGGTGGAGTATACGCAGCAGGTTATCATGGCGTCGTAAATTAACGTAAGGTGATCGGTCAGATTTCGACTGGCCTGAGACTGATGACAAACCTACAAAACTGCCTGATGCGCTTCGCTTATCAGGCCTACGTAGTTTATGCAATATATTGAATTTGCATGGTCTTGTAGGCCAGATAAGACGTTCACGTCGCATCCGGCATGAACTCAGCACTTTGTCAAAAATCTAACCTACTTTTAATTCAGGGAATTACCGCAAAGCCCACATACCATCATGCAACGTAACAAAACTCAGGCACGTTCCCCTCGCCCCGAGAAAATAGCATTAATGCGCCCAGCGCCAGCATAAAAATTTTGAGCGGTGGTGTTGGCGTGATAATACAAACTAATAATACCGGCAAGTCCGACACCCAGCATGTAACCACCGCCAAAATTGCGCCAGTATGGGGATGCCGAAAAGTCATTACAACGAGGTCAAAATCCATTTCTGTTTTGCATTATTCTTCCATTCTTTTGAATGGTGAAGTGCTCCCGTGTTATACTTATGGACACTTTTCGAAATGATGGCGGAAAAACGGGACCGCTGGCCCCGTTCTGGCTGACCGGTGAACTTACAATCTCACCGGATCGATATGCCAGATATGATCGGCGTACTCTTTGATGGTACGGTCAGAAGAAGTAGCCCATATTGGCAATGTTCAACATCGCTTTGCGGTCCCACTCTTCCTGAAGCTCGTAGGGACATCGACTTTATCTGACAATCGACATAGCTGCGATAATCCATGAGTACGTATTGATCGCCGCCAGTTGATCAGCAGTCAACAGGTCGCGATAGCGACCGGATCTTCCGGACCCACCGCTGCCGATTTGCGTCGGCACCTGGTGCAGCTCCTCATCTTTCTCGTAGTATTCACGCGGTTGTAGCCCTGACGACGCGCAGTTCTTCCACTTCCGCTGTGTTACCAAAATAAGATATTGTCAGCGCCGACATGATCAGCATCTCGTATTCACCTCAACGCGATAGTCAGCGCACCGTTAAGCGCAAACTTCACGTATTACTGGTGCCGGAAGCTGCGTCCCTGCCAGCGAAATCTGTTCGAACAGATCTGCCGCCGGAATGATCAGCTGCGCCGGTTCTTGTAGTTCAGGATGAACACGACTTTCAACATCGCCAATCTGCGGATCGTTGTTGATCACTTTCGCTACGTCATTGATCAAATGAATAATGTGCTTCGCCATGTAATAGGCCGAAACCGCCTTACCGCCAAAATATTCACGCGCGGCCCACTTCGCATCGGTCGGCCTTGATGCGGTTGCGTAATCACATGCAACACATTCATCAATTGACGTTTGTATTCGTGAATACGTTTGGTTGTACATCGAACAACGCCTTTGGATTCACCACCACATTCAGCTGCTGGGCGATATACTCTGCCAGACGCTTTTGTTCTCCAGCCATGATGCACAGCGTGATTAACCATTGGGAAATCACAGTGTTGTTTTTGCAACTCATTAAGCAGGCTAAGGTCGGTGCCAGTTGCGGCCAGGTGTTCGTCCAGCACGGCTGAAGCGATGGGTTCGCTACCGCCAGCCAGCGACGCGGCGTACACCGTTGGTGACGTTGGTGAAACTGACCGGGAAGATTTCGCAAAGTCGGCAAACAACGTTGCACCATCAGTTAGAGTGCAGTTCCGATACACCGTTAACTGTGGCTCACAACAACCGCCAGCCGGGCCATACGCACGACGACCGTTGGATTCATCAATGATCGACGCCCGTCCCAGCGCAGATCGGTATCGTTCAGTTCTACTGTTCCTGCAAGTTTCAGAAATAGTCGTTGATTTCAAGATGATCTGCAGGTGACGCGGCAAATTTTACCAGCATATCAACCGGCGGGTTTCCAGCGCCTCGCTCATCAGCGTGGTTAGTGTAGGAAGACCTGACAACACACCTCAAACGCGTCGTCCAGCTAAATTGGTGCTCATCGATCAGCAGACGCATCCTCTCAGGAATCGACAATTGCGGATGGGTATCATTGAGATGATCGCGATTTTATCCGCCAGGTTATCGTAGGTTTTATGCAACTGATAATGGCGGCAAATGTCCCTGAATGGTCGAAACCGGGAAGTATTCTGACGCAGGCACAGCTCACGCCAGGTAAATTGAGTCATCCGGATACAACCGCGAGATACGTTCTCGAGTGGTTTTATCTTCTTCCACTGCCGCGAAGTAGTCACCTGGTTGAATTACCGAGTTAATTTCGCTAAACGCGCACTCCACAAACGCAGCGTGTTGGTCGCACGTCGGTGTCGTAACCGAGGATTATCGTAAGCGACTCCCAGAATCTCTTCGGTTTTCAATCCAGCGCGTTTTACCTTCCTGCTGAATGCACGTACGCCAAAACGGACTTTATATGGCGCGTGTTGTGGCGTTTGAATTCCACAGGTTACCGTATTCCAGCCAGTAGTCTGGCGACTCTTTCTGGCTACCGTTAACGATGTTCTGCTTGAACATACCCATGGTCATAGCGGATGCCGTAACCACGCCCAGCAGCCCTAACGTCGCCAGAATCAAAGGAAGCAAGCCGCCAGACGTCCCAGGCACCGTGCGAGGCACGGATCGACACCTTCATCAATCAGCTCTTCAGTTAACCCATCGCTTCCAGTGCGCCCTGTACATCTTCGTAAATTCTAGCGACAACATGGCGTTGGAGGCTTGTTAATCAAAACTCCATCGACTGGTATTAAACCTGATGAGTTTCTTACGACAACTGGGCACGGTTTGAACGTAACCAGCGCTCCGGACGATCGCGCACAGCAAATAACGTTGCGTTGTGATCATGTTGTGGCGACGACCGGGTTCCTTTCAATCGTAAACATCAGCTTGTAAGCGATAGAGTGTAGGCTTCTACGCGCTAAGCGTGGGCGATAGATATGTAAACGGAGCGATATATAAACCACAAACTAACACAAGCGATAGTAGCTCACGGTACGACTTCGCCGCGACCTGCCAGCTAAAATCCATTGCCATAGCCTGACGTTTGCACAAACCGCCACAGTGAAGGACGGGACCACAGTACAAAAGCACGTCCGAATAGCCCGTAACAGCGACCGGGCATTACTATCTTCAAGACAAAGCCACAGCGACGCCATCTGCAAGGTTCTCGAGAACGATCCAGAAACCGTGTCATAAACCCACCGGTGCACATAACGGCAACGTACCATGCTTCAATCCATAAAGTTGCGTTAAGCCGCACGGTTCAAGCGGCTGGGCACCAGGGCCAGCGTCCGCGCCGCCCTTAATGCGATGCGAAAATGCTTCGTGACCATCAGGCGCCCACCCTGACCGGGGTATTCCGCTGCCGCCGCAAAACCTTCCTGCAACGCACCGGATCGCCCGCGCCGAGGTAGCATAACGCCCTGCTCCAGAAGACCCGGTAGGCTTCCAGCACCGGGTCAGACGCTCTGGCTGGTCGGGCGGCTCACCACCGCAAAAGCGGCACTTTATCGTCAACGCCGCCCCATTGCGATTTGTAACTGGCGCTTATTTTCCGCTTTATCTTCAACGTATCGCGGGTGTATGAGGCCAACAGTAAGTCGATCTCTGGGCTCCGTTTTCTCCCACGCCGTTCAGTACGCCGGAAGACGCCCTTCACGGTGACAGCACGTAGACCTTCGCGCCACCGTGAGCAAACTGCGGTCGGTGATCTCGCGAGCGTGGTTGGACGACCGCCGTAATGCTTGATCGGCATGGTACAGACCGGCCTTCAGAAACT";
    herr_t status = load_from_raw(path, templateModelFile, sequence, "Analyses/SignalAlign_Basecall_1D_000/Events");
    CuAssertIntEquals(testCase, 0, (int) status);

    hid_t fast5_handle = fast5_open(path);

    size_t n_events = 1220;
    basecalled_event dst_buf[n_events];
    size_t dst_offset[9]=  {HOFFSET(basecalled_event, start),
                            HOFFSET(basecalled_event, length),
                            HOFFSET(basecalled_event, mean),
                            HOFFSET(basecalled_event, stdv),
                            HOFFSET(basecalled_event, raw_start),
                            HOFFSET(basecalled_event, raw_length),
                            HOFFSET(basecalled_event, model_state),
                            HOFFSET(basecalled_event, move),
                            HOFFSET(basecalled_event, p_model_state)};

    size_t dst_sizes[9] = {sizeof(dst_buf[0].start),
                           sizeof(dst_buf[0].length),
                           sizeof(dst_buf[0].mean),
                           sizeof(dst_buf[0].stdv),
                           sizeof(dst_buf[0].raw_start),
                           sizeof(dst_buf[0].raw_length),
                           sizeof(dst_buf[0].model_state),
                           sizeof(dst_buf[0].move),
                           sizeof(dst_buf[0].p_model_state)};

    size_t dst_size = sizeof(basecalled_event);

    //TODO
    //H5TBread_table( fast5_handle, "Analyses/SignalAlign_Basecall_1D_000/Events", dst_size, dst_offset, dst_sizes, dst_buf );

    CuAssertIntEquals(testCase, 0, (int) dst_buf[0].raw_start);
    CuAssertDblEquals(testCase, 7.000000, dst_buf[0].raw_length, 0.001);
    CuAssertDblEquals(testCase, 92.086693, dst_buf[0].mean, 0.001);
    CuAssertDblEquals(testCase, 3.655048, dst_buf[0].stdv, 0.001);
    CuAssertStrEquals(testCase, "AACCT", dst_buf[0].model_state);
    CuAssertIntEquals(testCase, 0, dst_buf[0].move);

    CuAssertIntEquals(testCase, 7, (int) dst_buf[1].raw_start);
    CuAssertDblEquals(testCase, 15.000000, dst_buf[1].raw_length, 0.001);
    CuAssertDblEquals(testCase, 87.082771, dst_buf[1].mean, 0.001);
    CuAssertDblEquals(testCase, 1.637721, dst_buf[1].stdv, 0.001);
    CuAssertStrEquals(testCase, "ACCTA", dst_buf[1].model_state);
    CuAssertIntEquals(testCase, 1, dst_buf[1].move);

    //TODO
    //H5TBdelete_record (fast5_handle, "Analyses/SignalAlign_Basecall_1D_000/Events", 0, n_events);
    H5Ldelete( fast5_handle, "Analyses/SignalAlign_Basecall_1D_000/Events", H5P_DEFAULT );
}



CuSuite *eventAlignerTestSuite(void) {
    st_system("pwd");
    printf("\ntest suite location\n");

    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_fast5_get_raw_read_name);
    SUITE_ADD_TEST(suite, test_fast5_get_raw_read_group);
    SUITE_ADD_TEST(suite, test_fast5_get_fixed_string_attribute);
    SUITE_ADD_TEST(suite, test_fast5_get_read_id);
    SUITE_ADD_TEST(suite, test_fast5_get_experiment_type);
    SUITE_ADD_TEST(suite, test_fast5_read_float_attribute);
    SUITE_ADD_TEST(suite, test_fast5_get_channel_params);
    SUITE_ADD_TEST(suite, test_fast5_get_raw_samples);
    //SUITE_ADD_TEST(suite, test_fast5_set_event_table);
    SUITE_ADD_TEST(suite, test_fast5_get_start_time);
    SUITE_ADD_TEST(suite, test_fast5_set_basecall_event_table);

    SUITE_ADD_TEST(suite, test_set_NanoporeReadAdjustmentParameters);
    SUITE_ADD_TEST(suite, test_estimate_scalings_using_mom);
    SUITE_ADD_TEST(suite, test_update_SignalMachineWithNanoporeParameters);
    SUITE_ADD_TEST(suite, test_sonlib_lists);
    SUITE_ADD_TEST(suite, test_adaptive_banded_simple_event_align);
    SUITE_ADD_TEST(suite, test_event_table_to_basecalled_table);

    SUITE_ADD_TEST(suite, test_alignment_to_base_event_map);
    //TODO
//    SUITE_ADD_TEST(suite, test_load_from_raw_rna);
//    SUITE_ADD_TEST(suite, test_load_from_raw_dna);

    return suite;
}

//int main(int argc, char *argv[]) {
//    // collect output and create a new test suite
//    CuString *output = CuStringNew();
//    CuSuite *suite = CuSuiteNew();
//    // add and run this test suite
//    CuSuiteAddSuite(suite, eventAlignerTestSuite());
//    CuSuiteRun(suite);
//    CuSuiteSummary(suite, output);
//    CuSuiteDetails(suite, output);
//    printf("%s\n", output->buffer);
//    CuStringDelete(output);
//    int good = suite->failCount > 0;
//    return good;
//}