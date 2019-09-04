//
// Created by Andrew Bailey on 6/25/18.
//

#include <stdlib.h>
#include "hdf5.h"

#include <eventAligner.h>
#include "CuTest.h"
#include "signalMachine.h"
#include "kseq.h"
#include <zlib.h>


#define HOME "../" //this is based on where travis runs tests from
//#define HOME "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/" //this is based on where travis runs tests from

#define EVENT_LOCATION "/Analyses/UnittestEvents"

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
    basecalled_event_table* b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    CuAssertIntEquals(testCase, 7, (int) b_et->event[1].raw_start);
    CuAssertIntEquals(testCase, 15, (int)  b_et->event[1].raw_length);
    CuAssertDblEquals(testCase, 87.082771, b_et->event[1].mean, 0.001);
    CuAssertDblEquals(testCase, 1.637721, b_et->event[1].stdv, 0.001);
    CuAssertDblEquals(testCase, 77.195221, b_et->event[1].start, 0.0001);
    CuAssertDblEquals(testCase, 0.004980, b_et->event[1].length, 0.0001);
    CuAssertDblEquals(testCase, 0.0, b_et->event[1].p_model_state, 0.0001);

    free(b_et->event);
    free(b_et);
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
    strcpy(basecalled_et.event[0].model_state, "GATTA");
    basecalled_et.event[0].move = 0;

    basecalled_et.event[1].raw_start = 1;
    basecalled_et.event[1].raw_length = 2;
    basecalled_et.event[1].mean = 200.0;
    basecalled_et.event[1].stdv = 20.0;
    basecalled_et.event[1].start = 0.2;
    basecalled_et.event[1].length = 0.2;
    basecalled_et.event[1].p_model_state = 0.5;
    strcpy(basecalled_et.event[1].model_state, "TTACA");
    basecalled_et.event[1].move = 1;

    basecalled_event_table* bet_ptr = &basecalled_et;

    herr_t status = fast5_set_basecall_event_table(fast5_handle, EVENT_LOCATION, bet_ptr);
    fast5_close(fast5_handle);
    CuAssertIntEquals(testCase, 0, (int) status);

    //open and read
    fast5_handle = fast5_open(path);
    basecalled_event dst_buf[2];
    status = fast5_get_basecall_events(fast5_handle, EVENT_LOCATION, dst_buf);
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

    H5Ldelete( fast5_handle, EVENT_LOCATION, H5P_DEFAULT );

    fast5_close(fast5_handle);

}

static void test_alignment_to_base_event_map(CuTest *testCase){
    char* templateModelFile = stString_concat(HOME, "models/testModelR9p4_5mer_acegt_template.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling_MeanOnly,
                                                  stateMachine3_loadTransitionsFromFile, NULL);

    char* path = stString_concat(HOME, "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    fast5_raw_scaling channel_params = fast5_get_channel_params(fast5_handle);
    raw_table rt = fast5_get_raw_samples(fast5_handle, channel_params);
    int trim_start = 200;
    int trim_end = 10;
    int varseg_chunk = 100;
    float varseg_thresh = 0.0;
    trim_and_segment_raw(rt, trim_start, trim_end, varseg_chunk, varseg_thresh);
    const detector_param *ed_params = &event_detection_defaults;

    char* sequence = "TGCATGCCGTTTCCGTTACGTATTGCTAATCACGGGTGACGCCGTTTTCGCGCATTGAGCGAATCAGCAAAACCATCGCTAAACCACGGCTAACCCGGCGATGTGTGCTCCGTTCTCATCGACATCCCAAACCGTCAAACCATCCGGCGACAATCAGATCAGCGCAAAGATAATTAACCCACGTTGCAGGTAAATGCCACTTTGCGGATCGCGTTCGCCGTAGCCAGACGTAGCCCATCAGCACGCCACCACGCCAAGCCCGCCAAACCGGCCGCTGAATTTTGCTGCACATAGCCGCTTAACAGGGCGCTGATGGCGTAATGACGAGGTGGCTTACCGCTACCGAGGCGTTTTCCACCGCACCGCCGAGATACCACCACCAAGAGCAGGTTAAAGGGATCAGCGAAATTGCATTAACGCGTGGGTGGTAACGCCAGAACTCAAATTTCGTGTTGGATCCGAATGGCGGGGCCAGCCATAACGTCACTTCCTGATCGCCGGGGTACATGGCAATAAACACCACCACGCAGGCGATCATCATCACCCGGTTACCGGACCTGCGCGTTCACACCAAGGCGGCAAAGGATAACGGCGATAATGCAGGCCCTGCCGACCATGGCCTGCCTGCCAGCTCGCCGCCAGATAACACGCGCGGATCTGCCGGTTTCGAAAACGCGCCAGCTCCGCCCGTACGCTATCGGCACAGGACTCATCCGCCAGCCAGACATCGCTTTGGTTATGTTGTTGGTCGTGGGGATAATACCCTGCGTCGCCATGTAATCAACAAACGCCTGCGCCACGCGGGGATTGTAAAAGTGTCATCAGCATCGTTGCTGTCATATTCCACAAGGGACAGTATAAAGCGTTACGCGCCGTACGCCACCTCTGCGGAAACTGACGTTGCCGGGCTTCAAGCCGCCGTCAATGCTATGAACCACATCGTAGCCCTGTTGCAGCAGATGCTGCGCCGCGCCTTTGCTGCTGTGCCGTGATAACACATCACCATCACCGGGTGTCAAAGTCGTTATCACGCATAAAAGCGCCCAGCGTGTCGTTGGTTAAATGGAAAGCCTGCACCACTGTCCATTGCGAAACTCTGTGGATCCGCCTTGAATATCAGAGCCAGCACCGCCTCTTTCCTGCAACTTCTGGTGCGCGTCGGCAGCGTTAATGCATTGAACTGATCCATGCGTCTCTCTTTCTTTGACAAGTGGGCAGAATTACCGCACAGTTTACGTCGAAGCGGCAGATAAACGCCATAATGTTGCCATATCATAAAATGTTTTCAATGTTACCCCGCGATTCTTTGCTAATATGTTCGATAACGAACATTTATGAGCTTTAACGAAAGTGAATGAGGGCGGCATGACCAAAGATCTGATTGTGATGGGGGCGGCATCAATGGTGCTGGTATCGCGGCAGACGCCGCTGGACGCGGACCTCCGTGCTGATGCTGGGCGCAGGATCTCGCTTGCGGCCTCTTCCGCCAGTTCAAACTCATTCGGTGGCCTGCATACCTTGAGCACTATAATTCCGCTTTGGTCAGCGAGGCGCTGAACGTGAAGTGCTGCTGAAAATGGCCCCGCATATCGCGCCTTCCCGATGCGTTTCGCCTGCCATCGTCCGCATCTGCGCCCGGCGTGGATGATTCGCATTGGTCTGTTTATGTACGATCATCTGGGTAAACGCACCAGCTTGCCGGGATCAACTGGTTTGCGTTTTGGCGCAAATTCAGTGTTGAAATTAGCGCGGATTCCAGATATTCTGACTGTTGGTGGGCGACGCCCGTCTGGTACTCGCCAACGCCCAGATGGTGTGGTGCGTAAAGGCGGCGAAGTGCTACTCGGACTCGCGCCACCTCTGCTCGCCGCGAAACGGCCTGTGGATTGTGGAAGCGGAAGATCGATACCGGCAAAAATATAGCTGGCAAGCGCGCGGCAGTTAACGCGCCACCGGCCCGTGGGGTGAAACAGTTCTTCGACGACGGGATGCATCTGCCTTCGCCTTATGGCATTCGCCTGATCAAGGCAGCCATATTGTGGTGCCGCGCGTGCATGCAAGCAAGCCTACATTCTGCAAAACGAGATAAACGTATTGTGTTCGTGATCCCGTGGATGGACGAGTTTTCCATCATCGGCACTACCGATGTCGAGTACAAAGGCGAATCGAAAGCGGTGAAGATTGAAGAGTGAAATCAATTACCTGCTGAATAACACGCACTTTAAAAAGCCAGTTAAGCCATTGACGATATCGTCTGGACCTACTCCGGTGTGCGTCCGCTGTGTGATGATAGGTCGACTCGCCGCAGGCTATTACCCGTGATTACACCCTTGATATTCGCGGTGAAATGGCAAGCACCGCTGCTGTCGGTATTCGGCGGTAAGCTGACCACCTACCGAAAACTGGCGGAACATGCGCTGGAAAACTAACGCCGTATTATCAGGGTATTGGCCCGGCATGGACGAAAGAGAGTGTGCTACCGGGTGGCGCCATTGAAGCGACCGCGACGTTAATACCGCTCGCCTGCGCCGCCGCTATCCGTTCCTGACTGAATCGCTGGCGCGTCGTACGCTCGCACTTACGGCAGCAACAGCGAGCTGCTGCTCGGCAATGCAGGAACGGTAAGCGATCTCAGGGAAGATTTCGGTCATGGGTTCTACAAGCGGAGCTGAAATACGGTGGATGGGTCCACCGCCGACGACGCCCTGTAGTCGCACAAGTAAGGCATGTGGCTAAATGCGGATCAACAATCTCGTGTGAGTCGGTGGCTGGTGGAGTATACGCAGCAGGTTATCATGGCGTCGTAAATTAACGTAAGGTGATCGGTCAGATTTCGACTGGCCTGAGACTGATGACAAACCTACAAAACTGCCTGATGCGCTTCGCTTATCAGGCCTACGTAGTTTATGCAATATATTGAATTTGCATGGTCTTGTAGGCCAGATAAGACGTTCACGTCGCATCCGGCATGAACTCAGCACTTTGTCAAAAATCTAACCTACTTTTAATTCAGGGAATTACCGCAAAGCCCACATACCATCATGCAACGTAACAAAACTCAGGCACGTTCCCCTCGCCCCGAGAAAATAGCATTAATGCGCCCAGCGCCAGCATAAAAATTTTGAGCGGTGGTGTTGGCGTGATAATACAAACTAATAATACCGGCAAGTCCGACACCCAGCATGTAACCACCGCCAAAATTGCGCCAGTATGGGGATGCCGAAAAGTCATTACAACGAGGTCAAAATCCATTTCTGTTTTGCATTATTCTTCCATTCTTTTGAATGGTGAAGTGCTCCCGTGTTATACTTATGGACACTTTTCGAAATGATGGCGGAAAAACGGGACCGCTGGCCCCGTTCTGGCTGACCGGTGAACTTACAATCTCACCGGATCGATATGCCAGATATGATCGGCGTACTCTTTGATGGTACGGTCAGAAGAAGTAGCCCATATTGGCAATGTTCAACATCGCTTTGCGGTCCCACTCTTCCTGAAGCTCGTAGGGACATCGACTTTATCTGACAATCGACATAGCTGCGATAATCCATGAGTACGTATTGATCGCCGCCAGTTGATCAGCAGTCAACAGGTCGCGATAGCGACCGGATCTTCCGGACCCACCGCTGCCGATTTGCGTCGGCACCTGGTGCAGCTCCTCATCTTTCTCGTAGTATTCACGCGGTTGTAGCCCTGACGACGCGCAGTTCTTCCACTTCCGCTGTGTTACCAAAATAAGATATTGTCAGCGCCGACATGATCAGCATCTCGTATTCACCTCAACGCGATAGTCAGCGCACCGTTAAGCGCAAACTTCACGTATTACTGGTGCCGGAAGCTGCGTCCCTGCCAGCGAAATCTGTTCGAACAGATCTGCCGCCGGAATGATCAGCTGCGCCGGTTCTTGTAGTTCAGGATGAACACGACTTTCAACATCGCCAATCTGCGGATCGTTGTTGATCACTTTCGCTACGTCATTGATCAAATGAATAATGTGCTTCGCCATGTAATAGGCCGAAACCGCCTTACCGCCAAAATATTCACGCGCGGCCCACTTCGCATCGGTCGGCCTTGATGCGGTTGCGTAATCACATGCAACACATTCATCAATTGACGTTTGTATTCGTGAATACGTTTGGTTGTACATCGAACAACGCCTTTGGATTCACCACCACATTCAGCTGCTGGGCGATATACTCTGCCAGACGCTTTTGTTCTCCAGCCATGATGCACAGCGTGATTAACCATTGGGAAATCACAGTGTTGTTTTTGCAACTCATTAAGCAGGCTAAGGTCGGTGCCAGTTGCGGCCAGGTGTTCGTCCAGCACGGCTGAAGCGATGGGTTCGCTACCGCCAGCCAGCGACGCGGCGTACACCGTTGGTGACGTTGGTGAAACTGACCGGGAAGATTTCGCAAAGTCGGCAAACAACGTTGCACCATCAGTTAGAGTGCAGTTCCGATACACCGTTAACTGTGGCTCACAACAACCGCCAGCCGGGCCATACGCACGACGACCGTTGGATTCATCAATGATCGACGCCCGTCCCAGCGCAGATCGGTATCGTTCAGTTCTACTGTTCCTGCAAGTTTCAGAAATAGTCGTTGATTTCAAGATGATCTGCAGGTGACGCGGCAAATTTTACCAGCATATCAACCGGCGGGTTTCCAGCGCCTCGCTCATCAGCGTGGTTAGTGTAGGAAGACCTGACAACACACCTCAAACGCGTCGTCCAGCTAAATTGGTGCTCATCGATCAGCAGACGCATCCTCTCAGGAATCGACAATTGCGGATGGGTATCATTGAGATGATCGCGATTTTATCCGCCAGGTTATCGTAGGTTTTATGCAACTGATAATGGCGGCAAATGTCCCTGAATGGTCGAAACCGGGAAGTATTCTGACGCAGGCACAGCTCACGCCAGGTAAATTGAGTCATCCGGATACAACCGCGAGATACGTTCTCGAGTGGTTTTATCTTCTTCCACTGCCGCGAAGTAGTCACCTGGTTGAATTACCGAGTTAATTTCGCTAAACGCGCACTCCACAAACGCAGCGTGTTGGTCGCACGTCGGTGTCGTAACCGAGGATTATCGTAAGCGACTCCCAGAATCTCTTCGGTTTTCAATCCAGCGCGTTTTACCTTCCTGCTGAATGCACGTACGCCAAAACGGACTTTATATGGCGCGTGTTGTGGCGTTTGAATTCCACAGGTTACCGTATTCCAGCCAGTAGTCTGGCGACTCTTTCTGGCTACCGTTAACGATGTTCTGCTTGAACATACCCATGGTCATAGCGGATGCCGTAACCACGCCCAGCAGCCCTAACGTCGCCAGAATCAAAGGAAGCAAGCCGCCAGACGTCCCAGGCACCGTGCGAGGCACGGATCGACACCTTCATCAATCAGCTCTTCAGTTAACCCATCGCTTCCAGTGCGCCCTGTACATCTTCGTAAATTCTAGCGACAACATGGCGTTGGAGGCTTGTTAATCAAAACTCCATCGACTGGTATTAAACCTGATGAGTTTCTTACGACAACTGGGCACGGTTTGAACGTAACCAGCGCTCCGGACGATCGCGCACAGCAAATAACGTTGCGTTGTGATCATGTTGTGGCGACGACCGGGTTCCTTTCAATCGTAAACATCAGCTTGTAAGCGATAGAGTGTAGGCTTCTACGCGCTAAGCGTGGGCGATAGATATGTAAACGGAGCGATATATAAACCACAAACTAACACAAGCGATAGTAGCTCACGGTACGACTTCGCCGCGACCTGCCAGCTAAAATCCATTGCCATAGCCTGACGTTTGCACAAACCGCCACAGTGAAGGACGGGACCACAGTACAAAAGCACGTCCGAATAGCCCGTAACAGCGACCGGGCATTACTATCTTCAAGACAAAGCCACAGCGACGCCATCTGCAAGGTTCTCGAGAACGATCCAGAAACCGTGTCATAAACCCACCGGTGCACATAACGGCAACGTACCATGCTTCAATCCATAAAGTTGCGTTAAGCCGCACGGTTCAAGCGGCTGGGCACCAGGGCCAGCGTCCGCGCCGCCCTTAATGCGATGCGAAAATGCTTCGTGACCATCAGGCGCCCACCCTGACCGGGGTATTCCGCTGCCGCCGCAAAACCTTCCTGCAACGCACCGGATCGCCCGCGCCGAGGTAGCATAACGCCCTGCTCCAGAAGACCCGGTAGGCTTCCAGCACCGGGTCAGACGCTCTGGCTGGTCGGGCGGCTCACCACCGCAAAAGCGGCACTTTATCGTCAACGCCGCCCCATTGCGATTTGTAACTGGCGCTTATTTTCCGCTTTATCTTCAACGTATCGCGGGTGTATGAGGCCAACAGTAAGTCGATCTCTGGGCTCCGTTTTCTCCCACGCCGTTCAGTACGCCGGAAGACGCCCTTCACGGTGACAGCACGTAGACCTTCGCGCCACCGTGAGCAAACTGCGGTCGGTGATCTCGCGAGCGTGGTTGGACGACCGCCGTAATGCTTGATCGGCATGGTACAGACCGGCCTTCAGAAACT";
    event_table et = detect_events(rt, *ed_params);
    stList* kmer_list = build_kmer_list(sequence, sM->kmerLength, FALSE);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(kmer_list, *sM, et);
    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    stList* event_alignment = adaptive_banded_simple_event_align(et, sM, kmer_list);
    float start_time = fast5_get_start_time(fast5_handle);
    basecalled_event_table* b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    alignment_to_base_event_map(event_alignment, b_et, kmer_list, sM);
    char* prev_kmer = b_et->event[0].model_state;
    int move;
    int k = (int) sM->kmerLength;
    for (int i = 1; i < b_et->n; i++){
        move = b_et->event[i].move;
        if (strlen(prev_kmer) > 0) {
            CuAssertStrEquals(testCase, stString_getSubString(prev_kmer, move, k-move),
                              stString_getSubString(b_et->event[i].model_state, 0, k-move));
        }
        prev_kmer = b_et->event[i].model_state;
    }
    fast5_close(fast5_handle);

}

static void test_rna_alignment_to_base_event_map(CuTest *testCase){
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

    event_table et = detect_events(rt, *ed_params);
    event_table et2 = reverse_events(et);
    et = et2;
    stList* kmer_list = build_kmer_list(sequence, sM->kmerLength, TRUE);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(kmer_list, *sM, et);
    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    stList* event_alignment = adaptive_banded_simple_event_align(et, sM, kmer_list);
    float start_time = fast5_get_start_time(fast5_handle);
    basecalled_event_table* b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    rna_alignment_to_base_event_map(event_alignment, b_et, kmer_list, sM);
    basecalled_event_table *bet2 = reverse_basecalled_events(b_et);
    free(b_et);
    b_et = bet2;

    char* prev_kmer = b_et->event[0].model_state;
    int move;
    int k = (int) sM->kmerLength;
    for (int i = 1; i < b_et->n; i++){

//        printf("%i\n", b_et->event[i].move);
//        printf("%s\n", b_et->event[i].model_state);
//        printf("%" PRId64 "\n\n", b_et->event[i].raw_start);

        move = b_et->event[i].move;
        if (strlen(prev_kmer) > 0) {
            CuAssertStrEquals(testCase, stString_getSubString(prev_kmer, move, k-move),
                              stString_getSubString(b_et->event[i].model_state, 0, k-move));
        }
        prev_kmer = b_et->event[i].model_state;
    }

    fast5_close(fast5_handle);

}

static void test_estimate_scalings_using_mom(CuTest *testCase){
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling,
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
    stList* kmer_list = build_kmer_list(sequence, sM->kmerLength, TRUE);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(kmer_list, *sM, et);
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
    event_table et = detect_events(rt, *ed_params);
    event_table et2 = reverse_events(et);
    et = et2;
    stList* kmer_list = build_kmer_list(sequence, sM->kmerLength, TRUE);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(kmer_list, *sM, et);
    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    stList* something = adaptive_banded_simple_event_align(et, sM, kmer_list);
    struct AlignedPair *test_end = (struct AlignedPair*) stList_pop(something);
    CuAssertIntEquals(testCase, 453, test_end->ref_pos);
    CuAssertIntEquals(testCase, 1219, test_end->read_pos);
    fast5_close(fast5_handle);
}

static void test_load_from_raw_rna(CuTest *testCase) {
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    char* sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA";
    herr_t status = load_from_raw(path, templateModelFile, sequence, EVENT_LOCATION, true);
    CuAssertIntEquals(testCase, 0, (int) status);

    hid_t fast5_handle = fast5_open(path);

    int new_event_n = 1220; //how many events are generated (empirically determined)
    basecalled_event new_events[new_event_n + 1];
    fast5_get_basecall_events(fast5_handle, EVENT_LOCATION, new_events);

    int old_event_n = 1719; //how many events were called
    basecalled_event old_events[old_event_n + 1];
    fast5_get_basecall_events(fast5_handle, "/Analyses/Basecall_1D_000/BaseCalled_template/Events", old_events);

    // RNA strands are sequenced 3' -> 5' so, we need the reverse of the string
    CuAssertStrEquals(testCase, "AACCT", new_events[0].model_state);
    CuAssertStrEquals(testCase, "CCTAC", new_events[new_event_n - 1].model_state);
    CuAssertStrEquals(testCase, "AACCT", old_events[0].model_state);
    CuAssertStrEquals(testCase, "CCTAC", old_events[old_event_n - 1].model_state);

    status = H5Ldelete( fast5_handle, EVENT_LOCATION, H5P_DEFAULT );
    CuAssertIntEquals(testCase, 0, (int) status);
    fast5_close(fast5_handle);
}

static void test_load_from_raw_dna(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/1D/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch112_read108_strand.fast5");
    char* templateModelFile = stString_concat(HOME, "models/testModelR9p4_5mer_acegt_template.model");
    char* sequence = "TGCATGCCGTTTCCGTTACGTATTGCTAATCACGGGTGACGCCGTTTTCGCGCATTGAGCGAATCAGCAAAACCATCGCTAAACCACGGCTAACCCGGCGATGTGTGCTCCGTTCTCATCGACATCCCAAACCGTCAAACCATCCGGCGACAATCAGATCAGCGCAAAGATAATTAACCCACGTTGCAGGTAAATGCCACTTTGCGGATCGCGTTCGCCGTAGCCAGACGTAGCCCATCAGCACGCCACCACGCCAAGCCCGCCAAACCGGCCGCTGAATTTTGCTGCACATAGCCGCTTAACAGGGCGCTGATGGCGTAATGACGAGGTGGCTTACCGCTACCGAGGCGTTTTCCACCGCACCGCCGAGATACCACCACCAAGAGCAGGTTAAAGGGATCAGCGAAATTGCATTAACGCGTGGGTGGTAACGCCAGAACTCAAATTTCGTGTTGGATCCGAATGGCGGGGCCAGCCATAACGTCACTTCCTGATCGCCGGGGTACATGGCAATAAACACCACCACGCAGGCGATCATCATCACCCGGTTACCGGACCTGCGCGTTCACACCAAGGCGGCAAAGGATAACGGCGATAATGCAGGCCCTGCCGACCATGGCCTGCCTGCCAGCTCGCCGCCAGATAACACGCGCGGATCTGCCGGTTTCGAAAACGCGCCAGCTCCGCCCGTACGCTATCGGCACAGGACTCATCCGCCAGCCAGACATCGCTTTGGTTATGTTGTTGGTCGTGGGGATAATACCCTGCGTCGCCATGTAATCAACAAACGCCTGCGCCACGCGGGGATTGTAAAAGTGTCATCAGCATCGTTGCTGTCATATTCCACAAGGGACAGTATAAAGCGTTACGCGCCGTACGCCACCTCTGCGGAAACTGACGTTGCCGGGCTTCAAGCCGCCGTCAATGCTATGAACCACATCGTAGCCCTGTTGCAGCAGATGCTGCGCCGCGCCTTTGCTGCTGTGCCGTGATAACACATCACCATCACCGGGTGTCAAAGTCGTTATCACGCATAAAAGCGCCCAGCGTGTCGTTGGTTAAATGGAAAGCCTGCACCACTGTCCATTGCGAAACTCTGTGGATCCGCCTTGAATATCAGAGCCAGCACCGCCTCTTTCCTGCAACTTCTGGTGCGCGTCGGCAGCGTTAATGCATTGAACTGATCCATGCGTCTCTCTTTCTTTGACAAGTGGGCAGAATTACCGCACAGTTTACGTCGAAGCGGCAGATAAACGCCATAATGTTGCCATATCATAAAATGTTTTCAATGTTACCCCGCGATTCTTTGCTAATATGTTCGATAACGAACATTTATGAGCTTTAACGAAAGTGAATGAGGGCGGCATGACCAAAGATCTGATTGTGATGGGGGCGGCATCAATGGTGCTGGTATCGCGGCAGACGCCGCTGGACGCGGACCTCCGTGCTGATGCTGGGCGCAGGATCTCGCTTGCGGCCTCTTCCGCCAGTTCAAACTCATTCGGTGGCCTGCATACCTTGAGCACTATAATTCCGCTTTGGTCAGCGAGGCGCTGAACGTGAAGTGCTGCTGAAAATGGCCCCGCATATCGCGCCTTCCCGATGCGTTTCGCCTGCCATCGTCCGCATCTGCGCCCGGCGTGGATGATTCGCATTGGTCTGTTTATGTACGATCATCTGGGTAAACGCACCAGCTTGCCGGGATCAACTGGTTTGCGTTTTGGCGCAAATTCAGTGTTGAAATTAGCGCGGATTCCAGATATTCTGACTGTTGGTGGGCGACGCCCGTCTGGTACTCGCCAACGCCCAGATGGTGTGGTGCGTAAAGGCGGCGAAGTGCTACTCGGACTCGCGCCACCTCTGCTCGCCGCGAAACGGCCTGTGGATTGTGGAAGCGGAAGATCGATACCGGCAAAAATATAGCTGGCAAGCGCGCGGCAGTTAACGCGCCACCGGCCCGTGGGGTGAAACAGTTCTTCGACGACGGGATGCATCTGCCTTCGCCTTATGGCATTCGCCTGATCAAGGCAGCCATATTGTGGTGCCGCGCGTGCATGCAAGCAAGCCTACATTCTGCAAAACGAGATAAACGTATTGTGTTCGTGATCCCGTGGATGGACGAGTTTTCCATCATCGGCACTACCGATGTCGAGTACAAAGGCGAATCGAAAGCGGTGAAGATTGAAGAGTGAAATCAATTACCTGCTGAATAACACGCACTTTAAAAAGCCAGTTAAGCCATTGACGATATCGTCTGGACCTACTCCGGTGTGCGTCCGCTGTGTGATGATAGGTCGACTCGCCGCAGGCTATTACCCGTGATTACACCCTTGATATTCGCGGTGAAATGGCAAGCACCGCTGCTGTCGGTATTCGGCGGTAAGCTGACCACCTACCGAAAACTGGCGGAACATGCGCTGGAAAACTAACGCCGTATTATCAGGGTATTGGCCCGGCATGGACGAAAGAGAGTGTGCTACCGGGTGGCGCCATTGAAGCGACCGCGACGTTAATACCGCTCGCCTGCGCCGCCGCTATCCGTTCCTGACTGAATCGCTGGCGCGTCGTACGCTCGCACTTACGGCAGCAACAGCGAGCTGCTGCTCGGCAATGCAGGAACGGTAAGCGATCTCAGGGAAGATTTCGGTCATGGGTTCTACAAGCGGAGCTGAAATACGGTGGATGGGTCCACCGCCGACGACGCCCTGTAGTCGCACAAGTAAGGCATGTGGCTAAATGCGGATCAACAATCTCGTGTGAGTCGGTGGCTGGTGGAGTATACGCAGCAGGTTATCATGGCGTCGTAAATTAACGTAAGGTGATCGGTCAGATTTCGACTGGCCTGAGACTGATGACAAACCTACAAAACTGCCTGATGCGCTTCGCTTATCAGGCCTACGTAGTTTATGCAATATATTGAATTTGCATGGTCTTGTAGGCCAGATAAGACGTTCACGTCGCATCCGGCATGAACTCAGCACTTTGTCAAAAATCTAACCTACTTTTAATTCAGGGAATTACCGCAAAGCCCACATACCATCATGCAACGTAACAAAACTCAGGCACGTTCCCCTCGCCCCGAGAAAATAGCATTAATGCGCCCAGCGCCAGCATAAAAATTTTGAGCGGTGGTGTTGGCGTGATAATACAAACTAATAATACCGGCAAGTCCGACACCCAGCATGTAACCACCGCCAAAATTGCGCCAGTATGGGGATGCCGAAAAGTCATTACAACGAGGTCAAAATCCATTTCTGTTTTGCATTATTCTTCCATTCTTTTGAATGGTGAAGTGCTCCCGTGTTATACTTATGGACACTTTTCGAAATGATGGCGGAAAAACGGGACCGCTGGCCCCGTTCTGGCTGACCGGTGAACTTACAATCTCACCGGATCGATATGCCAGATATGATCGGCGTACTCTTTGATGGTACGGTCAGAAGAAGTAGCCCATATTGGCAATGTTCAACATCGCTTTGCGGTCCCACTCTTCCTGAAGCTCGTAGGGACATCGACTTTATCTGACAATCGACATAGCTGCGATAATCCATGAGTACGTATTGATCGCCGCCAGTTGATCAGCAGTCAACAGGTCGCGATAGCGACCGGATCTTCCGGACCCACCGCTGCCGATTTGCGTCGGCACCTGGTGCAGCTCCTCATCTTTCTCGTAGTATTCACGCGGTTGTAGCCCTGACGACGCGCAGTTCTTCCACTTCCGCTGTGTTACCAAAATAAGATATTGTCAGCGCCGACATGATCAGCATCTCGTATTCACCTCAACGCGATAGTCAGCGCACCGTTAAGCGCAAACTTCACGTATTACTGGTGCCGGAAGCTGCGTCCCTGCCAGCGAAATCTGTTCGAACAGATCTGCCGCCGGAATGATCAGCTGCGCCGGTTCTTGTAGTTCAGGATGAACACGACTTTCAACATCGCCAATCTGCGGATCGTTGTTGATCACTTTCGCTACGTCATTGATCAAATGAATAATGTGCTTCGCCATGTAATAGGCCGAAACCGCCTTACCGCCAAAATATTCACGCGCGGCCCACTTCGCATCGGTCGGCCTTGATGCGGTTGCGTAATCACATGCAACACATTCATCAATTGACGTTTGTATTCGTGAATACGTTTGGTTGTACATCGAACAACGCCTTTGGATTCACCACCACATTCAGCTGCTGGGCGATATACTCTGCCAGACGCTTTTGTTCTCCAGCCATGATGCACAGCGTGATTAACCATTGGGAAATCACAGTGTTGTTTTTGCAACTCATTAAGCAGGCTAAGGTCGGTGCCAGTTGCGGCCAGGTGTTCGTCCAGCACGGCTGAAGCGATGGGTTCGCTACCGCCAGCCAGCGACGCGGCGTACACCGTTGGTGACGTTGGTGAAACTGACCGGGAAGATTTCGCAAAGTCGGCAAACAACGTTGCACCATCAGTTAGAGTGCAGTTCCGATACACCGTTAACTGTGGCTCACAACAACCGCCAGCCGGGCCATACGCACGACGACCGTTGGATTCATCAATGATCGACGCCCGTCCCAGCGCAGATCGGTATCGTTCAGTTCTACTGTTCCTGCAAGTTTCAGAAATAGTCGTTGATTTCAAGATGATCTGCAGGTGACGCGGCAAATTTTACCAGCATATCAACCGGCGGGTTTCCAGCGCCTCGCTCATCAGCGTGGTTAGTGTAGGAAGACCTGACAACACACCTCAAACGCGTCGTCCAGCTAAATTGGTGCTCATCGATCAGCAGACGCATCCTCTCAGGAATCGACAATTGCGGATGGGTATCATTGAGATGATCGCGATTTTATCCGCCAGGTTATCGTAGGTTTTATGCAACTGATAATGGCGGCAAATGTCCCTGAATGGTCGAAACCGGGAAGTATTCTGACGCAGGCACAGCTCACGCCAGGTAAATTGAGTCATCCGGATACAACCGCGAGATACGTTCTCGAGTGGTTTTATCTTCTTCCACTGCCGCGAAGTAGTCACCTGGTTGAATTACCGAGTTAATTTCGCTAAACGCGCACTCCACAAACGCAGCGTGTTGGTCGCACGTCGGTGTCGTAACCGAGGATTATCGTAAGCGACTCCCAGAATCTCTTCGGTTTTCAATCCAGCGCGTTTTACCTTCCTGCTGAATGCACGTACGCCAAAACGGACTTTATATGGCGCGTGTTGTGGCGTTTGAATTCCACAGGTTACCGTATTCCAGCCAGTAGTCTGGCGACTCTTTCTGGCTACCGTTAACGATGTTCTGCTTGAACATACCCATGGTCATAGCGGATGCCGTAACCACGCCCAGCAGCCCTAACGTCGCCAGAATCAAAGGAAGCAAGCCGCCAGACGTCCCAGGCACCGTGCGAGGCACGGATCGACACCTTCATCAATCAGCTCTTCAGTTAACCCATCGCTTCCAGTGCGCCCTGTACATCTTCGTAAATTCTAGCGACAACATGGCGTTGGAGGCTTGTTAATCAAAACTCCATCGACTGGTATTAAACCTGATGAGTTTCTTACGACAACTGGGCACGGTTTGAACGTAACCAGCGCTCCGGACGATCGCGCACAGCAAATAACGTTGCGTTGTGATCATGTTGTGGCGACGACCGGGTTCCTTTCAATCGTAAACATCAGCTTGTAAGCGATAGAGTGTAGGCTTCTACGCGCTAAGCGTGGGCGATAGATATGTAAACGGAGCGATATATAAACCACAAACTAACACAAGCGATAGTAGCTCACGGTACGACTTCGCCGCGACCTGCCAGCTAAAATCCATTGCCATAGCCTGACGTTTGCACAAACCGCCACAGTGAAGGACGGGACCACAGTACAAAAGCACGTCCGAATAGCCCGTAACAGCGACCGGGCATTACTATCTTCAAGACAAAGCCACAGCGACGCCATCTGCAAGGTTCTCGAGAACGATCCAGAAACCGTGTCATAAACCCACCGGTGCACATAACGGCAACGTACCATGCTTCAATCCATAAAGTTGCGTTAAGCCGCACGGTTCAAGCGGCTGGGCACCAGGGCCAGCGTCCGCGCCGCCCTTAATGCGATGCGAAAATGCTTCGTGACCATCAGGCGCCCACCCTGACCGGGGTATTCCGCTGCCGCCGCAAAACCTTCCTGCAACGCACCGGATCGCCCGCGCCGAGGTAGCATAACGCCCTGCTCCAGAAGACCCGGTAGGCTTCCAGCACCGGGTCAGACGCTCTGGCTGGTCGGGCGGCTCACCACCGCAAAAGCGGCACTTTATCGTCAACGCCGCCCCATTGCGATTTGTAACTGGCGCTTATTTTCCGCTTTATCTTCAACGTATCGCGGGTGTATGAGGCCAACAGTAAGTCGATCTCTGGGCTCCGTTTTCTCCCACGCCGTTCAGTACGCCGGAAGACGCCCTTCACGGTGACAGCACGTAGACCTTCGCGCCACCGTGAGCAAACTGCGGTCGGTGATCTCGCGAGCGTGGTTGGACGACCGCCGTAATGCTTGATCGGCATGGTACAGACCGGCCTTCAGAAACT";
    hid_t fast5_handle = fast5_open(path);

    herr_t status = load_from_raw(path, templateModelFile, sequence, EVENT_LOCATION, false);
    CuAssertIntEquals(testCase, 0, (int) status);

    int new_event_n = 11020; //how many events are generated (empirically determined)
    basecalled_event new_events[new_event_n + 1];
    fast5_get_basecall_events(fast5_handle, EVENT_LOCATION, new_events);

    int old_event_n = 10922; //how many events were called
    basecalled_event old_events[old_event_n + 1];
    fast5_get_basecall_events(fast5_handle, "/Analyses/Basecall_1D_000/BaseCalled_template/Events", old_events);

    CuAssertStrEquals(testCase, "TGCAT", new_events[0].model_state);
    CuAssertStrEquals(testCase, "AAACT", new_events[new_event_n - 1].model_state);
    CuAssertStrEquals(testCase, "TGCAT", old_events[0].model_state);
    CuAssertStrEquals(testCase, "AAACT", old_events[old_event_n - 1].model_state);

    status = H5Ldelete( fast5_handle, EVENT_LOCATION, H5P_DEFAULT );
    CuAssertIntEquals(testCase, 0, (int) status);
    fast5_close(fast5_handle);
}

static void test_reverse_events(CuTest *testCase) {
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

    event_table et = detect_events(rt, *ed_params);
    event_table et2 = reverse_events(et);
    CuAssertIntEquals(testCase, (int) et.n, (int) et2.n);
    CuAssertIntEquals(testCase, (int) et.start, (int) et2.start);
    CuAssertIntEquals(testCase, (int) et.end, (int) et2.end);
    size_t n = et.n;
    for (int i = 0; i < n; i++){
        event_t event = et.event[n - 1 - i];
        event_t event2 = et2.event[i];
        CuAssertIntEquals(testCase, (int) event2.start, (int) event.start);
        CuAssertIntEquals(testCase, event2.pos, event.pos);
        CuAssertIntEquals(testCase, event2.state, event.state);
        CuAssertDblEquals(testCase, (double) event2.length, (double) event.length, 0.00001);
        CuAssertDblEquals(testCase, (double) event2.mean, (double) event.mean, 0.00001);
        CuAssertDblEquals(testCase, (double) event2.stdv, (double) event.stdv, 0.00001);
    }
    fast5_close(fast5_handle);

}

static void test_reverse_basecalled_events(CuTest *testCase){
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


    event_table et = detect_events(rt, *ed_params);
    event_table et2 = reverse_events(et);
    et = et2;

    char* templateModelFile = stString_concat(HOME, "/models/testModelR9p4_5mer_acgt_RNA.model");
    StateMachine *sM = stateMachine3_loadFromFile(templateModelFile, threeState, emissions_kmer_getGapProb,
                                                  emissions_signal_strawManGetKmerEventMatchProbWithDescaling_MeanOnly,
                                                  stateMachine3_loadTransitionsFromFile, NULL);
    const char* sequence = "CAUCCUGCCCUGUGUUAUCCAGUUAUGAGAUAAAAAAUGAAUAUAAGAGUGCUUGUCAUUAUAAAAGUUUUCCUUUUUAUUACCAUCCAAGCCACCAGCUGCCAGCCACCAGCAGCCAGCUGCCAGCACUAGCUUUUUUUUUUUAGCACUUAGUAUUUAGCAGCAUUUAUUAACAGGUACUUUAAGAAUGAUGAAGCAUUGUUUUAAUCUCACUGACUAUGAAGGUUUUAGUUUCUGCUUUUGCAAUUGUGUUUGUGAAAUUUGAAUACUUGCAGGCUUUGUAUGUGAAUAAUUUUAGCGGCUGGUUGGAGAUAAUCCUACGGGAAUUACUUAAAACUGUGCUUUAACUAAAAUGAAUGAGCUUUAAAAUCCCUCCUCCUACUCCAUCAUCAUCCCACUAUUCAUCUUAUCUCAUUAUCAUCAACCUAUCCCACAUCCCUAUCACCACAGCAAUCCAA";

    stList* kmer_list = build_kmer_list(sequence, sM->kmerLength, TRUE);

    NanoporeReadAdjustmentParameters test_squiggle_scalings = estimate_scalings_using_mom(kmer_list, *sM, et);
    update_SignalMachineWithNanoporeParameters(test_squiggle_scalings, sM);
    stList* something = adaptive_banded_simple_event_align(et, sM, kmer_list);
    float start_time = fast5_get_start_time(fast5_handle);
    basecalled_event_table* b_et = event_table_to_basecalled_table(&et, channel_params, start_time);
    rna_alignment_to_base_event_map(something, b_et, kmer_list, sM);
    basecalled_event_table *bet2 = reverse_basecalled_events(b_et);



    CuAssertIntEquals(testCase, (int) bet2->n, (int) b_et->n);
    CuAssertIntEquals(testCase, (int) bet2->start, (int) b_et->start);
    CuAssertIntEquals(testCase, (int) bet2->end, (int) b_et->end);
    CuAssertIntEquals(testCase, (int) bet2->aln_n, (int) b_et->aln_n);

    size_t n = b_et->n;
    for (int i = 0; i < n; i++){
        basecalled_event event = b_et->event[n - 1 - i];
        basecalled_event event2 = bet2->event[i];
        CuAssertIntEquals(testCase, (int) event2.raw_start, (int) event.raw_start);
        CuAssertIntEquals(testCase, (int) event2.raw_length, (int) event.raw_length);
        CuAssertIntEquals(testCase, event2.move, event.move);
        CuAssertDblEquals(testCase, event2.start, event.start, 0.00001);
        CuAssertDblEquals(testCase, event2.length, event.length, 0.00001);
        CuAssertDblEquals(testCase, event2.mean, event.mean, 0.00001);
        CuAssertDblEquals(testCase, event2.stdv, event.stdv, 0.00001);
        CuAssertDblEquals(testCase, event2.p_model_state, event.p_model_state, 0.00001);
        CuAssertStrEquals(testCase, event2.model_state, event.model_state);

    }
    fast5_close(fast5_handle);

}

static void test_build_kmer_list(CuTest *testCase){
    // DNA - ATGCATGC -> ATGCA, TGCAT, GCATG, CATGC
    // RNA - AUGCAUGC -> ACGTA, TACGT, GTACG, CGTAC
    const char* dna_sequence = "ATGCATGC";
    char* expected_dna_kmers[4] = {"ATGCA", "TGCAT", "GCATG", "CATGC"};
    stList* dna_kmers = build_kmer_list(dna_sequence, 5, FALSE);
    for (int x = 0; x < 4; x++){
        CuAssertStrEquals(testCase, expected_dna_kmers[x], stList_get(dna_kmers, x));
    }

    const char* rna_sequence = "AUGCAUGC";
    char* expected_rna_kmers[4] = {"ACGTA", "TACGT", "GTACG", "CGTAC"};
    stList* rna_kmers = build_kmer_list(rna_sequence, 5, TRUE);
    for (int x = 0; x < 4; x++){
        CuAssertStrEquals(testCase, expected_rna_kmers[x], stList_get(rna_kmers, x));
    }
}

static void test_fast5_create_group(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);

    herr_t status = fast5_create_group(fast5_handle, "/Analyses/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    status = fast5_create_group(fast5_handle, "/Analyses/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/Analyses/FakePath");
    status = fast5_create_group(fast5_handle, "/Analyses/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/Analyses/FakePath");
    status = fast5_create_group(fast5_handle, "/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/FakePath");
//    clean up
    fast5_close(fast5_handle);

}

static void test_fast5_create_all_groups(CuTest *testCase) {
    char *path = stString_concat(HOME,
                                 "tests/minion_test_reads/RNA_edge_cases/DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_61_ch_151_strand.fast5");
    hid_t fast5_handle = fast5_open(path);

    herr_t status = fast5_create_all_groups(fast5_handle, "/Analyses/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    status = fast5_create_all_groups(fast5_handle, "/Analyses/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/Analyses/FakePath");
    status = fast5_create_all_groups(fast5_handle, "/Analyses/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/Analyses/FakePath");
    status = fast5_create_all_groups(fast5_handle, "/FakePath");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/FakePath");
    status = fast5_create_all_groups(fast5_handle, "/FakePath/depth");
    CuAssertIntEquals(testCase, 0, status);
    H5Gunlink(fast5_handle, "/FakePath/depth");
    H5Gunlink(fast5_handle, "/FakePath");

//    clean up
    fast5_close(fast5_handle);
}

static void test_fast5_get_string(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/embedded_files/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    char* fastq = fast5_get_string(fast5_handle, "Analyses/Basecall_1D_000/BaseCalled_template/Fastq");

    CuAssertStrEquals(testCase, "@5cc8", stString_getSubString(fastq, 0, 5));

//    clean up
    fast5_close(fast5_handle);
    free(fastq);
}

static void test_fast5_get_fastq(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/embedded_files/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    char* fastq = fast5_get_fastq(fast5_handle);

    CuAssertStrEquals(testCase, "@5cc8", stString_getSubString(fastq, 0, 5));
//    clean up
    fast5_close(fast5_handle);
    free(fastq);
    free(path);
}

static void test_hdf5_group_exists(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/embedded_files/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5");
    hid_t fast5_handle = fast5_open(path);

    bool test = hdf5_group_exists(fast5_handle, "/Analyses/Basecall_1D_000");
    CuAssertTrue(testCase, test);
    test = hdf5_group_exists(fast5_handle, "/Analyses/Basecall_1D_0");
    CuAssertTrue(testCase, !test);

    fast5_close(fast5_handle);
}

static void test_write_fastqs_to_file(CuTest *testCase){
    char* dir = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/");
    char* out_file = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/test.fastq");

    int pass = write_fastqs_to_file(dir, out_file);
    CuAssertTrue(testCase, pass == 0);
    if (pass == 0){
        CuAssertTrue(testCase, stFile_exists(out_file));
        stFile_rmrf(out_file);
    } else {
        if (stFile_exists(out_file)){
            stFile_rmrf(out_file);
        }
    }
//    clean up
    free(out_file);
    free(dir);
}

static void test_check_file_ext(CuTest *testCase){
    char* file_path = "asdfasdf.dfsdf";
    char* file_path2 = "asd.fasdf.dfsdf";
    bool test;

    test = check_file_ext(file_path, "dfsdf");
    CuAssertTrue(testCase, test);
    test = check_file_ext(file_path2, "dfsdf");
    CuAssertTrue(testCase, test);

}

static void test_path_join_two_strings(CuTest *testCase){
    char* file_path = "asdf/asdf/asdf";
    char* file_path2 = "/asdf/asdf/asdf/";
    char* test;

    test = path_join_two_strings(file_path, "dfsdf");
    CuAssertStrEquals(testCase, test, "asdf/asdf/asdf/dfsdf");
    test = path_join_two_strings(file_path2, "dfsdf");
    CuAssertStrEquals(testCase, test, "/asdf/asdf/asdf/dfsdf");
    free(test);
}

static void test_parse_fastq_string(CuTest *testCase){
    char* path = stString_concat(HOME, "tests/minion_test_reads/embedded_files/LomanLabz_PC_20161025_FNFAB42699_MN17633_sequencing_run_20161025_E_coli_native_450bps_82361_ch92_read1108_strand.fast5");
    hid_t fast5_handle = fast5_open(path);
    char* fastq_str = fast5_get_fastq(fast5_handle);
    fastq_entry *fastq_data = parse_fastq_string(fastq_str);
    CuAssertStrEquals(testCase, "5cc8", stString_getSubString(fastq_data->name, 0, 4));
    fast5_close(fast5_handle);
    free(fastq_str);
//    free(fastq_data);
    fastq_entry_destruct(fastq_data);

}

static void test_write_readdb_file1(CuTest *testCase){
    char* dir = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/");
    char* out_file = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/test.readdb");
    int pass = write_readdb_file1(dir, out_file);
    CuAssertTrue(testCase, pass == 0);
    if (pass == 0){
        CuAssertTrue(testCase, stFile_exists(out_file));
        stFile_rmrf(out_file);
    } else {
        if (stFile_exists(out_file)){
            stFile_rmrf(out_file);
        }
    }
//    clean up
    free(out_file);
    free(dir);

}


static void test_write_fastq_and_readdb_file1(CuTest *testCase){
    char* dir = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/");
    char* out_file = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/test.readdb");
    char* out_file2 = stString_concat(HOME, "tests/minion_test_reads/canonical_ecoli_R9/test.fastq");

    int pass = write_fastq_and_readdb_file1(dir, out_file2, out_file);
    CuAssertTrue(testCase, pass == 0);
    if (pass == 0){
        CuAssertTrue(testCase, stFile_exists(out_file));
        CuAssertTrue(testCase, stFile_exists(out_file2));
        stFile_rmrf(out_file);
        stFile_rmrf(out_file2);
    } else {
        if (stFile_exists(out_file)){
            stFile_rmrf(out_file);
        }
        if (stFile_exists(out_file2)){
            stFile_rmrf(out_file2);
        }
    }
//    clean up
    free(out_file);
    free(out_file2);
    free(dir);
}



CuSuite *eventAlignerTestSuite(void) {

    CuSuite *suite = CuSuiteNew();

    SUITE_ADD_TEST(suite, test_fast5_get_raw_read_name);
    SUITE_ADD_TEST(suite, test_fast5_get_raw_read_group);
    SUITE_ADD_TEST(suite, test_fast5_get_fixed_string_attribute);
    SUITE_ADD_TEST(suite, test_fast5_get_read_id);
    SUITE_ADD_TEST(suite, test_fast5_get_experiment_type);
    SUITE_ADD_TEST(suite, test_fast5_read_float_attribute);
    SUITE_ADD_TEST(suite, test_fast5_get_channel_params);
    SUITE_ADD_TEST(suite, test_fast5_get_raw_samples);
    SUITE_ADD_TEST(suite, test_fast5_get_start_time);
    SUITE_ADD_TEST(suite, test_fast5_set_basecall_event_table);
    SUITE_ADD_TEST(suite, test_fast5_create_group);
    SUITE_ADD_TEST(suite, test_fast5_create_all_groups);

    SUITE_ADD_TEST(suite, test_set_NanoporeReadAdjustmentParameters);
    SUITE_ADD_TEST(suite, test_estimate_scalings_using_mom);
    SUITE_ADD_TEST(suite, test_update_SignalMachineWithNanoporeParameters);
    SUITE_ADD_TEST(suite, test_sonlib_lists);
    SUITE_ADD_TEST(suite, test_adaptive_banded_simple_event_align);
    SUITE_ADD_TEST(suite, test_event_table_to_basecalled_table);

    SUITE_ADD_TEST(suite, test_rna_alignment_to_base_event_map);
    SUITE_ADD_TEST(suite, test_alignment_to_base_event_map);

    SUITE_ADD_TEST(suite, test_load_from_raw_rna);
    SUITE_ADD_TEST(suite, test_load_from_raw_dna);

    SUITE_ADD_TEST(suite, test_reverse_events);
    SUITE_ADD_TEST(suite, test_reverse_basecalled_events);
    SUITE_ADD_TEST(suite, test_build_kmer_list);

    SUITE_ADD_TEST(suite, test_fast5_get_string);
    SUITE_ADD_TEST(suite, test_hdf5_group_exists);
    SUITE_ADD_TEST(suite, test_fast5_get_fastq);
    SUITE_ADD_TEST(suite, test_check_file_ext);
    SUITE_ADD_TEST(suite, test_path_join_two_strings);
    SUITE_ADD_TEST(suite, test_write_fastqs_to_file);
    SUITE_ADD_TEST(suite, test_parse_fastq_string);
    SUITE_ADD_TEST(suite, test_write_readdb_file1);
    SUITE_ADD_TEST(suite, test_write_fastq_and_readdb_file1);

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
//
//KSEQ_INIT(gzFile, gzread)
//
//int main(int argc, char *argv[])
//{
//    gzFile fp;
//    kseq_t *seq;
//    int l;
//    if (argc == 1) {
//        fprintf(stderr, "Usage: %s <in.seq>\n", argv[0]);
//        return 1;
//    }
//    fp = gzopen(argv[1], "r"); // STEP 2: open the file handler
//    seq = kseq_init(fp); // STEP 3: initialize seq
//    while ((l = kseq_read(seq)) >= 0) { // STEP 4: read sequence
//        printf("name: %s\n", seq->name.s);
//        if (seq->comment.l) printf("comment: %s\n", seq->comment.s);
//        printf("seq: %s\n", seq->seq.s);
//        if (seq->qual.l) printf("qual: %s\n", seq->qual.s);
//    }
//    printf("return value: %d\n", l);
//    kseq_destroy(seq); // STEP 5: destroy seq
//    gzclose(fp); // STEP 6: close the file handler
//    return 0;
//}
