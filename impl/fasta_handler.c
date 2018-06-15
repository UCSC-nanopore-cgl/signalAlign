//
// Created by Andrew Bailey on 4/20/18. 4/20 ... blaze it.
//

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "htslib/faidx.h"
#include "sonLib.h"
#include "signalMachineUtils.h"
#include "htslib/hfile.h"




char *fastaHandler_getSubSequence(char *fastaReferencePath, int64_t start, int64_t end, bool strand,
                                const char *sequence_name) {
    // index reference
    faidx_t *fasta_index = fai_load(fastaReferencePath);
    int length;

    if (strand) {

        char *sequence = faidx_fetch_seq(fasta_index, sequence_name, start, end-1, &length);
        if (length == -2) {
            st_errAbort("[signalMachine] ERROR %d: sequence name: %s is not in reference fasta: %s \n",
                        length, sequence_name, fastaReferencePath);
        }
        fai_destroy(fasta_index);
        return sequence;

    }

    char *sequence = faidx_fetch_seq(fasta_index, sequence_name, end, start-1, &length);
    if (length == -2) {
        st_errAbort("[signalMachine] ERROR %d: sequence name: %s is not in reference fasta: %s \n",
                    length, sequence_name, fastaReferencePath);
    }
    fai_destroy(fasta_index);

    return sequence;


}


ReferenceSequence *fastaHandler_ReferenceSequenceConstructFull(char *forward_fastaReferencePath,
                                                              char *backward_fastaReferencePath,
                                                              struct PairwiseAlignment *pA,
                                                              const char *sequence_name) {
    ReferenceSequence *R = st_malloc(sizeof(ReferenceSequence));
    char *backward_sequence = NULL;

    R->A = referenceSequence_copyPairwiseAlignment(pA);

    char *forward_sequence = fastaHandler_getSubSequence(forward_fastaReferencePath, R->A->start1, R->A->end1,
                                                       R->A->strand1, sequence_name);

    if (forward_sequence == NULL) {
        st_errAbort("[signalMachine] ERROR: Unable to fetch reference sequence.  \n");
    }
    // if backward_fastaReference is needed use it to create backward sequence
    if (backward_fastaReferencePath){
        backward_sequence = signalUtils_stringReverse(fastaHandler_getSubSequence(
                backward_fastaReferencePath, R->A->start1, R->A->end1, R->A->strand1, sequence_name));
    } else {
        backward_sequence = signalUtils_stringReverse(stString_ComplementString(forward_sequence));
    }

    R->trimmedForwardSequence = forward_sequence;
    R->trimmedBackwardSequence = backward_sequence;
    // entire reference sequence is not needed
    R->reference = stString_copy("reference");
    R->complementOfReference = stString_copy("complementOfReference");

    R->getTemplateTargetSequence = referenceSequence_getTemplateTarget;
    R->getComplementTargetSequence = referenceSequence_getComplementTarget;

    R->initialized = TRUE;

    return R;
}



void build_fai_index_file() {
//    int bugs = 100;
//    double bug_rate = 1.2;
    char fast_path[] = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/test_sequences/fake_rna_reversed.fa";
//    fprintf(stdout, "You have %d bugs at the imaginary rate of %f.\n", bugs, bug_rate);
    fprintf(stdout, "fasta path: %s. \n", fast_path);
//    int success = fai_build(fast_path);
//    if(success==0){
//        fprintf(stdout, "Indexed Fasta\n");
//    }

    char *test1 = "asdf";
    char *test2 = stString_copy(test1);


    char *index_path = stString_concat(fast_path, ".fai");
    fprintf(stdout, "sequence name %s \n", index_path);
    hFILE *fp = NULL;

    fp = hopen(index_path, "rb");
    if (fp == 0) {
        fprintf(stdout, "NULL SHIT I THINK \n");
    }
    faidx_t *fasta_thing = fai_load3(fast_path, index_path, NULL, FAI_CREATE);


    fprintf(stdout, "loaded Fasta \n");

    const char *sequence_name = faidx_iseq(fasta_thing, 0);
    fprintf(stdout, "sequence name %s \n", sequence_name);

    int length = 11;
    int *l = &length;

//    fprintf(stdout, "length %d \n", length);

    char *nucleotide_sequence = faidx_fetch_seq(fasta_thing, sequence_name, 0, 10, l);
    fprintf(stdout, "Forward Sequence: %s \n", nucleotide_sequence);

//    char *rev_nucleotide_sequence = faidx_fetch_seq(fasta_thing, sequence_name, 10, 0, l);
//    fprintf(stdout, "No seg fault: %s \n", rev_nucleotide_sequence);

    char *revers_comp = stString_reverseComplementString(nucleotide_sequence);
    fprintf(stdout, "rev-comp sequence: %s \n", revers_comp);
    char *comp = stString_ComplementString(nucleotide_sequence);
    fprintf(stdout, "comp sequence: %s \n", comp);

    signalUtils_stringReverse(nucleotide_sequence);
    fprintf(stdout, "Reversed sequence: %s \n", nucleotide_sequence);
//    char c = 'C';
    char a = stString_reverseComplementChar('C');
    fprintf(stdout, "ReverseCompChar: %c:%c \n", a, 'C');
    char b = stString_reverseComplementChar('t');
    fprintf(stdout, "ReverseCompChar: %c:%c \n", b, 't');
    char c = stString_reverseComplementChar('G');
    fprintf(stdout, "ReverseCompChar: %c:%c \n", c, 'g');
    char d = stString_reverseComplementChar('a');
    fprintf(stdout, "ReverseCompChar: %c:%c \n", d, 'a');


//    char exonerateCigarFile[] =
//            "/Users/andrewbailey/data/signal_align_output/tempFiles_alignment/tempFiles_DEAMERNANOPORE_20170922_FAH26525_MN16450_sequencing_run_MA_821_R94_NA12878_mRNA_09_22_17_67136_read_36_ch_218_strand/temp_cigar_7d31de25-8c15-46d8-a08c-3d5043258c89.txt";
//    char forward_path[] =
//            "/Users/andrewbailey/data/signal_align_output/tempFiles_alignment/rna_fake_reversed.None.forward.txt";
//    char backward_path[] =
//            "/Users/andrewbailey/data/signal_align_output/tempFiles_alignment/rna_fake_reversed.None.backward.txt";

//    ReferenceSequence *R = st_malloc(sizeof(ReferenceSequence));
//    FILE *fileHandleIn = fopen(exonerateCigarFile, "r");
    // parse input CIGAR to get anchors
//    struct PairwiseAlignment *pA;
//    pA = cigarRead(fileHandleIn);
//    ReferenceSequence *R1;
//    ReferenceSequence *R2;

//    R1 = fastaIndex_ReferenceSequenceConstructFull(fast_path, pA, sequence_name);
//
//    R2 = signalUtils_ReferenceSequenceConstructFull(forward_path, backward_path, pA);
//    signalUtils_ReferenceSequenceDestruct(R2);
//    signalUtils_ReferenceSequenceDestruct(R1);
    free(index_path);
    fai_destroy(fasta_thing);
}