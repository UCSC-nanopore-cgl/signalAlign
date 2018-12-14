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
                                                              const char *sequence_name,
                                                              bool rna) {
    ReferenceSequence *R = st_malloc(sizeof(ReferenceSequence));
    char *backward_sequence = NULL;
    if (rna) {
        listReverse(pA->operationList);
    }

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
    if (rna){
        char *tmp = backward_sequence;
//        printf("%s\n", tmp);
        backward_sequence = stString_ReverseString(forward_sequence);
        forward_sequence = stString_ReverseString(tmp);
        int64_t tmp2 = R->A->start1;
        R->A->start1 = R->A->end1;
        R->A->end1 = tmp2;
        pA->end1 = R->A->end1;
        pA->start1 = R->A->start1;
        pA->strand1 = !pA->strand1;
        R->A->strand1 = !R->A->strand1;

//        printf("%s\n", forward_sequence);
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

