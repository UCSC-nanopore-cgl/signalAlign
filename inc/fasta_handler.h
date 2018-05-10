//
// Created by Andrew Bailey on 4/20/18.
//

#ifndef NANOPORE_RNN_FASTA_HANDLER_H
#define NANOPORE_RNN_FASTA_HANDLER_H


#include <stdio.h>
#include <stdlib.h>
#include "htslib/faidx.h"
#include "pairwiseAlignment.h"
#include "signalMachineUtils.h"
#include "sonLib.h"


void build_fai_index_file();

char *fastaHandler_getSubSequence(char *fastaReferencePath, int64_t start, int64_t end, bool strand,
                                char *sequence_name);


ReferenceSequence *fastaHandler_ReferenceSequenceConstructFull(char *forward_fastaReferencePath,
                                                             char *backward_fastaReferencePath,
                                                             struct PairwiseAlignment *pA,
                                                             const char *sequence_name);


#endif //NANOPORE_RNN_FASTA_HANDLER_H
