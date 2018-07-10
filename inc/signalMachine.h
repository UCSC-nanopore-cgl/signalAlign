//
// Created by Andrew Bailey on 6/26/18.
//

#ifndef NANOPORE_RNN_SIGNALMACHINE_H
#define NANOPORE_RNN_SIGNALMACHINE_H

#include "signalMachineUtils.h"


StateMachine *buildStateMachine(const char *modelFile, NanoporeReadAdjustmentParameters npp, StateMachineType type,
                                NanoporeHDP *nHdp);

#endif //NANOPORE_RNN_SIGNALMACHINE_H
