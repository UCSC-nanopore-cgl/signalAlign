#include <stdio.h>
#include <stdint.h>
#include <stdlib.h>
#include "nanopore.h"
#include "pairwiseAligner.h"


void printNanoporeReadAdjustmentParameters(NanoporeReadAdjustmentParameters *npRead) {
    printf("scale %lf\n", npRead->scale);
    printf("shift %lf\n", npRead->shift);
    printf("var %lf\n", npRead->var);
    printf("scale_sd %lf\n", npRead->scale_sd);
    printf("shift_sd %lf\n", npRead->shift_sd);
    printf("var_sd %lf\n", npRead->var_sd);
    printf("drift %lf\n", npRead->drift);
}


static NanoporeRead *nanopore_NanoporeReadConstruct(int64_t readLength,
                                                    int64_t nbTemplateEvents,
                                                    int64_t nbComplementEvents,
                                                    int64_t templateReadLength,
                                                    int64_t complementReadLength) {
    NanoporeRead *npRead = st_malloc(sizeof(NanoporeRead));

    npRead->readLength = readLength;
    npRead->nbTemplateEvents = nbTemplateEvents;
    npRead->nbComplementEvents = nbComplementEvents;
    npRead->templateReadLength = templateReadLength;
    npRead->complementReadLength = complementReadLength;

    npRead->twoDread = st_malloc((npRead->readLength + 1) * sizeof(char));
    npRead->templateRead = st_malloc((npRead->templateReadLength + 1) * sizeof(char));
    npRead->complementRead= st_malloc((npRead->complementReadLength+ 1) * sizeof(char));

    // the map contains the index of the event corresponding to each kmer in the read sequence so
    // the length of the map has to be the same as the read sequence, not the number of events
    // there can be events that aren't mapped to a kmer
    npRead->templateEventMap = st_malloc(npRead->readLength * sizeof(int64_t));
    npRead->templateStrandEventMap = st_malloc(npRead->templateReadLength * sizeof(int64_t));
    npRead->templateEvents = st_malloc(npRead->nbTemplateEvents * NB_EVENT_PARAMS * sizeof(double));

    npRead->complementEventMap = st_malloc(npRead->readLength * sizeof(int64_t));
    npRead->complementStrandEventMap = st_malloc(npRead->complementReadLength * sizeof(int64_t));
    npRead->complementEvents = st_malloc(npRead->nbComplementEvents * NB_EVENT_PARAMS * sizeof(double));

    npRead->templateModelState = st_malloc(npRead->nbTemplateEvents * sizeof(int64_t));
    npRead->templatePModel = st_malloc(npRead->nbTemplateEvents * sizeof(double));

    npRead->complementModelState = st_malloc(npRead->nbComplementEvents * sizeof(int64_t));
    npRead->complementPModel = st_malloc(npRead->nbComplementEvents * sizeof(double));

    npRead->scaled = TRUE;

    // initialize to 0.0 these params are not in the ONT models
    npRead->templateParams.shift_sd = 0.0;
    npRead->complementParams.shift_sd = 0.0;

    return npRead;
}

static inline double nanopore_getEventMean(double *events, int64_t index) {
    //return *(events + (index * NB_EVENT_PARAMS));
    return events[(index * NB_EVENT_PARAMS)];
}

static inline double nanopore_getEventSd(double *events, int64_t index) {
    //return *(events + ((index * NB_EVENT_PARAMS) + sizeof(double)));
    return events[(1 + (index * NB_EVENT_PARAMS))];
}

static inline double nanopore_getEventDuration(double *events, int64_t index) {
    //return *(events + ((1 + index * NB_EVENT_PARAMS) + (sizeof(double))));
    return events[2 + (index * NB_EVENT_PARAMS)];
}

static inline double nanopore_getEventDeltaTime(double *events, int64_t index) {
    //return *(events + ((1 + index * NB_EVENT_PARAMS) + (sizeof(double))));
    return events[3 + (index * NB_EVENT_PARAMS)];
}

static void nanopore_descaleEvents(int64_t nb_events, double *events, double scale, double shift) {
    for (int64_t i = 0; i < nb_events; i += NB_EVENT_PARAMS) {
        events[i] = (events[i] - shift) / scale;
    }
}

EventKmerTuple *nanopore_eventKmerTupleConstruct(double mean, double sd, double deltaTime, int64_t kmerIndex) {
    EventKmerTuple *t = st_malloc(sizeof(EventKmerTuple));
    t->eventMean = mean;
    t->eventSd = sd;
    t->kmerIndex = kmerIndex;
    t->deltaTime = deltaTime;
    return t;
}

NanoporeReadAdjustmentParameters *nanopore_readAdjustmentParametersConstruct() {
    NanoporeReadAdjustmentParameters *params = st_malloc(sizeof(NanoporeReadAdjustmentParameters));
    params->scale = 0.0;
    params->scale_sd = 0.0;
    params->shift = 0.0;
    params->shift_sd = 0.0;
    params->var = 0.0;
    params->var_sd = 0.0;
    params->drift = 0.0;
    return params;
}


static NanoporeReadAdjustmentParameters const NanoporeReadAdjustmentParameters_default = {
        .scale = 1.0,
        .shift = 0.0,
        .drift = 0.0,
        .var = 1.0,
        .scale_sd = 1.0,
        .var_sd = 1.0,
        .shift_sd = 0.0
};

NanoporeReadAdjustmentParameters set7_NanoporeReadAdjustmentParameters(double _shift, double _scale, double _drift,
                                                                       double _var, double _scale_sd, double _var_sd,
                                                                       double _shift_sd){
    NanoporeReadAdjustmentParameters default_scaling = NanoporeReadAdjustmentParameters_default;
    // direct
    default_scaling.shift = _shift;
    default_scaling.scale = _scale;
    default_scaling.drift = _drift;
    default_scaling.var = _var;
    default_scaling.scale_sd = _scale_sd;
    default_scaling.var_sd = _var_sd;
    default_scaling.shift_sd = _shift_sd;

    return default_scaling;
}

NanoporeReadAdjustmentParameters set4_NanoporeReadAdjustmentParameters(double _shift, double _scale, double _drift, double _var) {

    NanoporeReadAdjustmentParameters set4_scalings = set7_NanoporeReadAdjustmentParameters(_shift, _scale, _drift, _var, 1.0, 1.0, 0.0);
    return set4_scalings;
}



NanoporeRead *nanopore_loadNanoporeReadFromFile(const char *nanoporeReadFile) {
    FILE *fH = fopen(nanoporeReadFile, "r");
    // line 1: all tab-seperated
    // 0 alignment read length
    // 1 # template events
    // 2 # complement events
    // 3 length of template read
    // 4 length of complement read
    // 5 template scale
    // 6 template shift
    // 7 template var
    // 8 template scale_sd
    // 9 template var_sd
    // 0 template drift
    // 1 complement scale
    // 2 complement shift
    // 3 complement var
    // 4 complement scale_sd
    // 5 complement var_sd
    // 6 complement drift

    char *string = stFile_getLineFromFile(fH);
    stList *tokens = stString_split(string);
    int64_t npRead_readLength, npRead_nbTemplateEvents, npRead_nbComplementEvents,
            templateReadLength, complementReadLength;

    if (stList_length(tokens) != 18) {
        st_errAbort("nanopore_loadNanoporeReadFromFile: incorrect line 1 for file %s got %"PRId64"tokens should get 18",
                    nanoporeReadFile, stList_length(tokens));

    }

    int64_t j = sscanf(stList_get(tokens, 0), "%"SCNd64, &npRead_readLength);
    if (j != 1) {
        st_errAbort("error parsing nanopore read length\n");
    }
    j = sscanf(stList_get(tokens, 1), "%"SCNd64, &npRead_nbTemplateEvents);
    if (j != 1) {
        st_errAbort("error parsing nanopore template event lengths\n");
    }
    j = sscanf(stList_get(tokens, 2), "%"SCNd64, &npRead_nbComplementEvents);
    if (j != 1) {
        st_errAbort("error parsing nanopore complement event lengths\n");
    }
    j = sscanf(stList_get(tokens, 3), "%"SCNd64, &templateReadLength);
    if (j != 1) {
        st_errAbort("error parsing nanopore template read length\n");
    }
    j = sscanf(stList_get(tokens, 4), "%"SCNd64, &complementReadLength);
    if (j != 1) {
        st_errAbort("error parsing nanopore complement read length\n");
    }

    NanoporeRead *npRead = nanopore_NanoporeReadConstruct(npRead_readLength,
                                                          npRead_nbTemplateEvents,
                                                          npRead_nbComplementEvents,
                                                          templateReadLength,
                                                          complementReadLength);

    // template params
    j = sscanf(stList_get(tokens, 5), "%lf", &(npRead->templateParams.scale));
    if (j != 1) {
        st_errAbort("error parsing nanopore template scale\n");
    }
    j = sscanf(stList_get(tokens, 6), "%lf", &(npRead->templateParams.shift));
    if (j != 1) {
        st_errAbort("error parsing nanopore template shift\n");
    }
    j = sscanf(stList_get(tokens, 7), "%lf", &(npRead->templateParams.var));
    if (j != 1) {
        st_errAbort("error parsing nanopore template var\n");
    }
    j = sscanf(stList_get(tokens, 8), "%lf", &(npRead->templateParams.scale_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore template scale_sd\n");
    }
    j = sscanf(stList_get(tokens, 9), "%lf", &(npRead->templateParams.var_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore template var_sd\n");
    }
    j = sscanf(stList_get(tokens, 10), "%lf", &(npRead->templateParams.drift));
    if (j != 1) {
        st_errAbort("error parsing nanopore template drift\n");
    }
    // complement params
    j = sscanf(stList_get(tokens, 11), "%lf", &(npRead->complementParams.scale));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement scale\n");
    }
    j = sscanf(stList_get(tokens, 12), "%lf", &(npRead->complementParams.shift));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement shift\n");
    }
    j = sscanf(stList_get(tokens, 13), "%lf", &(npRead->complementParams.var));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement var\n");
    }
    j = sscanf(stList_get(tokens, 14), "%lf", &(npRead->complementParams.scale_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement scale_sd\n");
    }
    j = sscanf(stList_get(tokens, 15), "%lf", &(npRead->complementParams.var_sd));
    if (j != 1) {
        st_errAbort("error parsing nanopore template var_sd\n");
    }
    j = sscanf(stList_get(tokens, 16), "%lf", &(npRead->complementParams.drift));
    if (j != 1) {
        st_errAbort("error parsing nanopore complement drift\n");
    }
    j = sscanf(stList_get(tokens, 17), "%i", &(npRead->twoD));
    if (j != 1) {
        st_errAbort("error parsing 'has twoD'\n");
    }

    // cleanup line 1
    free(string);
    stList_destruct(tokens);

    // line 2 [2D read sequence] (alignment table sequence) \n
    string = stFile_getLineFromFile(fH);
    if (npRead->twoD) {
        j = sscanf(string, "%s", npRead->twoDread);
        if (j != 1) {
            st_errAbort("error parsing read from npRead file, got %s\n", string);
        }
        npRead->twoDread[npRead->readLength] = '\0';
    } 
    free(string);

    // line 3 [template read] \n
    string = stFile_getLineFromFile(fH);
    j = sscanf(string, "%s", npRead->templateRead);
    if (j != 1) {
        st_errAbort("error parsing template read from npRead file got %s \n", string);
    }
    npRead->templateRead[npRead->templateReadLength] = '\0';
    free(string);

    // line 4 [template strand event map] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->templateReadLength) {
        st_errAbort(
                "template event map is not the correct length, should be %"PRId64", got %"PRId64"",
                npRead->templateReadLength,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < npRead->templateReadLength; i++) {
        j = sscanf(stList_get(tokens, i), "%"SCNd64, &(npRead->templateStrandEventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in template eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);
    // line 5 [complement read] \n
    string = stFile_getLineFromFile(fH);    
    if (npRead->twoD) {
        j = sscanf(string, "%s", npRead->complementRead);
        if (j != 1) {
            st_errAbort("error parsing complement read from npRead file got %s \n", string);
        }
        npRead->complementRead[npRead->complementReadLength] = '\0';
    }
    free(string);

    // line 6 [complement strand event map] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->complementReadLength) {
        st_errAbort(
                "template event map is not the correct length, should be %"PRId64", got %"PRId64"",
                npRead->complementReadLength,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < npRead->complementReadLength; i++) {
        j = sscanf(stList_get(tokens, i), "%"SCNd64, &(npRead->complementStrandEventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in template eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 7 [template 2D event map] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->readLength) {
        st_errAbort(
                "template event map is not the correct length, should be %lld, got %lld",
                npRead->readLength,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < npRead->readLength; i++) {
        j = sscanf(stList_get(tokens, i), "%"SCNd64, &(npRead->templateEventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in template eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 8 [template events (mean, stddev, length, start)] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbTemplateEvents * NB_EVENT_PARAMS)) {
        st_errAbort(
                "incorrect number of template events, should be %lld, got %lld",
                npRead->nbTemplateEvents,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < (npRead->nbTemplateEvents * NB_EVENT_PARAMS); i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->templateEvents[i]));
        if (j != 1) {
            st_errAbort("error loading in template events\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 9 [complement 2D event map] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check for correctness
    if (stList_length(tokens) != npRead->readLength) {
        st_errAbort(
                "complament event map is not the correct length, should be %lld, got %lld",
                npRead->readLength,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < npRead->readLength; i++) {
        j = sscanf(stList_get(tokens, i), "%"SCNd64, &(npRead->complementEventMap[i]));
        if (j != 1) {
            st_errAbort("error loading in complement eventMap\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 10 [complement events (mean, stddev, length)] \n
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbComplementEvents * NB_EVENT_PARAMS)) {
        st_errAbort(
                "incorrect number of complement events, should be %lld, got %lld",
                npRead->nbComplementEvents,
                stList_length(tokens));
    }
    for (int64_t i = 0; i < (npRead->nbComplementEvents * NB_EVENT_PARAMS); i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->complementEvents[i]));
        if (j != 1) {
            st_errAbort("error loading in complement events\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 11 model_state (template)
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    char *modelState;
    // check
    if (stList_length(tokens) != (npRead->nbTemplateEvents)) {
        st_errAbort("Got incorrect number of model states (kmers) got %"PRId64"\n");
    }
    for (int64_t i = 0; i < npRead->nbTemplateEvents; i++) {
        modelState = (char *)stList_get(tokens, i);
        npRead->templateModelState[i] = emissions_discrete_getKmerIndexFromPtr(modelState);
    }
    free(string);
    stList_destruct(tokens);

    // line 12 pModel (template)
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbTemplateEvents)) {
        st_errAbort("Got incorrect number of model probs got %"PRId64"\n");
    }
    for (int64_t i = 0; i < npRead->nbTemplateEvents; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->templatePModel[i]));
        if (j != 1) {
            st_errAbort("error loading in template model P(model_state)\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    // line 13 model_state (complement)
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbComplementEvents)) {
        st_errAbort("Got incorrect number of model states (kmers) got %"PRId64"\n");
    }
    for (int64_t i = 0; i < npRead->nbComplementEvents; i++) {
        modelState = (char *)stList_get(tokens, i);
        npRead->complementModelState[i] = emissions_discrete_getKmerIndexFromPtr(modelState);
    }
    free(string);
    stList_destruct(tokens);

    // line 14 pModel (complement)
    string = stFile_getLineFromFile(fH);
    tokens = stString_split(string);
    // check
    if (stList_length(tokens) != (npRead->nbComplementEvents)) {
        st_errAbort("Got incorrect number of model probs got %"PRId64"\n");
    }
    for (int64_t i = 0; i < npRead->nbComplementEvents; i++) {
        j = sscanf(stList_get(tokens, i), "%lf", &(npRead->complementPModel[i]));
        if (j != 1) {
            st_errAbort("error loading in template model P(model_state)\n");
        }
    }
    free(string);
    stList_destruct(tokens);

    fclose(fH);

    /*

    // needed
    npRead->readLength;
    npRead->templateReadLength;
    npRead->complementReadLength;
    npRead->nbTemplateEvents;
    npRead->nbComplementEvents;

    // needed
    npRead->twoD;
    npRead->templateEvents;
    npRead->templateStrandEventMap;
    npRead->templateEventMap;

    // needed if 2d
    npRead->complementEventMap;
    npRead->complementStrandEventMap;
    npRead->complementEvents;

    // needed if ESTIMATE_PARAMS
    npRead->templateRead;
    npRead->templateParams;
    npRead->complementRead;
    npRead->complementParams;

    // unused
    npRead->scaled;
    npRead->complementModelState;
    npRead->complementPModel;
    npRead->templateModelState;
    npRead->templatePModel;
    npRead->twoDread;

    // to verify:
    npRead->complementEventMap = NULL;
    npRead->complementStrandEventMap = NULL;
    npRead->complementEvents = NULL;
    npRead->complementModelState = NULL;
    npRead->complementPModel = NULL;
    npRead->templateModelState = NULL;
    npRead->templatePModel = NULL;
    npRead->scaled = NULL;
    npRead->twoDread = NULL;
    npRead->complementRead = NULL;

     */

    return npRead;
}

stList *nanopore_remapAnchorPairs(stList *anchorPairs, int64_t *eventMap) {
    stList *mappedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *pair = stList_get(anchorPairs, i);
        stList_append(mappedPairs,
                      stIntTuple_construct2(stIntTuple_get(pair, 0), eventMap[stIntTuple_get(pair, 1)]));
    }

    return mappedPairs;
}

stList *nanopore_remapAnchorPairsWithOffset(stList *unmappedPairs, int64_t *eventMap, int64_t mapOffset) {
    stList *mappedPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    for (int64_t i = 0; i < stList_length(unmappedPairs); i++) {
        stIntTuple *pair = stList_get(unmappedPairs, i);

        stList_append(mappedPairs,
                      stIntTuple_construct2(stIntTuple_get(pair, 0), eventMap[stIntTuple_get(pair, 1)] -
                              eventMap[mapOffset]));
    }

    return mappedPairs;
}

stList *nanopore_getAnchorKmersToEventsMap(stList *anchorPairs, double *eventSequence, char *nucleotideSequence) {
    stList *mapOfEventsToKmers = stList_construct3(0, &free);
    // loop over the remappedAnchors and make (event_mean, event_sd, reference_kmer_index) tuples
    for (int64_t i = 0; i < stList_length(anchorPairs); i++) {
        stIntTuple *pair = stList_get(anchorPairs, i);
        if (stIntTuple_length(pair) != 2) {
            st_errAbort("nanopore_getAnchorKmersToEventsMap: remapped tuple incorrect length\n");
        }
        int64_t ref_idx = stIntTuple_get(pair, 0);
        int64_t eventIndex = stIntTuple_get(pair, 1);
        int64_t kmerIndex = emissions_discrete_getKmerIndexFromPtr(nucleotideSequence + ref_idx);
        double eventMean = nanopore_getEventMean(eventSequence, eventIndex);
        double eventSd = nanopore_getEventSd(eventSequence, eventIndex);
        double eventDeltaTime = nanopore_getEventDeltaTime(eventSequence, eventIndex);
        EventKmerTuple *t = nanopore_eventKmerTupleConstruct(eventMean, eventSd, eventDeltaTime, kmerIndex);
        stList_append(mapOfEventsToKmers, t);
    }
    return mapOfEventsToKmers;
}

// This function basically reads off the 1D alignment table, and will give back eventKmerTuples that have kmers with
// whatever length the base caller used
static stList *nanopore_makeEventTuplesFromOneDRead(int64_t *kmerIndices, double *events, double *probs,
                                                    int64_t numberOfEvents, double threshold) {
    stList *map = stList_construct3(0, &free);
    for (int64_t i = 0; i < numberOfEvents; i++) {
        double p = probs[i];
        if (p < threshold) {
            continue;
        } else {
            int64_t kmerIndex = kmerIndices[i];
            double eventMean = nanopore_getEventMean(events, i);
            double eventSd = nanopore_getEventSd(events, i);
            double eventDuration = nanopore_getEventDuration(events, i);
            EventKmerTuple *t = nanopore_eventKmerTupleConstruct(eventMean, eventSd, eventDuration, kmerIndex);
            stList_append(map, t);
        }
    }
    return map;
}

stList *nanopore_getTemplateOneDAssignments(NanoporeRead *npRead, double threshold) {
    return nanopore_makeEventTuplesFromOneDRead(npRead->templateModelState, npRead->templateEvents,
                                                npRead->templatePModel, npRead->nbTemplateEvents, threshold);
}

stList *nanopore_getComplementOneDAssignments(NanoporeRead *npRead, double threshold) {
    return nanopore_makeEventTuplesFromOneDRead(npRead->complementModelState, npRead->complementEvents,
                                                npRead->complementPModel, npRead->nbComplementEvents, threshold);
}

stList *nanopore_getOneDAssignmentsFromRead(int64_t *strandEventMap, double *events, char *strandRead,
                                            int64_t readLength, char *alphabet, int64_t alphabetSize,
                                            int64_t kmerLength) {
    // map of event<>kmerIndex
    stList *map = stList_construct3(0, &free);

    // number of kmers there will be (this appears as rows in the .fast5 files, so I called it rows)
    int64_t rows = readLength - (kmerLength - 1);

    // see note below
    int64_t previousEventIndex = -1;
    for (int64_t i = 0; i < rows; i++) {
        int64_t eventIndex = strandEventMap[i];
        int64_t kmerIndex = kmer_id(strandRead + i, alphabet, alphabetSize, kmerLength);
        double eventMean = nanopore_getEventMean(events, eventIndex);
        double eventSd = nanopore_getEventSd(events, eventIndex);
        double eventDeltaTime = nanopore_getEventDeltaTime(events, eventIndex);

        // every event is mapped to it's highest probability position in the read so repeat events are just to fill in
        // skipped kmers (which will have lower prob). so only collect event<>kmer matches for the first event (highest
        // prob)
        if (eventIndex > previousEventIndex) {
            //st_uglyf("seqIndex: %lld eventIndex: %lld mean: %f sd: %f dur: %f kmerIndex: %lld\n",
            //         i, eventIndex, eventMean, eventSd, eventDeltaTime, kmerIndex);
            EventKmerTuple *t = nanopore_eventKmerTupleConstruct(eventMean, eventSd, eventDeltaTime, kmerIndex);
            stList_append(map, t);
            previousEventIndex = eventIndex;
        }
    }
    return map;
}

static void nanopore_adjustEventsForDriftP(double *events, double drift, int64_t nbEvents) {
    for (int64_t i = 0; i < nbEvents; i++) {
        double deltaTime = nanopore_getEventDeltaTime(events, i);
        events[i * NB_EVENT_PARAMS] = (events[i * NB_EVENT_PARAMS] - (deltaTime * drift));
    }
}

void nanopore_adjustEventsForDrift(NanoporeRead *npRead) {
    nanopore_adjustEventsForDriftP(npRead->templateEvents, npRead->templateParams.drift, npRead->nbTemplateEvents);
    nanopore_adjustEventsForDriftP(npRead->complementEvents, npRead->complementParams.drift,
                                   npRead->nbComplementEvents);
}

void nanopore_adjustTemplateEventsForDrift(NanoporeRead *npRead) {
    nanopore_adjustEventsForDriftP(npRead->templateEvents, npRead->templateParams.drift, npRead->nbTemplateEvents);
}

void nanopore_adjustComplementEventsForDrift(NanoporeRead *npRead) {
    nanopore_adjustEventsForDriftP(npRead->complementEvents, npRead->complementParams.drift,
                                   npRead->nbComplementEvents);
}

void nanopore_dontAdjustEvents(NanoporeRead *npRead) {
    (void) npRead;
    return;
}

void nanopore_descaleNanoporeRead(NanoporeRead *npRead) {
    nanopore_descaleEvents(npRead->nbTemplateEvents, npRead->templateEvents,
                           npRead->templateParams.scale,
                           npRead->templateParams.shift);
    nanopore_descaleEvents(npRead->nbComplementEvents, npRead->complementEvents,
                           npRead->complementParams.scale,
                           npRead->complementParams.shift);
    npRead->scaled = FALSE;
}

void nanopore_nanoporeReadDestruct(NanoporeRead *npRead) {
    free(npRead->twoDread);
    free(npRead->templateEventMap);
    free(npRead->templateEvents);
    free(npRead->complementEventMap);
    free(npRead->complementEvents);
    free(npRead->complementRead);
    free(npRead->templateRead);
    free(npRead->templateStrandEventMap);
    free(npRead->complementStrandEventMap);
    free(npRead->templateModelState);
    free(npRead->complementModelState);
    free(npRead);
}

double rand_uniform2(double a, double b) {
    return a + ((double) rand()) / ((double) RAND_MAX / (b - a));
}

// solve by Gaussian elimination
void nanopore_lineq_solve(double *A, double *b, double *x_out, int64_t n) {
    double* aux_matrix = (double*) malloc(sizeof(double) * n * n);

    int64_t idx;
    for (int64_t i = 0; i < n; i++) {
        x_out[i] = b[i];
        for (int64_t j = 0; j < n; j++) {
            idx = i * n + j;
            aux_matrix[idx] = A[idx];
        }
    }

    double factor;
    for (int64_t i = 0; i < n; i++) {
        if (fabs(aux_matrix[i * n + i]) < MACHEP) {
            int64_t swap = i + 1;
            while (aux_matrix[swap * n + i] < MACHEP) {
                swap++;
                if (swap >= n) {
                    fprintf(stderr, "Matrix is not invertible.\n");
                    exit(EXIT_FAILURE);
                }
            }
            double temp;
            for (int64_t j = 0; j < n; j++) {
                temp = aux_matrix[i * n + j];
                aux_matrix[i * n + j] = aux_matrix[swap * n + j];
                aux_matrix[swap * n + j] = temp;
            }
            temp = x_out[i];
            x_out[i] = x_out[swap];
            x_out[swap] = x_out[i];
        }

        factor = 1.0 / aux_matrix[i * n + i];
        x_out[i] *= factor;
        for (int64_t j = 0; j < n; j++) {
            aux_matrix[i * n + j] *= factor;
        }

        for (int64_t sub = i + 1; sub < n; sub++) {
            factor = aux_matrix[sub * n + i];
            x_out[sub] -= factor * x_out[i];
            for (int64_t j = 0; j < n; j++) {
                aux_matrix[sub * n + j] -= factor * aux_matrix[i * n + j];
            }
        }

    }

    for (int64_t i = n - 1; i >= 0; i--) {
        for (int64_t sub = i - 1; sub >= 0; sub--) {
            factor = aux_matrix[sub * n + i];
            x_out[sub] -= factor * x_out[i];
            for (int64_t j = 0; j < n; j++) {
                aux_matrix[sub * n + j] -= factor * aux_matrix[i * n + j];
            }
        }
    }

    free(aux_matrix);
}

// weighted least squares to compute normalization params
void nanopore_compute_mean_scale_params(double *model, stList *kmerToEventMap, NanoporeReadAdjustmentParameters *params,
                                        bool drift_out, bool var_out) {
    double* XWX = NULL;
    double* XWy = NULL;
    double* beta = NULL;
    if (stList_length(kmerToEventMap) == 0) {
        st_errAbort("Cannot get scale params with no assignments\n");
    }
    if (drift_out) {
        XWX = (double*) calloc(9, sizeof(double));
        XWy = (double*) calloc(3, sizeof(double));

        for (int i = 0; i < stList_length(kmerToEventMap); i++) {
            EventKmerTuple *t = stList_get(kmerToEventMap, i);
            double event = t->eventMean;
            double time = t->deltaTime;
            int64_t id = t->kmerIndex;

            double level_mean = model[(id * MODEL_PARAMS)];
            double level_sd = model[(id * MODEL_PARAMS + 1)];

            // weights should technically include variance parameter, but will only results in
            // some inefficiency, no bias
            double inv_var = 1.0 / (level_sd * level_sd);

            double scaled_mean = level_mean * inv_var;
            double scaled_time = time * inv_var;

            XWX[0] += inv_var;
            XWX[1] += scaled_mean;
            XWX[2] += scaled_time;
            XWX[4] += scaled_mean * level_mean;
            XWX[5] += scaled_mean * time;
            XWX[8] += scaled_time * time;

            XWy[0] += inv_var * event;
            XWy[1] += scaled_mean * event;
            XWy[2] += scaled_time * event;
        }
        XWX[3] = XWX[1];
        XWX[6] = XWX[2];
        XWX[7] = XWX[5];

        beta = (double*) malloc(sizeof(double) * 3);

        nanopore_lineq_solve(XWX, XWy, beta, 3);

        params->shift = beta[0];
        params->scale = beta[1];
        params->drift = beta[2];

        if (var_out) {
            double dispersion = 0.0;
            for (int64_t i = 0; i < stList_length(kmerToEventMap); i++) {
                EventKmerTuple *t = stList_get(kmerToEventMap, i);

                int64_t id = t->kmerIndex;

                double level_mean = model[(id * MODEL_PARAMS)];
                double level_sd = model[(id * MODEL_PARAMS + 1)];
                double level_var = level_sd * level_sd;

                double event = t->eventMean;
                double time = t->deltaTime;

                double predicted_val = beta[0] + beta[1] * level_mean + beta[2] * time;
                double residual = event - predicted_val;
                dispersion += (residual * residual) / level_var;
            }
            params->var = sqrt(dispersion / stList_length(kmerToEventMap));
        }
    }
    else {
        XWX = (double*) calloc(4, sizeof(double));
        XWy = (double*) calloc(2, sizeof(double));

        for (int i = 0; i < stList_length(kmerToEventMap); i++) {
            EventKmerTuple *t = stList_get(kmerToEventMap, i);

            double event = t->eventMean;
            int64_t id = t->kmerIndex;

            double level_mean = model[(id * MODEL_PARAMS)];

            double level_sd = model[(id * MODEL_PARAMS + 1)];
            // weights should technically include variance parameter, but will only results in
            // some inefficiency, no bias
            double inv_var = 1.0 / (level_sd * level_sd);

            double scaled_mean = level_mean * inv_var;

            XWX[0] += inv_var;
            XWX[1] += scaled_mean;
            XWX[3] += scaled_mean * level_mean;

            XWy[0] += inv_var * event;
            XWy[1] += scaled_mean * event;
        }
        XWX[2] = XWX[1];

        beta = (double*) malloc(sizeof(double) * 2);

        nanopore_lineq_solve(XWX, XWy, beta, 2);

        params->shift = beta[0];
        params->scale = beta[1];

        if (var_out) {
            double dispersion = 0.0;
            for (int64_t i = 0; i < stList_length(kmerToEventMap); i++) {
                EventKmerTuple *t = stList_get(kmerToEventMap, i);

                int64_t id = t->kmerIndex;

                double level_mean = model[(id * MODEL_PARAMS)];

                double level_sd = model[(id * MODEL_PARAMS + 1)];
                double level_var = level_sd * level_sd;

                double event = t->eventMean;

                double predicted_val = beta[0] + beta[1] * level_mean;
                double residual = event - predicted_val;
                dispersion += (residual * residual) / level_var;
            }
            params->var = sqrt(dispersion / stList_length(kmerToEventMap));
        }
    }
    free(XWX);
    free(XWy);
    free(beta);
}

void nanopore_compute_noise_scale_params(double *model, stList *kmerToEventMap,
                                         NanoporeReadAdjustmentParameters *params) {
    double* XWX = NULL;
    double* XWy = NULL;
    double* beta = NULL;
    if (stList_length(kmerToEventMap) == 0) {
        st_errAbort("Cannot get scale params with no assignments\n");
    }

    XWX = (double*) calloc(4, sizeof(double));
    XWy = (double*) calloc(2, sizeof(double));

    for (int i = 0; i < stList_length(kmerToEventMap); i++) {
        EventKmerTuple *t = stList_get(kmerToEventMap, i);

        double noise = t->eventSd;
        int64_t id = t->kmerIndex;

        double noise_mean = model[(id * MODEL_PARAMS + 2)];

        double noise_sd = model[(id * MODEL_PARAMS + 3)];
        // weights should technically include variance parameter, but will only results in
        // some inefficiency, no bias
        double inv_var = 1.0 / (noise_sd * noise_sd);

        double scaled_mean = noise_mean * inv_var;

        XWX[0] += inv_var;
        XWX[1] += scaled_mean;
        XWX[3] += scaled_mean * noise_mean;

        XWy[0] += inv_var * noise;
        XWy[1] += scaled_mean * noise;
    }
    XWX[2] = XWX[1];

    beta = (double*) malloc(sizeof(double) * 2);

    nanopore_lineq_solve(XWX, XWy, beta, 2);

    params->shift_sd = beta[0];
    params->scale_sd = beta[1];

    double dispersion = 0.0;
    for (int64_t i = 0; i < stList_length(kmerToEventMap); i++) {
        EventKmerTuple *t = stList_get(kmerToEventMap, i);

        int64_t id = t->kmerIndex;

        double nosie_mean = model[(id * MODEL_PARAMS + 2)];
        double noise_sd = model[(id * MODEL_PARAMS + 3)];

        double level_var = noise_sd * noise_sd;

        double event = t->eventSd;

        double predicted_val = beta[0] + beta[1] * nosie_mean;
        double residual = event - predicted_val;
        dispersion += (residual * residual) / level_var;
    }
    params->var_sd = sqrt(dispersion / stList_length(kmerToEventMap));

    free(XWX);
    free(XWy);
    free(beta);
}

void nanopore_convert_to_lognormal_params(int64_t alphabet_size, int64_t kmer_length, double *model,
                                          stList *kmerToEventMap) {
    int64_t num_kmers = 1;
    for (int64_t i = 0; i < alphabet_size; i++) {
        num_kmers *= kmer_length;
    }

    for (int64_t i = 0; i < num_kmers; i++) {
        double noise_sd = model[(i * MODEL_PARAMS + 3)];
        double untrans_var = noise_sd * noise_sd;
        double untrans_mean = model[(i * MODEL_PARAMS + 2)];

        double trans_var = log(1.0 + untrans_var / (untrans_mean * untrans_mean));
        double trans_mean = log(untrans_mean - trans_var / 2.0);

        //noise_sds[i] = sqrt(trans_var);
        model[(i * MODEL_PARAMS + 3)] = sqrt(trans_var);

        //noise_means[i] = trans_mean;
        model[(i * MODEL_PARAMS + 2)] = trans_mean;
    }

    for (int i = 0; i < stList_length(kmerToEventMap); i ++) {
        EventKmerTuple *t = stList_get(kmerToEventMap, i);
        t->eventSd = log(t->eventSd);
        //read_noises[i] = log(read_noises[i]);
    }
}
