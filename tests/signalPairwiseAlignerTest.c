#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <inttypes.h>
#include <stdbool.h>
#include <nanopore.h>
#include "stateMachine.h"
#include "CuTest.h"
#include "sonLib.h"
#include "pairwiseAligner.h"
#include "randomSequences.h"


// brute force probability formulae
static double test_standardNormalPdf(double x) {
    double pi = 3.141592653589793;
    double inv_sqTwoPi = 1 / (sqrt(2*pi));
    double result = inv_sqTwoPi * exp(-(x * x) / 2);
    return result;
}

static double test_normalPdf(double x, double mu, double sigma) {
    double pi = 3.141592653589793;
    double inv_sqTwoPi = 1 / (sqrt(2*pi));
    double c = inv_sqTwoPi * (1/sigma);
    double a = (x - mu) / sigma;
    double result = c * exp(-0.5f * a * a);
    return result;
}

static double test_inverseGaussianPdf(double x, double mu, double lambda) {
    double pi = 3.141592653589793;
    double c = lambda / (2 * pi * pow(x, 3.0));
    double sqt_c = sqrt(c);
    double xmmu = x - mu;
    double result = sqt_c * exp((-lambda * xmmu * xmmu)/
                                (2 * mu * mu * x));
    return result;
}

Sequence *makeTestKmerSequence() {
    char *s = "ATGXAXA"; // has 2 6mers
    int64_t lX = sequence_correctSeqLength(strlen(s), kmer, KMER_LENGTH);
    Sequence *seq = sequence_construct(lX, s, sequence_getKmer, kmer);
    return seq;
}

Sequence *makeKmerSequence(char *nucleotides) {
    int64_t lX = sequence_correctSeqLength(strlen(nucleotides), kmer, KMER_LENGTH);
    Sequence *seq = sequence_constructKmerSequence(lX, nucleotides,
                                                   sequence_getKmer, sequence_sliceNucleotideSequence,
                                                   kmer);
    return seq;
}

//////////////////////////////////////////////Function Tests/////////////////////////////////////////////////////////
/*  Old test that I'm keeping around in case I re-instate the duration model
static void test_poissonPosteriorProb(CuTest *testCase) {
    double event1[] = {62.784241, 0.664989, 0.00332005312085};
    double test0 = emissions_signal_getDurationProb(event1, 0);
    double test1 = emissions_signal_getDurationProb(event1, 1);
    double test2 = emissions_signal_getDurationProb(event1, 2);
    double test3 = emissions_signal_getDurationProb(event1, 3);
    double test4 = emissions_signal_getDurationProb(event1, 4);
    double test5 = emissions_signal_getDurationProb(event1, 5);
    CuAssertTrue(testCase, test0 < test1);
    CuAssertTrue(testCase, test1 > test2);
    CuAssertTrue(testCase, test2 > test3);
    CuAssertTrue(testCase, test3 > test4);
    CuAssertTrue(testCase, test4 > test5);
    //st_uglyf("0 - %f\n1 - %f\n2 - %f\n3 - %f\n4 - %f\n5 - %f\n", test0, test1, test2, test3, test4, test5);
}
*/

static void test_stateMachine3EmissionsPdfs(CuTest *testCase) {
    // test Gaussian distribution PDF
    double control = test_standardNormalPdf(0);
    double obs = emissions_signal_logGaussianProbabilityDensity(0, 0, 1);
    CuAssertDblEquals(testCase, log(control), obs, 0.001);
    CuAssertDblEquals(testCase, control, exp(obs), 0.001);
    st_randomSeed(5);
    double sigma = 0.0;
    while (sigma <= 0.1) {
        sigma = st_random();
    }
    double mu = rand_uniform2(0, 10);
    double posOrNeg = st_random();
    double x = 0.0;

    if (posOrNeg > 0.5) {
        x = mu - st_random();
    } else {
        x = mu + st_random();
    }
    control = test_normalPdf(x, mu, sigma);
    obs = emissions_signal_logGaussianProbabilityDensity(x, mu, sigma);
    CuAssertDblEquals(testCase, log(control), obs, 0.001);
    CuAssertDblEquals(testCase, control, exp(obs), 0.001);
    // test Inverse Gaussian PDF
    double lambda = st_random();
    control = test_inverseGaussianPdf(x, mu, lambda);
    obs = emissions_signal_logInverseGaussianProbabilityDensity(x, mu, lambda);
    CuAssertDblEquals(testCase, log(control), obs, 0.001);
    CuAssertDblEquals(testCase, control, exp(obs), 0.001);
}

static bool testDiagonalsEqual(Diagonal d1, Diagonal d2) {
    bool b = diagonal_equals(d1, d2);
    if (!b) {
        st_logCritical("Diagonals not equal: d1: %s, d2: %s \n", diagonal_getString(d1), diagonal_getString(d2));
    }
    return b;
}

static void test_logAdd(CuTest *testCase) {
    for (int64_t test = 0; test < 100000; test++) {
        double i = st_random();
        double j = st_random();
        double k = i + j;
        double l = exp(logAdd(log(i), log(j)));
        //st_logInfo("I got %f %f\n", k, l);
        CuAssertTrue(testCase, l < k + 0.001);
        CuAssertTrue(testCase, l > k - 0.001);
    }
}

static void test_genericSequenceTests(CuTest *testCase, Sequence *testSequence, int64_t length,
                                      char *tS) {
    CuAssertIntEquals(testCase, length, testSequence->length);
    CuAssertStrEquals(testCase, tS, testSequence->elements);
    for (int64_t i = 0; i < length; i++) {
        CuAssert(testCase, "Nucleotide match fail",
                 *(tS + i) == *(char *)(testSequence->get(testSequence->elements, i)));
    }
}


static void test_Sequence(CuTest *testCase) {
    int64_t length = 1000;
    char *tS = getRandomSequence(length);
    Sequence* testSequence = sequence_construct(length, tS, sequence_getKmer, nucleotide);
    test_genericSequenceTests(testCase, testSequence, length, tS);
    testSequence = sequence_construct2(length, tS, sequence_getKmer, sequence_sliceNucleotideSequence,
                                       nucleotide);
    sequence_destruct(testSequence);
}

static void test_referenceSequence(CuTest *testCase) {
    int64_t length = 1000;
    char *tS = getRandomSequence(length);

    // test construct Kmer sequence
    Sequence *testSequence = sequence_constructKmerSequence(length, tS, sequence_getKmer,
                                                            sequence_sliceNucleotideSequence,
                                                            kmer);
    test_genericSequenceTests(testCase, testSequence, length, tS);
//    test_referenceSequenceTests(testCase, testSequence);

    // test copy
    Sequence *copy = sequence_deepCopyNucleotideSequence(testSequence);
    test_genericSequenceTests(testCase, copy, length, tS);
//    test_referenceSequenceTests(testCase, copy);
    CuAssertStrEquals(testCase, testSequence->elements, copy->elements);
    free(copy);

    // test slicing
    int64_t r = st_randomInt(0, length);
    Sequence *slice = testSequence->sliceFcn(testSequence, 10, length - r);
    CuAssertStrEquals(testCase, ((char *)testSequence->elements + 10), slice->elements);
    CuAssert(testCase, "slice sequence type fail",testSequence->type == slice->type);

    sequence_deepDestruct(testSequence);
    sequence_destruct(slice);
}

static void test_eventSequence(CuTest *testCase) {
    double sY[24] = {
            58.743435, 0.887833, 0.0571, //ACGATA 0
            53.604965, 0.816836, 0.0571, //CGATAC 1
            58.432015, 0.735143, 0.0571, //GATACG 2
            63.684352, 0.795437, 0.0571, //ATACGG 3
            63.520262, 0.757803, 0.0571, //TACGGA 4
            58.921430, 0.812959, 0.0571, //ACGGAC 5
            59.895882, 0.740952, 0.0571, //CGGACA 6
            61.684303, 0.722332, 0.0571, //GGACAT 7
    };

    Sequence *eventSequence = sequence_constructEventSequence(8, sY);
    for (int64_t i = 0; i < eventSequence->length; i++) {
        int64_t index = i * NB_EVENT_PARAMS;
        double actual = *(double *)eventSequence->get(eventSequence->elements, i);
        CuAssertDblEquals(testCase, sY[index], actual, 0.0);
    }

    int64_t newStart = 1;
    int64_t sliceLength = 3;
    Sequence *slice = eventSequence->sliceFcn(eventSequence, newStart, sliceLength);
    for (int64_t i = 0; i < slice->length; i++) {
        int64_t index = (i + newStart) * NB_EVENT_PARAMS;
        double actual = *(double *)slice->get(slice->elements, i);
        CuAssertDblEquals(testCase, sY[index], actual, 0.0);
    }

    sequence_destruct(eventSequence);
    sequence_destruct(slice);

}

static void test_1dNanoporeRead(CuTest *testCase) {
    char *tempFile = stString_print("../tests/test_npReads/r9p4_oneD.npRead");
    CuAssertTrue(testCase, stFile_exists(tempFile));

    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(tempFile, "ACGT", 4, 5);
    CuAssertTrue(testCase, npRead->twoD == FALSE);

}

static void test_loadNanoporeRead(CuTest *testCase) {
    int64_t length = 4096;
    char *read = getRandomSequence(length);
    double param = 0.0;
    double prob = 1.0;
    char *tempFile = stString_print("./tempRead.npread");
    stList *kmers = path_listPotentialKmers(KMER_LENGTH, strlen(CANONICAL_NUCLEOTIDES), CANONICAL_NUCLEOTIDES);
    int twodFlag = 1;

    CuAssertTrue(testCase, !stFile_exists(tempFile));
    FILE *fH = fopen(tempFile, "w");

    // write line 1
    fprintf(fH, "%"PRId64"\t""%"PRId64"\t""%"PRId64"\t""%"PRId64"\t""%"PRId64""
                    "\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%f\t%i\n",
            length, length, length, length, length,
            param, param, param, param, param, param, param, param, param, param, param, param, twodFlag);

    // line 2 2D (alignment table) read sequence
    fprintf(fH, "%s\n", read);

    // line 3 template read
    fprintf(fH, "%s\n", read);

    // line 4 template 'strand' map
    for (int64_t i = 0; i < length; i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");

    // line 5 complement read
    fprintf(fH, "%s\n", read);

    // line 6 complement 'stand' map
    for (int64_t i = 0; i < length; i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");

    // line 7 template 2D 'event' map
    for (int64_t i = 0; i < length; i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");

    // line 8 template events
    for (int64_t i = 0; i < (length * NB_EVENT_PARAMS); i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");

    // line 9 complement 2D 'event' map
    for (int64_t i = 0; i < length; i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");

    // line 10 complement events
    for (int64_t i = 0; i < (length * NB_EVENT_PARAMS); i++) {
        fprintf(fH, "%"PRId64"\t", i);
    }
    fprintf(fH, "\n");

    // line 11 template model_state
    char *kmer;
    for (int64_t i = 0; i < length; i++) {
        kmer = (char *)stList_get(kmers, i);
        fprintf(fH, "%s\t", kmer);
    }
    fprintf(fH, "\n");

    // line 12 pModel template
    for (int64_t i = 0; i < length; i++) {
        fprintf(fH, "%f\t", prob);
    }
    fprintf(fH, "\n");

    // line 13 complement model_state
    for (int64_t i = 0; i < length; i++) {
        kmer = (char *)stList_get(kmers, i);
        fprintf(fH, "%s\t", kmer);
    }
    fprintf(fH, "\n");

    // line 14 pModel complement
    for (int64_t i = 0; i < length; i++) {
        fprintf(fH, "%f\t", prob);
    }
    fprintf(fH, "\n");

    fclose(fH);

    NanoporeRead *npRead = nanopore_loadNanoporeReadFromFile(tempFile, "ACGT", 4, KMER_LENGTH);
    CuAssertTrue(testCase, npRead->readLength == length);
    CuAssertTrue(testCase, npRead->templateReadLength == length);
    CuAssertTrue(testCase, npRead->complementReadLength == length);
    CuAssertTrue(testCase, npRead->nbTemplateEvents == length);
    CuAssertTrue(testCase, npRead->nbComplementEvents == length);

    CuAssertTrue(testCase, npRead->templateParams.scale == param);
    CuAssertTrue(testCase, npRead->templateParams.shift == param);
    CuAssertTrue(testCase, npRead->templateParams.var == param);
    CuAssertTrue(testCase, npRead->templateParams.scale_sd == param);
    CuAssertTrue(testCase, npRead->templateParams.var_sd == param);

    CuAssertTrue(testCase, npRead->complementParams.scale == param);
    CuAssertTrue(testCase, npRead->complementParams.shift == param);
    CuAssertTrue(testCase, npRead->complementParams.var == param);
    CuAssertTrue(testCase, npRead->complementParams.scale_sd == param);
    CuAssertTrue(testCase, npRead->complementParams.var_sd == param);
    CuAssertTrue(testCase, npRead->twoD);

    CuAssertStrEquals(testCase, npRead->twoDread, read);
    CuAssertStrEquals(testCase, npRead->templateRead, read);
    CuAssertStrEquals(testCase, npRead->complementRead, read);

    for (int64_t i = 0; i < length; i++) {
        CuAssertTrue(testCase, npRead->templateEventMap[i] == i);
        CuAssertTrue(testCase, npRead->complementEventMap[i] == i);
        CuAssertTrue(testCase, npRead->templateStrandEventMap[i] == i);
        CuAssertTrue(testCase, npRead->complementStrandEventMap[i] == i);
    }
    for (int64_t i = 0; i < (length * NB_EVENT_PARAMS); i++) {
        CuAssertTrue(testCase, npRead->templateEvents[i] == i);
        CuAssertTrue(testCase, npRead->complementEvents[i] == i);
    }

    for (int64_t i = 0; i < length; i++) {
        kmer = (char *)stList_get(kmers, i);
        int64_t index = kmer_id(kmer, CANONICAL_NUCLEOTIDES, 4, KMER_LENGTH);
        CuAssertIntEquals(testCase, index, npRead->templateModelState[i]);
        CuAssertIntEquals(testCase, index, npRead->complementModelState[i]);
    }
    for (int64_t i = 0; i < length; i++) {
        CuAssertDblEquals(testCase, npRead->templatePModel[i], prob, 0.0);
        CuAssertDblEquals(testCase, npRead->complementPModel[i], prob, 0.0);
    }

    nanopore_nanoporeReadDestruct(npRead);
    stFile_rmrf(tempFile);
    stList_destruct(kmers);
    free(read);
    free(tempFile);
}

static void test_getSplitPoints(CuTest *testCase) {
    int64_t matrixSize = 2000 * 2000;

    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);

    //Test a small region, which produces no splits
    int64_t lX = 3000;
    int64_t lY = 1000;
    stList *splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, lX, lY)));
    stList_destruct(splitPoints);

    //Test with one really big matrix with no anchors
    lX = 20000;
    lY = 25000;
    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 1, 1);
    CuAssertIntEquals(testCase, 0, stList_length(splitPoints));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 1, 0);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase,
                 stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(18000, 23000, lX, lY)));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 1);
    CuAssertIntEquals(testCase, 1, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 2000, 2000)));
    stList_destruct(splitPoints);

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);
    CuAssertIntEquals(testCase, 2, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 2000, 2000)));
    CuAssertTrue(testCase,
                 stIntTuple_equalsFn(stList_get(splitPoints, 1), stIntTuple_construct4(18000, 23000, lX, lY)));
    stList_destruct(splitPoints);

    //Now test with some more points
    stList_append(anchorPairs, stIntTuple_construct2(2000, 2000)); //This should not create a split
    stList_append(anchorPairs, stIntTuple_construct2(4002, 4001)); //This should cause a split
    stList_append(anchorPairs, stIntTuple_construct2(5000, 5000)); //This should not cause a split
    stList_append(anchorPairs, stIntTuple_construct2(8000, 6000)); //Neither should this (it is maximum sized)
    stList_append(anchorPairs, stIntTuple_construct2(9000, 9000)); //Or this (it is maximum sized)
    stList_append(anchorPairs, stIntTuple_construct2(10000, 14000)); //This should create a split
    stList_append(anchorPairs, stIntTuple_construct2(15000, 15000)); //This should also create a split
    stList_append(anchorPairs, stIntTuple_construct2(16000, 16000)); //This should not, but there will be a split with the end.

    splitPoints = getSplitPoints(anchorPairs, lX, lY, matrixSize, 0, 0);

    for (int64_t i = 0; i < stList_length(splitPoints); i++) {
        stIntTuple *j = stList_get(splitPoints, i);
        st_logInfo("I got split point: x1: %" PRIi64 " y1: %" PRIi64 " x2: %" PRIi64 " y2: %" PRIi64 "\n",
                   stIntTuple_get(j, 0), stIntTuple_get(j, 1), stIntTuple_get(j, 2), stIntTuple_get(j, 3));
    }

    CuAssertIntEquals(testCase, 5, stList_length(splitPoints));
    CuAssertTrue(testCase, stIntTuple_equalsFn(stList_get(splitPoints, 0), stIntTuple_construct4(0, 0, 3001, 3001)));
    CuAssertTrue(testCase,
                 stIntTuple_equalsFn(stList_get(splitPoints, 1), stIntTuple_construct4(3002, 3001, 9500, 11001)));
    CuAssertTrue(testCase,
                 stIntTuple_equalsFn(stList_get(splitPoints, 2), stIntTuple_construct4(9501, 12000, 12001, 14500)));
    CuAssertTrue(testCase,
                 stIntTuple_equalsFn(stList_get(splitPoints, 3), stIntTuple_construct4(13000, 14501, 18000, 18001)));
    CuAssertTrue(testCase,
                 stIntTuple_equalsFn(stList_get(splitPoints, 4), stIntTuple_construct4(18001, 23000, 20000, 25000)));

    stList_destruct(splitPoints);
    stList_destruct(anchorPairs);
}

static void test_bands(CuTest *testCase) {
    stList *anchorPairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
    ///stList_append(anchorPairs, stIntTuple_construct2( 0, 0));
    stList_append(anchorPairs, stIntTuple_construct2(1, 0));
    stList_append(anchorPairs, stIntTuple_construct2(2, 1));
    stList_append(anchorPairs, stIntTuple_construct2(3, 3));
    /////stList_append(anchorPairs, stIntTuple_construct2( 5, 4));
    //Start the traversal
    int64_t lX = 6, lY = 5;
    Band *band = band_construct(anchorPairs, lX, lY, 2);
    BandIterator *bandIt = bandIterator_construct(band);

    //Forward pass
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(2, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(3, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(4, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(5, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(6, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(11, 1, 1)));

    //Go backward
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(11, 1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(7, -3, 3)));
    //Now walk forward a bit
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(10, 0, 2)));
    //Now carry on back again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(10, 0, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(9, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(8, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(7, -3, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(6, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(5, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(4, -2, 4)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(3, -1, 3)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(2, -2, 2)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));

    //Now forward again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(0, 0, 0)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getNext(bandIt), diagonal_construct(1, -1, 1)));
    //Now back again
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(1, -1, 1)));
    CuAssertTrue(testCase, testDiagonalsEqual(bandIterator_getPrevious(bandIt), diagonal_construct(0, 0, 0)));

    //Cleanup
    bandIterator_destruct(bandIt);
    band_destruct(band);
    stList_destruct(anchorPairs);
}

static void test_diagonal(CuTest *testCase) {
    //Construct an example diagonal.
    int64_t xL = 10, yL = 20, xU = 30, yU = 0; //Coordinates of the upper and lower
    //pairs in x,y coordinates
    Diagonal d = diagonal_construct(xL + yL, xL - yL, xU - yU);
    CuAssertIntEquals(testCase, diagonal_getXay(d), xL + yL);
    CuAssertIntEquals(testCase, diagonal_getMinXmy(d), xL - yL);
    CuAssertIntEquals(testCase, diagonal_getMaxXmy(d), xU - yU);
    CuAssertIntEquals(testCase, diagonal_getWidth(d), (xU - yU - (xL - yL)) / 2 + 1);
    CuAssertIntEquals(testCase, diagonal_getXCoordinate(xL + yL, xL - yL), xL);
    CuAssertIntEquals(testCase, diagonal_getYCoordinate(xL + yL, xL - yL), yL);
    CuAssertTrue(testCase, diagonal_equals(d, d));
    CuAssertTrue(testCase, !diagonal_equals(d, diagonal_construct(0, 0, 0)));
    //A bogus diagonal is one such that |xay + xmy| % 2 != 0 or such that xmyR < xmyL.
    //Try constructing bogus diagonals, should throw exceptions.
    stTry
        {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry
        {
            diagonal_construct(10, 5, 5);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
    stTry
        {
            diagonal_construct(10, 6, 4);
            CuAssertTrue(testCase, 0);
        }
        stCatch(PAIRWISE_ALIGNMENT_EXCEPTION_ID)
            {
                st_logInfo(stExcept_getMsg(PAIRWISE_ALIGNMENT_EXCEPTION_ID));
            }stTryEnd
}

static void test_hdCellConstruct(CuTest *testCase) {
    char *ambigKmer = "ATGXAXAAAAAA";
    int64_t nbCytosines = 3;
    char *cytosines = "CEO";
    HDCell *cell = hdCell_construct(ambigKmer, 3, nbCytosines, cytosines, 6);
    Path *path = hdCell_getPath(cell, 0);
    Path *path2 = hdCell_getPath(cell, 8);
    CuAssertTrue(testCase, hdCell_getPath(cell, 9) == NULL);
    CuAssertIntEquals(testCase, 9, (int )cell->numberOfPaths);
    CuAssertStrEquals(testCase, path->kmer, "ATGCAC");
    CuAssertStrEquals(testCase, path2->kmer, "ATGOAO");
    hdCell_destruct(cell);
}

static void test_hdCellConstructWorstCase(CuTest *testCase) {
    char *ambigKmer = "XXXXXX";
    int64_t nbCytosines = 3;
    char *cytosines = "CEO";
    HDCell *cell = hdCell_construct(ambigKmer, 3, nbCytosines, cytosines, 6);
    Path *path = hdCell_getPath(cell, 0);
    Path *path2 = hdCell_getPath(cell, 728);
    CuAssertIntEquals(testCase, 729, (int )cell->numberOfPaths);
    CuAssertStrEquals(testCase, path->kmer, "CCCCCC");
    CuAssertStrEquals(testCase, path2->kmer, "OOOOOO");
    hdCell_destruct(cell);
}

static void test_dpDiagonal(CuTest *testCase) {
    // load model and make stateMachine
    char *testModelPath = stString_print("../../signalAlign/models/testModelR73_acegot_template.model");
    StateMachine *sM = getStateMachine3(testModelPath);

    Diagonal diagonal = diagonal_construct(3, -1, 1); // makes a diagonal with 2 cells

    Sequence *seq = makeTestKmerSequence();  // ATGXAXA

    DpDiagonal *dpDiagonal = dpDiagonal_construct(diagonal, sM->stateNumber, sM->kmerLength, seq);

    //Get cell
    HDCell *c1 = dpDiagonal_getCell(dpDiagonal, -1);
    CuAssertTrue(testCase, c1 != NULL);

    HDCell *c2 = dpDiagonal_getCell(dpDiagonal, 1); // gets cell 1
    CuAssertTrue(testCase, c2 != NULL);

    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, 3) == NULL);
    CuAssertTrue(testCase, dpDiagonal_getCell(dpDiagonal, -3) == NULL);

    dpDiagonal_initialiseValues(dpDiagonal, sM, sM->endStateProb); //Test initialise values
    double totalProb = LOG_ZERO;

    for (int64_t p = 0; p < c1->numberOfPaths; p++) {
        Path *path1 = hdCell_getPath(c1, p);
        Path *path2 = hdCell_getPath(c2, p);
        for (int64_t s = 0; s < path1->stateNumber; s++) {
            CuAssertDblEquals(testCase, path1->cells[s], sM->endStateProb(sM, s), 0.0);
            CuAssertDblEquals(testCase, path2->cells[s], sM->endStateProb(sM, s), 0.0);
            totalProb = logAdd(totalProb, 2 * path1->cells[s]);
            totalProb = logAdd(totalProb, 2 * path2->cells[s]);
        }
    }

    DpDiagonal *dpDiagonal2 = dpDiagonal_clone(dpDiagonal);
    CuAssertTrue(testCase, dpDiagonal_equals(dpDiagonal, dpDiagonal2));
    //Check it runs
    CuAssertDblEquals(testCase, totalProb, dpDiagonal_dotProduct(dpDiagonal, dpDiagonal2), 0.001);
    free(testModelPath);
    dpDiagonal_destruct(dpDiagonal);
    dpDiagonal_destruct(dpDiagonal2);

}

static void test_dpMatrix(CuTest *testCase) {
    int64_t lX = 3, lY = 2;
    char *s = "ATGXAATT";
    //         ATGCAA   0
    //          TGCAAT  1
    //           GCAATT 2
    Sequence *sX = makeKmerSequence(s);

    DpMatrix *dpMatrix = dpMatrix_construct(lX + lY, 5, KMER_LENGTH);

    // check initialization
    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    // make sure there aren't any fantom diagonals
    for (int64_t i = -1; i <= lX + lY + 10; i++) {
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
    }

    // make some diagonals in the dpMatrix, and check them, then make sure that
    // the number of active diagonals is correct.
    for (int64_t i = 0; i <= lX + lY; i++) {
        DpDiagonal *dpDiagonal = dpMatrix_createDiagonal(dpMatrix, diagonal_construct(i, -i, i), sX);
        CuAssertTrue(testCase, dpDiagonal == dpMatrix_getDiagonal(dpMatrix, i));
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i + 1);
    }

    // test for destroying diagonals
    for (int64_t i = lX + lY; i >= 0; i--) {
        dpMatrix_deleteDiagonal(dpMatrix, i);
        CuAssertTrue(testCase, dpMatrix_getDiagonal(dpMatrix, i) == NULL);
        CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), i);
    }

    // double check that they are gone
    CuAssertIntEquals(testCase, dpMatrix_getActiveDiagonalNumber(dpMatrix), 0);

    dpMatrix_destruct(dpMatrix);
}

static void checkBlastPairs(CuTest *testCase, stList *blastPairs, int64_t lX, int64_t lY, bool checkNonOverlapping) {
    //st_uglyf("I got %" PRIi64 " pairs to check\n", stList_length(blastPairs));
    int64_t pX = -1;
    int64_t pY = -1;
    for (int64_t i = 0; i < stList_length(blastPairs); i++) {
        stIntTuple *j = stList_get(blastPairs, i);
        CuAssertTrue(testCase, stIntTuple_length(j) == 2);

        int64_t x = stIntTuple_get(j, 0);
        int64_t y = stIntTuple_get(j, 1);

        CuAssertTrue(testCase, x >= 0);
        CuAssertTrue(testCase, y >= 0);
        CuAssertTrue(testCase, x < lX);
        CuAssertTrue(testCase, y < lY);
        if (checkNonOverlapping) {
            CuAssertTrue(testCase, x > pX);
            CuAssertTrue(testCase, y > pY);
        }
        pX = x;
        pY = y;
    }
}

static void test_getBlastPairs(CuTest *testCase) {
    /*
     * Test the blast heuristic to get the different pairs.
     */
    for (int64_t test = 0; test < 10; test++) {
        //Make a pair of sequences
        char *sX = getRandomSequence(st_randomInt(0, 10000));
        char *sY = evolveSequence(sX); //stString_copy(seqX);
        int64_t lX = strlen(sX), lY = strlen(sY);
        st_logInfo("Sequence X to align: %s END, seq length %" PRIi64 "\n", sX, lX);
        st_logInfo("Sequence Y to align: %s END, seq length %" PRIi64 "\n", sY, lY);
        int64_t trim = st_randomInt(0, 5);
        bool repeatMask = st_random() > 0.5;
        st_logInfo("Using random trim %" PRIi64 ", recursive %" PRIi64 " \n", trim, repeatMask);
        stList *blastPairs = getBlastPairs(sX, sY, trim, repeatMask);
        checkBlastPairs(testCase, blastPairs, lX, lY, 0);
        stList_destruct(blastPairs);
        free(sX);
        free(sY);
    }
}
/*
static void test_filterToRemoveOverlap(CuTest *testCase) {
    for (int64_t i = 0; i < 5; i++) {
        //Make random pairs
        int64_t lX = st_randomInt(0, 1000);
        int64_t lY = st_randomInt(0, 1000);
        stList *pairs = stList_construct3(0, (void (*)(void *)) stIntTuple_destruct);
        double acceptProb = st_random() + 0.1;
        for (int64_t x = 0; x < lX; x++) {
            for (int64_t y = 0; y < lY; y++) {
                if (st_random() > acceptProb) {
                    stList_append(pairs, stIntTuple_construct2(x, y));
                }
            }
        }

        //Now run filter pairs
        stList *nonoverlappingPairs = filterToRemoveOverlap(pairs);

        //Check non overlapping
        checkBlastPairs(testCase, nonoverlappingPairs, lX, lY, 1);

        //Now check maximal
        stList *nonoverlappingPairs2 = stList_construct();
        for (int64_t i = 0; i < stList_length(pairs); i++) {
            stIntTuple *pair = stList_get(pairs, i);
            int64_t x = stIntTuple_get(pair, 0);
            int64_t y = stIntTuple_get(pair, 1);
            bool nonOverlapping = 1;
            for (int64_t j = 0; j < stList_length(pairs); j++) {
                stIntTuple *pair2 = stList_get(pairs, j);
                int64_t x2 = stIntTuple_get(pair2, 0);
                int64_t y2 = stIntTuple_get(pair2, 1);
                if ((x2 <= x && y2 >= y) || (x2 >= x && y2 <= y)) {
                    nonOverlapping = 0;
                    break;
                }
            }
            if (nonOverlapping) {
                stList_append(nonoverlappingPairs2, pair);
            }
        }
        stSortedSet *nonOverlappingPairsSet = stList_getSortedSet(nonoverlappingPairs,
                                                                  (int (*)(const void *, const void *)) stIntTuple_cmpFn);
        stSortedSet *nonOverlappingPairsSet2 = stList_getSortedSet(nonoverlappingPairs2,
                                                                   (int (*)(const void *, const void *)) stIntTuple_cmpFn);
        st_uglyf("The non-overlapping set sizes are %" PRIi64 " %" PRIi64 "\n",
                    stSortedSet_size(nonOverlappingPairsSet), stSortedSet_size(nonOverlappingPairsSet2));

        CuAssertTrue(testCase, stSortedSet_equals(nonOverlappingPairsSet, nonOverlappingPairsSet2));

        //Cleanup
        stSortedSet_destruct(nonOverlappingPairsSet);
        stSortedSet_destruct(nonOverlappingPairsSet2);
        stList_destruct(nonoverlappingPairs2);
        stList_destruct(pairs);
        stList_destruct(nonoverlappingPairs);

    }
}
*/
static void test_getBlastPairsWithRecursion(CuTest *testCase) {
    /*
     * Test the blast heuristic to get the different pairs.
     */
    for (int64_t test = 0; test < 10; test++) {
        //Make a pair of sequences
        char *seqX = getRandomSequence(st_randomInt(0, 10000));
        char *seqY = evolveSequence(seqX); //stString_copy(seqX);
        int64_t lX = strlen(seqX), lY = strlen(seqY);

        PairwiseAlignmentParameters *p = pairwiseAlignmentBandingParameters_construct();

        stList *blastPairs = getBlastPairsForPairwiseAlignmentParameters(seqX, seqY, p);

        checkBlastPairs(testCase, blastPairs, lX, lY, 1);
        stList_destruct(blastPairs);
        free(seqX);
        free(seqY);
    }
}

static void test_create_ambig_bases(CuTest *testCase) {
  stHash* something = create_ambig_bases();
  char* ambig_1 = "X";
  char* jt = stHash_search(something, ambig_1);
  CuAssertStrEquals(testCase, "ACGT", jt);
  stHash_destruct(something);
}


static void test_stList_construct(CuTest *testCase) {
  stList *methyls = stList_construct();
  char* ambig_1 = "X";
  stList_append(methyls, ambig_1);

  char* kmer = stList_remove(methyls, 0);

  CuAssertStrEquals(testCase, "X", kmer);
  stList_destruct(methyls);
}


static void test_path_permutePattern(CuTest *testCase){
  stList *methylPatterns = stList_construct3(0, &free);
  char *arr = (char *) malloc(sizeof(char) * 2);
  path_permutePattern(methylPatterns, 0, 2, arr, 2, "AB");
  char *pattern = stList_get(methylPatterns, 0);
  char *pattern2 = stList_get(methylPatterns, 1);
  CuAssertStrEquals(testCase, "AA", pattern);
  CuAssertStrEquals(testCase, "AB", pattern2);
}

CuSuite *signalPairwiseAlignerTestSuite(void) {
  CuSuite *suite = CuSuiteNew();

  SUITE_ADD_TEST(suite, test_create_ambig_bases);
  SUITE_ADD_TEST(suite, test_path_permutePattern);
  SUITE_ADD_TEST(suite, test_stList_construct);

  SUITE_ADD_TEST(suite, test_bands);
  SUITE_ADD_TEST(suite, test_diagonal);
  SUITE_ADD_TEST(suite, test_logAdd);
  SUITE_ADD_TEST(suite, test_Sequence);
  SUITE_ADD_TEST(suite, test_referenceSequence);
  SUITE_ADD_TEST(suite, test_eventSequence);
  SUITE_ADD_TEST(suite, test_loadNanoporeRead);
  SUITE_ADD_TEST(suite, test_1dNanoporeRead);
  SUITE_ADD_TEST(suite, test_getSplitPoints);
  SUITE_ADD_TEST(suite, test_hdCellConstruct);
  SUITE_ADD_TEST(suite, test_hdCellConstructWorstCase);
  SUITE_ADD_TEST(suite, test_dpDiagonal);
  SUITE_ADD_TEST(suite, test_dpMatrix);
  SUITE_ADD_TEST(suite, test_getBlastPairs);
  SUITE_ADD_TEST(suite, test_getBlastPairsWithRecursion);
  SUITE_ADD_TEST(suite, test_stateMachine3EmissionsPdfs);

    //SUITE_ADD_TEST(suite, test_filterToRemoveOverlap);  // wonky
    

    //SUITE_ADD_TEST(suite, test_poissonPosteriorProb);
    return suite;
}
