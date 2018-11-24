cmake_minimum_required(VERSION 3.7)
project(signalAlign)

set(CMAKE_C_STANDARD 99)

add_custom_target(build_signalAlign
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

link_directories(/usr/local/lib)

include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
include_directories(/usr/local/include)
include_directories(/usr/include)

include_directories(/System/Library/Frameworks/Python.framework/Headers)

set(SA_LIBRARY -lz -lm -lhts -lhdf5
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/signalAlignLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
        -lz -ldl -lpthread )

set(SOURCE_FILES
        impl/stateMachine.c
        impl/signalMachineUtils.c
        impl/nanopore.c
        impl/continuousHmm.c
        impl/discreteHmm.c
        impl/hdp.c
        impl/hdp_math_utils.c
        impl/multipleAligner.c
        impl/nanopore.c
        impl/nanopore_hdp.c
        impl/pairwiseAligner.c
        impl/ranlib.c
        impl/rnglib.c
        impl/signalMachineUtils.c
        impl/fasta_handler.c
        impl/eventAligner.c
        impl/event_detection.c
        impl/scrappie_common.c
        impl/util.c
        inc/scrappie_stdlib.h
        inc/scrappie_structures.h
        inc/sse_mathfun.h
        inc/signalMachine.h)

set(SIGNALMACHINE_TEST tests/allTests.c
        tests/eventAlignerTests.c
        tests/fastaHandlerTests.c
        tests/hdpTests.c
        tests/nanoporeHdpTests.c
        tests/randomSequences.c
        tests/signalPairwiseAlignerTest.c
        tests/stateMachineTests.c
        tests/variableOrderPairwiseAlignerTests.c)

set(SIGNAL_MACHINE signalMachine.c)
set(DEBUGGING debugging.c)
#set(EVENTDETECT eventdetection/filters.c)
set(EVENTALIGNERTESTS tests/eventAlignerTests.c)
set(EXTRACT extract.c)


#add_executable(filters ${EVENTDETECT})

add_executable(signalMachine ${SOURCE_FILES} ${SIGNAL_MACHINE})
target_link_libraries(signalMachine ${SA_LIBRARY})

add_executable(kmerEventAlign ${SOURCE_FILES} kmerEventAlign.c)
target_link_libraries(kmerEventAlign ${SA_LIBRARY})

add_executable(eventAlignerTests ${SOURCE_FILES} ${EVENTALIGNERTESTS})
target_link_libraries(eventAlignerTests ${SA_LIBRARY})

add_executable(signalAlignLibTests ${SOURCE_FILES} ${SIGNALMACHINE_TEST})
target_link_libraries(signalAlignLibTests  ${SA_LIBRARY})

add_executable(extract ${SOURCE_FILES} ${EXTRACT})
target_link_libraries(extract ${SA_LIBRARY})
