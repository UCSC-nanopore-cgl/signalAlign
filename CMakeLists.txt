cmake_minimum_required(VERSION 3.7)
project(signalAlign)

set(CMAKE_C_STANDARD 99)

#if(APPLE)
#    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
#endif(APPLE)


add_custom_target(build_signalAlign
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(build_nanoporeParams
        COMMAND make
        COMMAND make nanoporeParams
        COMMAND pip install -r requirements.txt
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

add_custom_target(nanotensor_reqs
        COMMAND pip install -r requirements.txt
        COMMAND sudo python setup.py install
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})

link_directories(/usr/local/lib)
link_directories(${CMAKE_SOURCE_DIR}/inc)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/C/impl/)
link_directories(/usr/local/lib)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
link_directories(/usr/include/)

include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/impl/)
include_directories(/usr/local/lib)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
include_directories(/usr/include/)
include_directories(/System/Library/Frameworks/Python.framework/Headers)

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
        impl/
        )

set(SIGNAL_MACHINE signalMachine.c)
set(DEBUGGING debugging.c)

set(EVENTDETECT eventdetection/filters.c)


add_executable(filters ${EVENTDETECT})

add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
)

add_custom_target(
        signalAlign
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/signalAlignLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
        ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

add_executable(signalMachine ${SOURCE_FILES} ${SIGNAL_MACHINE})
target_link_libraries(signalMachine -lpthread -lrt -lm -L${CMAKE_CURRENT_SOURCE_DIR}/htslib hts
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/signalAlignLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
        ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        )
add_dependencies(signalMachine signalAlign)


#include_directories(${CMAKE_SOURCE_DIR}/include)
#add_executable(kmerEventAlignTests ${SOURCE_FILES} impl/eventAligner.c tests/eventAlignerTests.c)
#target_link_libraries(kmerEventAlignTests -lz -lrt -lm -L${CMAKE_CURRENT_SOURCE_DIR}/htslib hts
#        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
#        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/signalAlignLib.a
#        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
#        ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
#        ${CMAKE_CURRENT_SOURCE_DIR}/lib/libhdf5.a
#        -lz -ldl -lpthread
#        )
#add_dependencies(kmerEventAlignTests signalAlign)
