cmake_minimum_required(VERSION 3.7)
project(signalAlign)

set(CMAKE_C_STANDARD 99)

#if(APPLE)
#    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
#endif(APPLE)

add_custom_target(build_nanoporeParams
        COMMAND make
        COMMAND make nanoporeParams
        COMMAND pip install -r requirements.txt
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})# -C /Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign)

add_custom_target(build_signalAlign
        COMMAND make
#        COMMAND make nanoporeParams
#        COMMAND pip install -r requirements.txt
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})# -C /Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign)


add_custom_target(nanotensor_reqs
        COMMAND pip install -r requirements.txt
        COMMAND sudo python setup.py install
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR})# -C /Users/andrewbailey/CLionProjects/nanopore-RNN/signalAlign)

include(${CMAKE_ROOT}/Modules/ExternalProject.cmake)

#ExternalProject_Add(sonLib
#        SOURCE_DIR ${CMAKE_SOURCE_DIR}/signalAlign/sonLib/C/
#        PREFIX ${CMAKE_CURRENT_BINARY_DIR}/sonLib
#        CONFIGURE_COMMAND ""
#        BUILD_COMMAND make
#        PREFIX=${CMAKE_CURRENT_BINARY_DIR}/sonLib
#        )

#find_package(Boost)
#include_directories(signalAlign/sonLib/C/inc/)
#set(sonLib signalAlign/sonLib/C/inc/)
link_directories(/usr/local/lib)
link_directories(${CMAKE_SOURCE_DIR}/inc)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/C/impl/)
link_directories(/usr/local/lib)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
link_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
link_directories(/usr/include/)

#add_subdirectory(signalAlign
#        BUILD_COMMAND make)

include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/impl/)
include_directories(/usr/local/lib)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
include_directories(/usr/include/)
include_directories(/System/Library/Frameworks/Python.framework/Headers)

set(estimateNanoporeParams estimateNanoporeParams.c)


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
        impl/banded_alignment.c
        impl/
        /Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/allTests.c)

set(SIGNAL_MACHINE signalMachine.c)
set(DEBUGGING debugging.c)

set(EVENTDETECT eventdetection/filters.c)

add_executable(filters ${EVENTDETECT})


add_custom_command(
        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/
)


#add_custom_command(
#        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/signalAlign/htslib/libhts.a
#        COMMAND make
#        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/signalAlign/htslib
#)

add_custom_target(
        signalAlign
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/signalAlignLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
        ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a)

#add_custom_target(
#        htslib
#        ${CMAKE_CURRENT_SOURCE_DIR}/signalAlign/htslib/libhts.a
#)

#add_executable(estimateNanoporeParams ${estimateNanoporeParams})
add_executable(signalMachine ${SOURCE_FILES} ${SIGNAL_MACHINE})
target_link_libraries(signalMachine -lpthread -lrt -lm -L${CMAKE_CURRENT_SOURCE_DIR}/htslib hts
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/signalAlignLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
        ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        )
add_dependencies(signalMachine signalAlign)



add_executable(debugging ${SOURCE_FILES} ${DEBUGGING})
target_link_libraries(debugging)
target_link_libraries(debugging -L${CMAKE_CURRENT_SOURCE_DIR}/htslib hts ${CMAKE_CURRENT_SOURCE_DIR}/htslib/libhts.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a

        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a)
add_dependencies(debugging signalAlign)



#        signalAlign/sonLib/C/impl/pairwiseAlignment.c
#        signalAlign/sonLib/C/impl/avl.c
#        signalAlign/sonLib/C/impl/bioioC.c
#        signalAlign/sonLib/C/impl/commonC.c
#        signalAlign/sonLib/C/impl/fastCMaths.c
#        signalAlign/sonLib/C/impl/hashTableC.c
#        signalAlign/sonLib/C/impl/hashTableC_itr.c
#        signalAlign/sonLib/C/impl/lz4.c
#        signalAlign/sonLib/C/impl/lz4.h
#        signalAlign/sonLib/C/impl/lz4hc.c
#        signalAlign/sonLib/C/impl/lz4hc.h
#        signalAlign/sonLib/C/impl/pairwiseAlignment.c
#        signalAlign/sonLib/C/impl/sonLibAlign.c
#        signalAlign/sonLib/C/impl/sonLibCache.c
#        signalAlign/sonLib/C/impl/sonLibCommon.c
#        signalAlign/sonLib/C/impl/sonLibCompression.c
#        signalAlign/sonLib/C/impl/sonLibConnectivity.c
#        signalAlign/sonLib/C/impl/sonLibEulerTour.c
#        signalAlign/sonLib/C/impl/sonLibExcept.c
#        signalAlign/sonLib/C/impl/sonLibFile.c
#        signalAlign/sonLib/C/impl/sonLibGlobalsInternal.h
#        signalAlign/sonLib/C/impl/sonLibHash.c
#        signalAlign/sonLib/C/impl/sonLibKVDatabase.c
#        signalAlign/sonLib/C/impl/sonLibKVDatabase_BigRecordFile.c
#        signalAlign/sonLib/C/impl/sonLibKVDatabase_MySql.c
#        signalAlign/sonLib/C/impl/sonLibKVDatabase_TokyoCabinet.c
#        signalAlign/sonLib/C/impl/sonLibKVDatabaseConf.c
#        signalAlign/sonLib/C/impl/sonLibKVDatabasePrivate.h
#        signalAlign/sonLib/C/impl/sonLibList.c
#        signalAlign/sonLib/C/impl/sonLibListPrivate.h
#        signalAlign/sonLib/C/impl/sonLibNaiveConnectivity.c
#        signalAlign/sonLib/C/impl/sonLibRandom.c
#        signalAlign/sonLib/C/impl/sonLibSet.c
#        signalAlign/sonLib/C/impl/sonLibSortedSet.c
#        signalAlign/sonLib/C/impl/sonLibString.c
#        signalAlign/sonLib/C/impl/sonLibTreap.c
#        signalAlign/sonLib/C/impl/sonLibTree.c
#        signalAlign/sonLib/C/impl/sonLibTuples.c
#        signalAlign/sonLib/C/impl/stGraph.c
#        signalAlign/sonLib/C/impl/stMatrix.c
#        signalAlign/sonLib/C/impl/stPhylogeny.c
#        signalAlign/sonLib/C/impl/stPosetAlignment.c
#        signalAlign/sonLib/C/impl/stSafeC.c
#        signalAlign/sonLib/C/impl/stThreadPool.c)

#add_custom_command(
#        OUTPUT ${CMAKE_CURRENT_SOURCE_DIR}/signalAlign/sonLib/lib/signalAlignLib.a
#        COMMAND make
#        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/signalAlign/sonLib/
#)
#
#add_custom_target(
#    sonlib
#    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/signalAlign/sonLib/lib/signalAlignLib.a)

