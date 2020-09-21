############################################################################################################
cmake_minimum_required(VERSION 3.15)
project(signalAlign VERSION 0.2.1)

# set build shared libs option if not set
set(BUILD_SHARED_LIBS  "ON" CACHE STRING "Default build shared library")

if(NOT APPLE)
    set(CMAKE_INSTALL_RPATH $ORIGIN)
endif()
### Require out-of-source builds
file(TO_CMAKE_PATH "${PROJECT_BINARY_DIR}/CMakeLists.txt" LOC_PATH)
if(EXISTS "${LOC_PATH}")
    message(FATAL_ERROR "You cannot build in a source directory (or any directory with a CMakeLists.txt file). Please make a build subdirectory. Feel free to remove CMakeCache.txt and CMakeFiles.")
endif()

# add cmake dir
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" "${CMAKE_CURRENT_SOURCE_DIR}/cmake/modules")
#message("CMAKE_MODULE_PATH Prefix: ${CMAKE_MODULE_PATH}")
set(CMAKE_INSTALL_PREFIX .)
message(STATUS "CMakeInstall Prefix: ${CMAKE_INSTALL_PREFIX}")

############################################################################################################
# Find dependencies
############################################################################################################
include(${PROJECT_SOURCE_DIR}/cmake/AddDependencies.cmake)
############################################################################################################
# Create object library
############################################################################################################

add_library(sa_objlib OBJECT
        ${PROJECT_SOURCE_DIR}/impl/continuousHmm.c
        ${PROJECT_SOURCE_DIR}/impl/discreteHmm.c
        ${PROJECT_SOURCE_DIR}/impl/eventAligner.c
        ${PROJECT_SOURCE_DIR}/impl/event_detection.c
        ${PROJECT_SOURCE_DIR}/impl/fasta_handler.c
        ${PROJECT_SOURCE_DIR}/impl/hdp.c
        ${PROJECT_SOURCE_DIR}/impl/hdp_math_utils.c
        ${PROJECT_SOURCE_DIR}/impl/multipleAligner.c
        ${PROJECT_SOURCE_DIR}/impl/nanopore.c
        ${PROJECT_SOURCE_DIR}/impl/nanopore_hdp.c
        ${PROJECT_SOURCE_DIR}/impl/pairwiseAligner.c
        ${PROJECT_SOURCE_DIR}/impl/ranlib.c
        ${PROJECT_SOURCE_DIR}/impl/rnglib.c
        ${PROJECT_SOURCE_DIR}/impl/scrappie_common.c
        ${PROJECT_SOURCE_DIR}/impl/signalMachineUtils.c
        ${PROJECT_SOURCE_DIR}/impl/stateMachine.c
        ${PROJECT_SOURCE_DIR}/impl/util.c)

#${PROJECT_SOURCE_DIR}/inc/scrappie_stdlib.h
#${PROJECT_SOURCE_DIR}/inc/scrappie_structures.h
#${PROJECT_SOURCE_DIR}/inc/sse_mathfun.h
#${PROJECT_SOURCE_DIR}/inc/signalMachine.h
message(${HTSLIB_INCLUDE_DIRS})
target_link_libraries(sa_objlib PUBLIC ${signalalign_LINK_LIBRARIES})
target_include_directories(sa_objlib
        PUBLIC
        $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/signalalign>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/inc>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib>
        $<BUILD_INTERFACE:${signalalign_INCLUDE_DIRS}>
        )
# shared libraries need PIC
set_target_properties(sa_objlib PROPERTIES
        POSITION_INDEPENDENT_CODE 1)

target_compile_features(sa_objlib PUBLIC c_std_11)
target_compile_options(sa_objlib PRIVATE
#                -Wall
        #        -Werror
        #        -Wno-error=maybe-uninitialized
        #        -Wconversion
        #        -Wno-unused-resultON
        #        -Wpedantic
        $<$<CONFIG:RELEASE>:-g0>
        $<$<CONFIG:RELEASE>:-O3>
        $<$<CONFIG:DEBUG>:-O0>
        $<$<CONFIG:DEBUG>:-ggdb3>
        )
target_compile_definitions(sa_objlib PRIVATE
        $<$<CONFIG:RELEASE>:NDEBUG>
        $<$<CONFIG:RELEASE>:BOOST_DISABLE_ASSERTS>
        )
add_dependencies(sa_objlib sonLib_make)

############################################################################################################
# Create library
############################################################################################################
add_library(salib $<TARGET_OBJECTS:sa_objlib>)
target_link_libraries(salib PUBLIC
        sa_objlib)
set_target_properties(salib PROPERTIES
        OUTPUT_NAME "signalalign"
        PREFIX ""
        DEFINE_SYMBOL ""
        VERSION ${signalAlign_VERSION}
        SOVERSION ${signalAlign_VERSION_MAJOR}.${signalAlign_VERSION_MINOR})

############################################################################################################
# executable function
############################################################################################################

function(add_sa_executable executable source)
    add_executable(${executable} ${source})
    if (NOT APPLE AND NOT BUILD_SHARED_LIBS)
        #    message(STATUS "not APPLE")
        target_link_libraries(${executable} PUBLIC -static -static-libgcc -static-libstdc++)
    endif(NOT APPLE AND NOT BUILD_SHARED_LIBS)
    target_link_libraries(${executable} PUBLIC salib)
    add_custom_command(
            TARGET ${executable} POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_BINARY_DIR}/${executable}
            ${PROJECT_SOURCE_DIR}/bin/${executable})
endfunction()
############################################################################################################
# Create executables
############################################################################################################
add_sa_executable(signalMachine impl/signalMachine.c)
add_sa_executable(buildHdpUtil impl/buildHdpUtil.c)
add_sa_executable(compareDistributions impl/compareDistributions.c)
#add_sa_executable(estimateNanoporeParams impl/estimateNanoporeParams.c)
add_sa_executable(extract impl/extract.c)
add_sa_executable(kmerEventAlign impl/kmerEventAlign.c)
############################################################################################################
############################################################################################################
add_subdirectory(tests)

#SET(CMAKE_CXX_FLAGS "-g -Wall -std=c++11 -fopenmp -fsigned-char -DBOOST_NO_CXX11_SCOPED_ENUMS")
#set(CMAKE_VERBOSE_MAKEFILE ON)
#set(CMAKE_C_FLAGS_DEBUG  "${CMAKE_C_FLAGS_DEBUG} -O0 -g -Wall --pedantic -fopenmp -funroll-loops  -DNDEBUG")

#if(APPLE)
#    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
#endif(APPLE)

#include_directories($ENV{HOME}/anaconda3/envs/signalalign/include)
#include_directories(${CMAKE_SOURCE_DIR}/inc)
#include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
#include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
#include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
#include_directories(/usr/include)
#include_directories(/usr/local/include)
#include_directories(/System/Library/Frameworks/Python.framework/Headers)
#set(SIGNALMACHINE_TEST tests/allTests.c
#        tests/eventAlignerTests.c
#        tests/fastaHandlerTests.c
#        tests/hdpTests.c
#        tests/nanoporeHdpTests.c
#        tests/randomSequences.c
#        tests/signalPairwiseAlignerTest.c
#        tests/stateMachineTests.c
#        tests/variableOrderPairwiseAlignerTests.c)
#
#set(SIGNAL_MACHINE signalMachine.c)
##set(DEBUGGING debugging.c)
##set(EVENTDETECT eventdetection/filters.c)
#set(EVENTALIGNERTESTS tests/eventAlignerTests.c)
#set(EXTRACT extract.c)
#set(HDPUTIL buildHdpUtil.c)
#
#
##add_executable(filters ${EVENTDETECT})
#
#add_executable(signalMachine ${SOURCE_FILES} ${SIGNAL_MACHINE})
#target_link_libraries(signalMachine ${SA_LIBRARY})
#target_include_directories(signalMachine PUBLIC /usr/local/lib)
#
#add_executable(buildHdpUtil ${SOURCE_FILES} ${HDPUTIL})
#target_link_libraries(buildHdpUtil ${SA_LIBRARY})
#
#add_executable(kmerEventAlign ${SOURCE_FILES} kmerEventAlign.c)
#target_link_libraries(kmerEventAlign ${SA_LIBRARY})
#
#
#add_executable(signalAlignLibTests ${SOURCE_FILES} ${SIGNALMACHINE_TEST})
##target_link_libraries(signalAlignLibTests  ${SA_LIBRARY})
##add_dependencies(signalAlignLibTests sonLib_make)
#
#add_custom_command(
#        TARGET signalAlignLibTests POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_SOURCE_DIR}/bin/cPecanLastz_D
#        ${CMAKE_CURRENT_BINARY_DIR}/cPecanLastz_D)
#add_custom_command(
#        TARGET signalAlignLibTests POST_BUILD
#        COMMAND ${CMAKE_COMMAND} -E copy
#        ${CMAKE_SOURCE_DIR}/bin/cPecanLastz
#        ${CMAKE_CURRENT_BINARY_DIR}/cPecanLastz)
