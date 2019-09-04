cmake_minimum_required(VERSION 3.7)
project(signalAlign)


set(CMAKE_C_FLAGS_DEBUG  "${CMAKE_C_FLAGS_DEBUG} -O3 -g -Wall --pedantic -fopenmp -funroll-loops  -DNDEBUG")

#if(APPLE)
#    set(CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS "${CMAKE_SHARED_LIBRARY_CREATE_C_FLAGS} -undefined dynamic_lookup")
#endif(APPLE)

find_library(HDF5_LIBRARIES libhdf5.a
        HINTS ${CMAKE_SOURCE_DIR}/lib
        ${CMAKE_INSTALL_PREFIX}/lib
        /usr/lib)

find_path(HDF5_INCLUDE_DIRS hdf5.h
        HINTS ${HDF5_LIBRARIES}
        ${CMAKE_SOURCE_DIR}/include
        ${CMAKE_INSTALL_PREFIX}/include
        /usr/include)

message("CMakeInstall Prefix: ${CMAKE_INSTALL_PREFIX}")
message(HDF5_INCLUDE_DIRS: "${HDF5_INCLUDE_DIRS}")
message(HDF5_LIBRARIES: "${HDF5_LIBRARIES}")

find_library(HTSLIB_LIBRARY NAMES libhts.a
        HINTS ${CMAKE_SOURCE_DIR}/htslib
        ${CMAKE_INSTALL_PREFIX}/include)

find_path(HTSLIB_INCLUDE_DIR htslib/hts.h
        HINTS ${HTSLIB_LIBRARY}
        ${CMAKE_SOURCE_DIR}/htslib
        ${CMAKE_INSTALL_PREFIX}/include
        /usr/include)


set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY} )
set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR} )

message(HTSLIB_LIBRARIES: "${HTSLIB_LIBRARIES}")
message(HTSLIB_INCLUDE_DIRS: "${HTSLIB_INCLUDE_DIRS}")


include_directories(${HTSLIB_INCLUDE_DIRS})
include_directories(${HDF5_INCLUDE_DIRS})

link_directories(/usr/local/lib)
link_directories(/usr/lib)

include_directories($ENV{HOME}/anaconda3/envs/signalalign/include)
include_directories(${CMAKE_SOURCE_DIR}/inc)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/C/inc/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/cutest/)
include_directories(${CMAKE_SOURCE_DIR}/sonLib/externalTools/quicktree_1.1/include)
#include_directories(/usr/include)
#include_directories(/usr/local/include)
include_directories(/System/Library/Frameworks/Python.framework/Headers)


set(SA_LIBRARY
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a ${HTSLIB_LIBRARY} ${HDF5_LIBRARIES} -lm -lz -llzma -lbz2 -lcurl -lpthread -ldl)

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
set(HDPUTIL buildHdpUtil.c)


#add_executable(filters ${EVENTDETECT})

add_executable(signalMachine ${SOURCE_FILES} ${SIGNAL_MACHINE})
target_link_libraries(signalMachine ${SA_LIBRARY})
target_include_directories(signalMachine PUBLIC /usr/local/lib)

add_executable(buildHdpUtil ${SOURCE_FILES} ${HDPUTIL})
target_link_libraries(buildHdpUtil ${SA_LIBRARY})

add_executable(kmerEventAlign ${SOURCE_FILES} kmerEventAlign.c)
target_link_libraries(kmerEventAlign ${SA_LIBRARY})

add_executable(eventAlignerTests ${SOURCE_FILES} ${EVENTALIGNERTESTS})
target_link_libraries(eventAlignerTests ${SA_LIBRARY})

add_executable(signalAlignLibTests ${SOURCE_FILES} ${SIGNALMACHINE_TEST})
target_link_libraries(signalAlignLibTests  ${SA_LIBRARY})

#add_executable(extract ${SOURCE_FILES} ${EXTRACT})
#target_link_libraries(extract ${SA_LIBRARY})

#/usr/local/bin/gcc  -I/usr/local/include -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/inc -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/lib -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/cutest -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/quicktree_1.1/include -I/System/Library/Frameworks/Python.framework/Headers  -g -O3 -g -Wall --pedantic -fopenmp -funroll-loops -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk   -o CMakeFiles/signalMachine.dir/impl/hdp_math_utils.c.o   -c /Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/impl/hdp_math_utils.c
#gcc-5 -O3 -g -Wall --pedantic -fopenmp -funroll-loops -DNDEBUG   -I./include -I./htslib -fPIC -Iinc/ -I././sonLib/lib/  -c impl/*.c ./htslib/libhts.a ./lib/libhdf5.a -lz -llzma -lbz2 -lcurl -lpthread -lcrypto -ldl

#gcc-5 -O3 -g -Wall --pedantic -fopenmp -funroll-loops -DNDEBUG   -I./include -I./htslib -fPIC    -c impl/*.c ./htslib/libhts.a ./lib/libhdf5.a -lz -llzma -lbz2 -ldl
#/usr/local/bin/gcc  -I/usr/local/include  -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/cutest -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/quicktree_1.1/include -I/System/Library/Frameworks/Python.framework/Headers  -O3 -g -Wall --pedantic -fopenmp -funroll-loops -DNDEBUG -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk   -o CMakeFiles/buildHdpUtil.dir/impl/hdp_math_utils.c.o   -c /Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/impl/hdp_math_utils.c
#/r9_ecoli/small_batch_sa_output/tempFiles_trainModels/template_trained.hmm None /r9_ecoli/small_batch_sa_output/tempFiles_trainModels/template.singleLevelPrior2.nhdp None
#/usr/local/bin/gcc  -I/usr/local/include -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/inc -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/lib -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/cutest -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/quicktree_1.1/include -I/System/Library/Frameworks/Python.framework/Headers  -g -O3 -g -Wall --pedantic -fopenmp -funroll-loops  -DNDEBUG -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk   -o CMakeFiles/signalMachine.dir/impl/hdp_math_utils.c.o   -c /Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/impl/hdp_math_utils.c
#/usr/local/bin/gcc   -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/inc -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/lib -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/cutest -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/quicktree_1.1/include -I/System/Library/Frameworks/Python.framework/Headers  -g -O3 -g -Wall --pedantic -fopenmp -funroll-loops  -DNDEBUG -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk   -o CMakeFiles/signalMachine.dir/impl/hdp_math_utils.c.o   -c /Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/impl/hdp_math_utils.c
#/usr/local/bin/gcc  -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/htslib -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/include -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/inc -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/lib -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/cutest -I/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/sonLib/externalTools/quicktree_1.1/include -I/usr/local/include -I/System/Library/Frameworks/Python.framework/Headers  -g -O3 -g -Wall --pedantic -fopenmp -funroll-loops  -DNDEBUG -isysroot /Applications/Xcode.app/Contents/Developer/Platforms/MacOSX.platform/Developer/SDKs/MacOSX10.14.sdk   -o CMakeFiles/signalMachine.dir/impl/hdp_math_utils.c.o   -c /Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/impl/hdp_math_utils.c
