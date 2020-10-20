message(CHECK_START "Finding signalAlign Dependencies")
list(APPEND CMAKE_MESSAGE_INDENT "  ")

if(BUILD_SHARED_LIBS)
    message(STATUS "BUILDING SHARED LIBS")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .so .dylib ${CMAKE_FIND_LIBRARY_SUFFIXES})
else(NOT BUILD_SHARED_LIBS)
    message(STATUS "BUILDING STATIC LIBS")
    SET(CMAKE_FIND_LIBRARY_SUFFIXES .a ${CMAKE_FIND_LIBRARY_SUFFIXES})
endif(BUILD_SHARED_LIBS)

############################################################################################################
# find packages
############################################################################################################
find_package(HTSlib REQUIRED htslib)
find_package(Eigen3 REQUIRED NO_MODULE)
#set(HDF5_FIND_DEBUG TRUE)
set(HDF5_PREFER_PARALLEL TRUE)
find_package(HDF5 1.10.0 COMPONENTS C CXX REQUIRED)

############################################################################################################
# sonlib
############################################################################################################
set(sonLib_LIB
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/sonLib.a
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/lib/cuTest.a
        )
set(sonLib_inc
        ${CMAKE_CURRENT_SOURCE_DIR}/sonLib/externalTools/cutest/)

add_custom_command(
        OUTPUT ${sonLib_LIB}
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/sonLib
)
#
add_custom_target(sonLib_make DEPENDS ${sonLib_LIB})


############################################################################################################
# make dir
############################################################################################################
file(MAKE_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/bin/)
############################################################################################################
# lastz
############################################################################################################
set(lastz
        ${CMAKE_CURRENT_SOURCE_DIR}/bin/cPecanLastz
        ${CMAKE_CURRENT_SOURCE_DIR}/bin/cPecanLastz_D)
add_custom_command(
        OUTPUT ${lastz}
        COMMAND make
        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/externalTools
)
#
add_custom_target(lastz_make DEPENDS ${lastz})

############################################################################################################
# link libraries
############################################################################################################

set(signalalign_LINK_LIBRARIES
        ${sonLib_LIB}
        HTSlib::HTSlib
        ${CMAKE_DL_LIBS}
        m
        pthread)


if (NOT TARGET hdf5::hdf5)
    list(APPEND signalalign_LINK_LIBRARIES ${HDF5_LIBRARIES})
    set(signalalign_INCLUDE_DIRS
            ${HDF5_INCLUDE_DIRS})
    if(NOT APPLE)
        list(APPEND signalalign_LINK_LIBRARIES aec)
    endif(NOT APPLE)
else()
    list(APPEND signalalign_LINK_LIBRARIES hdf5::hdf5)
endif()
