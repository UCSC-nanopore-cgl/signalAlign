## FindHTSlib.cmake
#
# Finds the HTSlib library
#
# This will define the following variables
#
#    HTSlib_FOUND
#    HTSlib_INCLUDE_DIRS
#
# and the following imported targets
#
#     HTSlib::HTSlib
#
# Author Nicholas Knoblauch (borrowing heavily from): find libxml2

find_package(PkgConfig)
find_package(PkgConfig QUIET)
PKG_CHECK_MODULES(PC_HTSLIB htslib)
set(HTSLIB_DEFINITIONS ${PC_HTSLIB_CFLAGS_OTHER})

find_path(HTSLIB_INCLUDE_DIR NAMES hts.h
   HINTS
   ${PC_HTSLIB_INCLUDEDIR}
   ${PC_HTSLIB_INCLUDE_DIRS}
   PATH_SUFFIXES htslib
   )

# CMake 3.9 and below used 'HTSLIB_LIBRARIES' as the name of
# the cache entry storing the find_library result.  Use the
# value if it was set by the project or user.
if(DEFINED HTSLIB_LIBRARIES AND NOT DEFINED HTSLIB_LIBRARY)
  set(HTSLIB_LIBRARY ${HTSLIB_LIBRARIES})
endif()

find_library(HTSLIB_LIBRARY NAMES hts
   HINTS
   ${PC_HTSLIB_LIBDIR}
   ${PC_HTSLIB_LIBRARY_DIRS}
   )

#message( STATUS "HTSLIB_FOUND: ${PC_HTSLIB_FOUND}" ) # no output for this

if(PC_HTSLIB_VERSION)
#  message( STATUS "HTSLIB_VERSION: ${PC_HTSLIB_VERSION}" ) # no output for this
    set(HTSLIB_VERSION_STRING ${PC_HTSLIB_VERSION})
elseif(HTSLIB_INCLUDE_DIR AND EXISTS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h")
    file(STRINGS "${HTSLIB_INCLUDE_DIR}/htslib/hts.h" htslib_version_str
         REGEX "^#define[\t ]+HTS_VERSION[\t ]+\".*\"")

    string(REGEX REPLACE "^#define[\t ]+HTS_VERSION[\t ]+\"([^\"]*)\".*" "\\1"
           HTSLIB_VERSION_STRING "${htslib_version_str}")
    unset(htslib_version_str)
endif()

#set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR} ${PC_HTSLIB_INCLUDE_DIRS})
set(HTSLIB_INCLUDE_DIRS ${HTSLIB_INCLUDE_DIR})

set(HTSLIB_LIBRARIES ${HTSLIB_LIBRARY})

FIND_PACKAGE_HANDLE_STANDARD_ARGS(HTSlib
                                  REQUIRED_VARS HTSLIB_LIBRARY HTSLIB_INCLUDE_DIR
                                  VERSION_VAR HTSLIB_VERSION_STRING)

mark_as_advanced(HTSLIB_INCLUDE_DIR HTSLIB_LIBRARY)

find_package(BZip2 REQUIRED)
find_package(LibLZMA REQUIRED)
find_package(ZLIB REQUIRED)
find_package(CURL REQUIRED)

if(HTSlib_FOUND AND NOT TARGET HTSlib::HTSlib)
  add_library(HTSlib::HTSlib UNKNOWN IMPORTED)
  target_link_libraries(HTSlib::HTSlib INTERFACE
          lzma
          bz2
          curl
          z)
  set_target_properties(HTSlib::HTSlib PROPERTIES
          INTERFACE_INCLUDE_DIRECTORIES "${HTSLIB_INCLUDE_DIRS}"
          )
  set_property(TARGET HTSlib::HTSlib APPEND PROPERTY IMPORTED_LOCATION "${HTSLIB_LIBRARY}")
endif()

