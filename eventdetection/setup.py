#!/usr/bin/env python

from distutils.core import setup, Extension
import os
import numpy as np


event_detect = 'eventdetection'
pkg_path = os.path.join(os.path.dirname(__file__))
HOME = '/'.join(os.path.abspath(__file__).split("/")[:-2])

sonlib_include = os.path.join(HOME, "sonLib/C/inc")
sa_include = os.path.join(HOME, "inc")
h5_include = os.path.join(HOME, "include")
htsLib_include = os.path.join(HOME, "htslib/htslib")

main_lib = "/usr/local/lib"
son_Lib = os.path.join(HOME, "sonLib/lib")
h5_lib = os.path.join(HOME, "lib")
htsLib = os.path.join(HOME, "htslib")

signalAlign_a = os.path.join(HOME, "sonLib/lib/signalAlignLib.a")
h5_lib_a = os.path.join(h5_lib, "libhdf5.a")
son_Lib_a = os.path.join(HOME, "sonLib/lib/sonLib.a")
cu_test_a = os.path.join(HOME, "sonLib/lib/cuTest.a")

# libraries = ['z', 'sz', 'm', 'hts', 'dl']
# library_dirs = [main_lib, htsLib, son_Lib]
# extra_objects = [signalAlign_a, h5_lib_a, son_Lib_a]
# include_dirs = [event_detect, sa_headers, h5_include, sonlib_inc]
# lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm
# /usr/local
# -L${ttPrefix}/lib -Wl,-rpath,${ttPrefix}/lib -lkyototycoon -lkyotocabinet -lz -lbz2 -lpthread -lm -lstdc++
libraries = ['dl', 'z', 'sz', 'm', 'hts', 'pthread', 'gomp'] #'kyototycoon',
library_dirs = [h5_lib, son_Lib, htsLib]
library_dirs = []

extra_objects = [h5_lib_a, signalAlign_a, son_Lib_a]
include_dirs = [h5_include, sa_include, sonlib_include, htsLib_include]

c_compile_args = ['-pedantic', '-Wall', '-std=c99', '-mmacosx-version-min=10.11']
optimisation = ['-DNDEBUG', '-fstrict-aliasing', '-undefined', 'dynamic_lookup', '-fopenmp']

                # '-Wl,-rpath,/usr/local/lib', '-DHAVE_TOKYO_CABINET=1', '-DHAVE_KYOTO_TYCOON=1']
c_compile_args.extend(optimisation)
runtime_library_dirs = []

extensions = [Extension('kmeralign',
                        sources=[os.path.join(pkg_path, 'event_align_wrapper.c')],
                        include_dirs=include_dirs,
                        extra_compile_args=c_compile_args,
                        libraries=libraries,
                        library_dirs=library_dirs,
                        extra_objects=extra_objects)]

setup(name="kmeralign",
      version="1.0",
      description="Test Python wrapper on C function",
      author="Andrew Bailey",
      author_email="andbaile@ucsc.edu",
      ext_modules=extensions)
