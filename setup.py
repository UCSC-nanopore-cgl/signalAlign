#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
import os
import numpy as np

c_compile_args = ['-pedantic', '-Wall', '-std=c99']
optimisation = ['-DNDEBUG', '-fstrict-aliasing']
c_compile_args.extend(optimisation)

event_detect = 'eventdetection'
pkg_path = os.path.join(os.path.dirname(__file__), event_detect)

son_Lib = "sonLib/lib"

include_dirs = [event_detect, son_Lib]

extensions = []

extensions.append(Extension(
    'nanonetfilters',
    sources=[os.path.join(pkg_path, 'filters.c')],
    include_dirs=include_dirs,
    extra_compile_args=c_compile_args
))

extensions.append(
    Extension("signalalign.cparsers", sources=[os.path.join(pkg_path, 'cparsers.c')], include_dirs=[np.get_include()]))

HOME = os.path.dirname(__file__)

sonlib_include = os.path.join(HOME, "sonLib/C/inc")
sa_include = os.path.join(HOME, "inc")
h5_include = os.path.join(HOME, "include")
htsLib_include = os.path.join(HOME, "htslib/htslib")

signalAlign_a = os.path.join(HOME, "sonLib/lib/signalAlignLib.a")
h5_lib_a = os.path.join(HOME, "lib/libhdf5.a")
h5_hl_lib_a = os.path.join(HOME, "lib/libhdf5_hl.a")
son_Lib_a = os.path.join(HOME, "sonLib/lib/sonLib.a")
cu_test_a = os.path.join(HOME, "sonLib/lib/cuTest.a")

libraries = ['dl', 'z', 'sz', 'm', 'hts', 'pthread', 'gomp']

extra_objects = [h5_lib_a, signalAlign_a, son_Lib_a, h5_hl_lib_a]
include_dirs = [h5_include, sa_include, sonlib_include, htsLib_include]

c_compile_args = ['-pedantic', '-Wall', '-std=c99', '-mmacosx-version-min=10.11', '-DNDEBUG', '-fstrict-aliasing',
                  '-undefined', 'dynamic_lookup', '-fopenmp']

extensions.append(Extension('kmeralign',
                            sources=[os.path.join(pkg_path, 'event_align_wrapper.c')],
                            include_dirs=include_dirs,
                            extra_compile_args=c_compile_args,
                            libraries=libraries,
                            extra_objects=extra_objects))

setup(name="signalAlign",
      version="0.1.7",
      description="A library for signal-level analysis of ONT data",
      author="Art Rand / Andrew Bailey / Trevor Pesout",
      author_email="andbaile@ucsc.edu",
      url="https://github.com/UCSC-nanopore-cgl/signalAlign",
      package_dir={"": "src"},
      ext_modules=extensions,
      packages=find_packages("src"),
      install_requires=["numpy>=1.9.2",
                        "h5py>=2.2.1",
                        "pysam>=0.8.2.1",
                        "pandas>=0.23.1",
                        "sonLib>=1.1.0",
                        "PyYAML>=3.12",
                        "Cython>=0.26"]
      )
