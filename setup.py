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

### nanonet filters ###

extensions.append(Extension(
    'nanonetfilters',
    sources=[os.path.join(pkg_path, 'filters.c')],
    include_dirs=include_dirs,
    extra_compile_args=c_compile_args
))


### cparsers ###

extensions.append(
    Extension("signalalign.cparsers", sources=[os.path.join(pkg_path, 'cparsers.c')], include_dirs=[np.get_include()]))


### kmer align ###

HOME = os.path.abspath(os.path.dirname(__file__))

sonlib_include = os.path.join(HOME, "sonLib/C/inc")
sa_include = os.path.join(HOME, "inc")

signalAlign_a = os.path.join(HOME, "sonLib/lib/signalAlignLib.a")
son_Lib_a = os.path.join(HOME, "sonLib/lib/sonLib.a")
cu_test_a = os.path.join(HOME, "sonLib/lib/cuTest.a")

libraries = ['dl', 'z', 'm', 'pthread', 'gomp', 'hdf5']
extra_objects = [signalAlign_a, son_Lib_a]
include_dirs = [sa_include, sonlib_include]
c_compile_args = ['-pedantic', '-Wall', '-std=c99', '-DNDEBUG', '-fstrict-aliasing', '-fopenmp',
                  '-L{}'.format(os.path.join(HOME, 'lib'))]

extensions.append(Extension('kmeralign',
                            sources=[os.path.join(pkg_path, 'event_align_wrapper.c')],
                            include_dirs=include_dirs,
                            extra_compile_args=c_compile_args,
                            libraries=libraries,
                            extra_objects=extra_objects))

setup(name="signalAlign",
      version="0.2.0",
      description="A library for signal-level analysis of ONT data",
      author="Art Rand / Andrew Bailey / Trevor Pesout",
      author_email="andbaile@ucsc.edu",
      url="https://github.com/UCSC-nanopore-cgl/signalAlign",
      package_dir={"": "src"},
      # library_dirs=[os.path.join(HOME, "lib")],
      ext_modules=extensions,
      packages=find_packages("src"),
      install_requires=["numpy>=1.9.2",
                        "h5py>=2.2.1",
                        "pysam>=0.8.2.1",
                        "pandas>=0.23.1",
                        "sonLib>=1.1.0",
                        "PyYAML>=3.12",
                        "Cython>=0.26",
                        "scikit-learn==0.19.0",
                        "matplotlib==2.0.2",
                        "pathos==0.2.1",
                        "py3helpers==0.1.2"]
      )
