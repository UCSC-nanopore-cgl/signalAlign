#!/usr/bin/env python

from setuptools import setup, find_packages, Extension
import os
import numpy as np


c_compile_args = ['-pedantic', '-Wall', '-std=c99']
optimisation = ['-DNDEBUG', '-fstrict-aliasing']
c_compile_args.extend(optimisation)


event_detect = 'eventdetection'
pkg_path = os.path.join(os.path.dirname(__file__), event_detect)

include_dirs = [event_detect]

extensions = []

extensions.append(Extension(
    'nanonetfilters',
    sources=[os.path.join(pkg_path, 'filters.c')],
    include_dirs=include_dirs,
    extra_compile_args=c_compile_args
))

extensions.append(Extension("signalalign.cparsers", sources=[os.path.join(pkg_path, 'cparsers.c')], include_dirs=[np.get_include()]))


setup(name="signalAlign",
      version="0.1.6",
      description="A library for signal-level analysis of ONT data",
      author="Art Rand",
      author_email="arand@soe.ucsc.edu",
      url="https://github.com/UCSC-nanopore-cgl/signalAlign",
      package_dir={"": "src"},
      ext_modules=extensions,
      packages=find_packages("src"),
      install_requires=["numpy>=1.9.2",
                        "h5py>=2.2.1",
                        "pysam>=0.8.2.1",
                        "pandas>=0.18.1",
                        "sonLib>=1.1.0",
                        "PyYAML>=3.12",
                        "Cython>=0.26"]
      )
