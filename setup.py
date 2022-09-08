#!/usr/bin/env python
from setuptools import setup, find_packages  # , Extension
import re


def get_version():
    try:
        content = open("CMakeLists.txt", "r").read()
        version = re.search(r'project\(signalAlign VERSION (.*)\)', content).group(1)
        return version.strip()
    except RuntimeError:
        return None


setup(name="signalAlign",
      version=get_version(),
      description="A library for signal-level analysis of ONT data",
      author="Art Rand / Andrew Bailey / Trevor Pesout",
      author_email="andbaile@ucsc.edu",
      url="https://github.com/UCSC-nanopore-cgl/signalAlign",
      package_dir={"": "src"},
      # library_dirs=[os.path.join(HOME, "lib")],
      # ext_modules=extensions,
      packages=find_packages("src"),
      scripts=["src/signalalign/train/trainModels.py",
               "src/signalalign/scripts/runSignalAlign.py",
               "src/signalalign/remove_sa_analyses.py",
               "src/signalalign/filter_reads.py",
               "src/signalalign/visualization/compare_trained_models.py",
               "src/signalalign/visualization/plot_accuracy_vs_alignment_deviation.py",
               "src/signalalign/visualization/plot_breaks_in_alignments.py",
               "src/signalalign/visualization/plot_em_model_distributions.py",
               "src/signalalign/visualization/plot_kmer_distributions.py",
               "src/signalalign/visualization/plot_labelled_read.py",
               "src/signalalign/visualization/plot_multiple_variant_accuracy.py",
               "src/signalalign/visualization/plot_variant_accuracy.py",
               "src/signalalign/visualization/sequencing_summary.py",
               "src/signalalign/visualization/verify_load_from_raw.py"],
      install_requires=["numpy>=1.9.2",
                        "h5py>=2.2.1,<3.0.0",
                        "pysam>=0.8.2.1,<0.16.0",
                        "pandas>=0.23.1",
                        "sonLib>=1.1.0",
                        "PyYAML>=3.12",
                        "Cython>=0.26",
                        "scikit-learn>=0.19.0,<=0.20.3",
                        "matplotlib>=2.0.2",
                        "pathos==0.2.9",
                        "scipy>=1.5.0"
                        "py3helpers[seq_tools]>=0.5.0",
                        'embed>=0.0.5']
      )
