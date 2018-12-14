#!/usr/bin/env python
"""Plot speeds of maximum expected accuracy methods"""
########################################################################
# File: plot_kmer_frequency.py
#  executable: plot_kmer_frequency.py
#
# Author: Andrew Bailey
# History: Created 08/28/18
########################################################################

from __future__ import print_function
import sys
import os
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.patches as mplpatches
import numpy as np
from collections import defaultdict
from timeit import default_timer as timer
from signalalign.signalAlignment import SignalAlignment
from py3helpers.utils import list_dir


def main():
    dna_tsv = "/Users/andrewbailey/CLionProjects/nanopore-RNN/submodules/signalAlign/tests/minion_test_reads/delete_me_after_debugging/6deaf971-6506-4e37-b486-cdf5e9d416ac.sm.forward.tsv"
    data = SignalAlignment.read_in_signal_align_tsv(dna_tsv, file_type='full')


if __name__ == "__main__":
    main()
    raise SystemExit
