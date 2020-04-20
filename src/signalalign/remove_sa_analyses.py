#!/usr/bin/env python
"""Remove embedded signalalign analyses from files"""
########################################################################
# File: remove_sa_analyses.py
#  executable: remove_sa_analyses.py
#
# Author: Andrew Bailey
# History: 02/06/19 Created
########################################################################

import os
from py3helpers.utils import list_dir
from py3helpers.multiprocess import *
from argparse import ArgumentParser
from signalalign.fast5 import Fast5
import numpy as np


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--directory', '-d', required=True, action='store',
                        dest='dir', type=str, default=None,
                        help="Path to directory of fast5 files")
    parser.add_argument('--analysis', required=False, action='store_true',
                        dest='analysis', default=False,
                        help="Remove all analysis files")
    parser.add_argument('--basecall', required=False, action='store_true',
                        dest='basecall', default=False,
                        help="Remove all basecall files")
    parser.add_argument('--signalalign', required=False, action='store_true',
                        dest='signalalign', default=False,
                        help="Remove all signalalign files")
    parser.add_argument('--threads', required=False, action='store',
                        dest='threads', default=1, type=int,
                        help="number of threads to run")

    args = parser.parse_args()
    return args


def remove_sa_analyses(fast5):
    """Remove signalalign analyses from a fast5 file"""
    assert os.path.exists(fast5), "Fast5 path does not exist".format(fast5)
    fh = Fast5(fast5, read='r+')
    counter = 0
    for analyses in [x for x in list(fh["Analyses"].keys()) if "SignalAlign" in x]:
        fh.delete(os.path.join("Analyses", analyses))
        counter += 1
    fh.close()
    return counter


def remove_basecall_analyses(fast5):
    """Remove basecall analyses from a fast5 file"""
    assert os.path.exists(fast5), "Fast5 path does not exist".format(fast5)
    fh = Fast5(fast5, read='r+')
    counter = 0
    for analyses in [x for x in list(fh["Analyses"].keys()) if "Basecall" in x]:
        fh.delete(os.path.join("Analyses", analyses))
        counter += 1
    fh.close()
    return counter


def remove_analyses(fast5):
    """Remove analyses from a fast5 file"""
    assert os.path.exists(fast5), "Fast5 path does not exist".format(fast5)
    fh = Fast5(fast5, read='r+')
    counter = 0
    for analyses in [x for x in list(fh["Analyses"].keys())]:
        fh.delete(os.path.join("Analyses", analyses))
        counter += 1
    fh.delete("Analyses")
    fh.close()
    return counter


def main():
    args = parse_args()

    function_to_run = None
    if args.analysis:
        function_to_run = remove_analyses
    else:
        if args.signalalign or not args.basecall:
            function_to_run = remove_sa_analyses
        elif args.basecall:
            function_to_run = remove_basecall_analyses
    assert function_to_run is not None, "Must select --analysis, --signalalign or --basecall."

    service = BasicService(function_to_run, service_name="forward_multiprocess_aggregate_all_variantcalls")
    files = list_dir(args.dir, ext="fast5")
    total, failure, messages, output = run_service(service.run, files,
                                                   {}, ["fast5"], worker_count=args.threads)
    print("Deleted {} analysis datasets deleted from {} files".format(np.asarray(output).sum(), len(files)))


if __name__ == '__main__':
    main()
