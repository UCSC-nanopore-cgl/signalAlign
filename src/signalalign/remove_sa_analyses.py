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
from argparse import ArgumentParser
from signalalign.fast5 import Fast5


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--directory', '-d', required=True, action='store',
                        dest='dir', type=str, default=None,
                        help="Path to json config file")

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


def main():
    args = parse_args()
    total = 0
    files = 0
    errors = 0
    for f5_file in list_dir(args.dir, ext="fast5"):
        try:
            total += remove_sa_analyses(f5_file)
            files += 1
        except KeyError as e:
            errors += 1
            print("FAILED {}: {}".format(f5_file, e))

    print("Deleted {} SignalAlign analysis datasets deleted from {} files".format(total, files))


if __name__ == '__main__':
    main()
