#!/usr/bin/env python
"""Plot the number of breaks per number of bases for each read"""
########################################################################
# File: plot_accuracy_vs_alignment_deviation.py
#  executable: plot_accuracy_vs_alignment_deviation.py
#
# Author: Andrew Bailey
# History: Created 02/12/19
########################################################################

import os
import sys
from argparse import ArgumentParser
import numpy as np
import platform
import matplotlib as mpl

if os.environ.get('DISPLAY', '') == '':
    print('no display found. Using non-interactive Agg backend')
    mpl.use('Agg')
if platform.system() == "Darwin":
    mpl.use("macosx")
import matplotlib.pyplot as plt

plt.style.use('seaborn-deep')

from timeit import default_timer as timer
from py3helpers.utils import load_json, create_dot_dict, list_dir
from signalalign.alignedsignal import CreateLabels
from signalalign.utils import multithread


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to json config file")

    args = parser.parse_args()
    return args


def multiprocess_get_gaps_from_reads(in_dir, sa_number, gap_threshold, label, worker_count=1,
                                     alignment_threshold=0.5, debug=False):
    """Multiprocess for calculating the number of gaps in a signalalign alignment. Must have MEA alignment
    :param threshold: probability threshold
    :param sa_number: SA embedded number to grab sa variant data
    :param label: what is the correct nucleotide call for a given variant called read
    :param in_dir: input directory with subdirectories assumed to have fast5s in them
    :param worker_count: number of workers to use
    :return: True
    """
    if debug:
        output = []
        for f5_path in list_dir(in_dir, ext="fast5"):
            data = get_gap_lengths_and_read_length(f5_path, sa_number, label, gap_threshold,
                                                   alignment_threshold=alignment_threshold)
            if data is not None:
                output.append(data)

    filter_reads_args = {"sa_number": sa_number, "gap_threshold": gap_threshold, "label": label,
                         "alignment_threshold": alignment_threshold}
    total, failure, messages, output = multithread.run_service2(
        get_distance_from_guide_service, list_dir(in_dir, ext="fast5"),
        filter_reads_args, ["f5_path"], worker_count)
    return output


def get_distance_from_guide_service(work_queue, done_queue, service_name="get_gap_lengths_and_read_length_service"):
    """
    Service used by the multithread module in signal align to filter reads from a large number of directories
    :param work_queue: arguments to be done
    :param done_queue: errors and returns to be put
    :param service_name: name of the service
    """
    # prep
    total_handled = 0
    failure_count = 0
    mem_usages = list()

    # catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                reads = get_gap_lengths_and_read_length(**f)
                done_queue.put(reads)
            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, multithread.current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, multithread.current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, multithread.current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(multithread.TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(multithread.FAILURE_KEY, failure_count))
        if len(mem_usages) > 0:
            done_queue.put("{}:{}".format(multithread.MEM_USAGE_KEY, ",".join(map(str, mem_usages))))


def get_gap_lengths_and_read_length(f5_path, sa_number, label, gap_threshold=30, alignment_threshold=0.5):
    new_data = []
    mea_data = []
    percent_correct = None
    try:
        cl_handle = CreateLabels(f5_path, kmer_index=2)
        mea_data = cl_handle.add_mea_labels(number=sa_number)
        new_data = find_gaps_in_alignment_data(mea_data, min_gap=gap_threshold)

        data = cl_handle.add_variant_data(number=sa_number)
        bools = data[label] > alignment_threshold
        bools = bools.to_frame(name='true_false')
        percent_correct = np.mean(bools)

        cl_handle.close()
    except KeyError as e:
        print("{}: {}".format(f5_path, e))
        cl_handle.close()

    return new_data, len(mea_data), percent_correct


def find_gaps_in_alignment_data(ref_event_alignment, min_gap=30):
    """Calculate the number of alignment gaps over the minimum number of signals allowed"""
    gaps = []
    prev_end = ref_event_alignment[0]["raw_start"] + ref_event_alignment[0]["raw_length"]
    for x in ref_event_alignment:
        if x["raw_start"] - prev_end > min_gap:
            gaps.append(x["raw_start"] - prev_end)

        prev_end = x["raw_start"] + x["raw_length"]
    return gaps


def plot_n_breaks_vs_n_events_scatter_plot(all_events, all_breaks, names, all_percent_correct, save_fig_path=None):
    plt.figure(figsize=(10, 8))
    panel1 = plt.axes([0.1, 0.1, .8, .8])
    panel1.set_xlabel('Time from start')
    panel1.set_ylabel('Delta from guide alignment (individual current readings)')
    colors = [[0.0, 0.0, 1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0],
              [1.0, 0.0, 1.0], [1.0, 0.5, 0.0], [0.0, 1.0, 1.0], [0.0, 0.0, 0.0]]
    color_selection = 0

    for n_events, breaks, name, percent_correct in zip(all_events, all_breaks, names, all_percent_correct):
        rgba_colors = np.tile(np.array(colors[color_selection]), (len(n_events), 1))
        rgba_colors = np.insert(rgba_colors, 3, percent_correct, axis=1)

        panel1.scatter(n_events, [len(x) for x in breaks if x is not None], label=name, color=rgba_colors)
        color_selection += 1

    plt.legend(loc="lower right")
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def plot_n_breaks_ratio_vs_accuracy_scatter_plot(all_events, all_breaks, names, all_percent_correct, save_fig_path=None):
    plt.figure(figsize=(10, 8))
    panel1 = plt.axes([0.1, 0.1, .8, .8])
    panel1.set_xlabel('Ratio of number of gaps in a read / number of events in a read')
    panel1.set_ylabel('Accuracy of variant calls')

    for n_events, breaks, name, percent_correct in zip(all_events, all_breaks, names, all_percent_correct):

        panel1.scatter([len(x)/n_events[i] for i, x in enumerate(breaks) if x is not None if n_events[i]], [percent_correct[i] for i, x in enumerate(breaks) if x is not None if n_events[i]], label=name)

    plt.legend(loc="lower right")
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def main(config=None):
    """Plot number of breaks in a read over some difference and plot compared to accuracy of the variant called read"""
    start = timer()
    if config is None:
        args = parse_args()
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    args = create_dot_dict(config)
    threshold = args.gap_threshold
    names = []
    all_event_lengths = []
    all_skips = []
    all_percent_correct = []
    data = []
    for experiment in args.plot:
        names.append(experiment.name)
        for sample in experiment.samples:
            data = multiprocess_get_gaps_from_reads(sample.embedded_fast5_dir, experiment.sa_number, label=sample.label,
                                                    gap_threshold=threshold,
                                                    worker_count=7, alignment_threshold=args.threshold,
                                                    debug=False)
        all_skips.append([x[0] for x in data])
        all_event_lengths.append([x[1] for x in data])
        all_percent_correct.append([x[2] for x in data])
    plot_n_breaks_vs_n_events_scatter_plot(all_event_lengths, all_skips, names, all_percent_correct,
                                           save_fig_path=os.path.join(args.save_fig_dir,
                                                                      "n_breaks_vs_n_events.png"))
    plot_n_breaks_ratio_vs_accuracy_scatter_plot(all_event_lengths, all_skips, names, all_percent_correct,
                                                 save_fig_path=os.path.join(args.save_fig_dir,
                                                                            "n_breaks_ratio_vs_accuracy.png"))
    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == '__main__':
    main()
