#!/usr/bin/env python
"""Plot the accuracy of a variant call as a function of distance from the guide alignment"""
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
import pandas as pd
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
from signalalign.alignedsignal import CreateLabels, get_distance_from_guide_alignment
from signalalign.utils import multithread


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--config', '-c', action='store',
                        dest='config', required=True, type=str, default=None,
                        help="Path to json config file")

    args = parser.parse_args()
    return args


def multiprocess_get_distance_from_guide(in_dir, sa_number, threshold, label, worker_count=1, debug=False):
    """Multiprocess for filtering reads but dont move the files
    :param threshold: probability threshold
    :param sa_number: SA embedded number to grab sa variant data
    :param label: what is the correct nucleotide call for a given variant called read
    :param in_dir: input directory with subdirectories assumed to have fast5s in them
    :param worker_count: number of workers to use
    :param debug: boolean option which will only use one process in order to fail if an error arises
    :return: True
    """
    # grab aligned segment
    if debug:
        output = []
        for f5_path in list_dir(in_dir, ext="fast5"):
            data = get_distance_from_guide_alignment_wrapper(f5_path, sa_number, threshold, label)
            if data is not None:
                output.append(data)
    else:
        filter_reads_args = {"sa_number": sa_number, "threshold": threshold, "label": label}
        total, failure, messages, output = multithread.run_service2(
            get_distance_from_guide_service, list_dir(in_dir, ext="fast5"),
            filter_reads_args, ["f5_path"], worker_count)
    return output


def get_distance_from_guide_service(work_queue, done_queue, service_name="get_distance_from_guide_service"):
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
                reads = get_distance_from_guide_alignment_wrapper(**f)
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


def get_distance_from_guide_alignment_wrapper(f5_path, sa_number, threshold, label):
    f5_data = None
    try:
        cl_handle = CreateLabels(f5_path, kmer_index=2)
        cl_handle.add_basecall_alignment_prediction(number=0)

        guide_data = pd.DataFrame([x for x in cl_handle.aligned_signal.prediction.values()][0])
        data = cl_handle.add_variant_data(number=sa_number)
        new_data = get_distance_from_guide_alignment(data, guide_data, reference_index_key="position",
                                                     minus_strand=cl_handle.aligned_signal.minus_strand)
        bools = new_data[label] > threshold
        bools = bools.to_frame(name='true_false')
        f5_data = pd.concat([new_data, bools], axis=1)
        cl_handle.close()
    except KeyError as e:
        print("{}: {}".format(f5_path, e))
        cl_handle.close()

    return f5_data


def plot_alignment_deviation(deviations, labels, save_fig_path=None, bins=None):
    """Plot deltas of alignments to the guide alignment along with the accuracy"""
    plt.figure(figsize=(6, 6))
    panel1 = plt.axes([0.2, 0.2, .6, .6])
    panel1.set_xlabel('Delta from guide alignment (individual current readings)')
    panel1.set_ylabel('Density')
    panel1.grid(color='black', linestyle='-', linewidth=1)
    if bins is None:
        bins = 10

    panel1.hist(deviations, bins=bins, ls='dashed', normed=True,
                label=labels, alpha=0.7)

    plt.legend(loc="lower right")
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def plot_violin_classication_alignment_deviation(deviations, labels, save_fig_path=None):
    """Plot deltas of alignments to the guide alignment along with the accuracy"""
    plt.figure(figsize=(10, 8))
    panel1 = plt.axes([0.1, 0.1, .8, .8])
    panel1.set_xlabel('Classifcation Outcome')
    panel1.set_ylabel('Delta from guide alignment (individual current readings)')
    panel1.get_xaxis().set_tick_params(direction='out')
    panel1.set_xticklabels(labels)
    panel1.set_xlim(0.25, len(labels) + 0.75)
    panel1.xaxis.set_ticks_position('bottom')
    panel1.set_xticks(np.arange(1, len(labels) + 1))

    # panel1.grid(color='black', linestyle='-', linewidth=1)
    panel1.violinplot(deviations)

    plt.legend(loc="lower right")
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def plot_deviation_vs_time_from_start(deltas, times, labels, save_fig_path=None):
    """Plot deltas of alignments vs time from start"""
    plt.figure(figsize=(10, 8))
    panel1 = plt.axes([0.1, 0.1, .8, .8])
    panel1.set_xlabel('Time from start')
    panel1.set_ylabel('Delta from guide alignment (individual current readings)')

    for x in zip(times, deltas, labels):
        panel1.scatter(x[0], x[1], label=x[2], alpha=0.6)

    plt.legend(loc="lower right")
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def plot_classification_accuracy_vs_deviation(all_data, labels, save_fig_path=None):
    """Plot the accuracy vs delta from guide alignment"""
    deltas, percents = get_percent_accuracy_vs_deltas(all_data, n_bins=20)
    plt.figure(figsize=(10, 8))
    panel1 = plt.axes([0.1, 0.1, .8, .8])
    panel1.set_xlabel('Delta from guide alignment (individual current readings)')
    panel1.set_ylabel('Accuracy of calls')

    for x in zip(percents, labels):
        panel1.bar(deltas, width=deltas[1]-deltas[0], height=x[0], label=x[1], alpha=0.6)

    plt.legend(loc="lower right")
    if save_fig_path:
        plt.savefig(save_fig_path)
    else:
        plt.show()
        pass


def get_percent_accuracy_vs_deltas(all_data, n_bins=20):
    deltas = np.linspace(min(all_data[0]['guide_delta']), max(all_data[0]["guide_delta"]), n_bins)
    all_percents = []
    for data in all_data:
        max_index = len(data)
        index = 0
        total = 0
        n_right = 0
        percents = []
        data.sort_values("guide_delta", inplace=True)
        # cum_sum_data = np.cumsum(data["true_false"])
        row = data.iloc[index]
        for delta in deltas[1:]:
            while row["guide_delta"] < delta:
                total += 1
                n_right += row["true_false"]
                index += 1
                if index >= max_index:
                    break
                row = data.iloc[index]
            if total > 0:
                percents.append(n_right/total)
            else:
                percents.append(0)
            total = 0
            n_right = 0

        if index >= max_index:
            percents.append(0)
        else:
            row = data.iloc[index]
            while row["guide_delta"] >= delta:
                total += 1
                n_right += row["true_false"]
                index += 1
                if index >= max_index:
                    break

                row = data.iloc[index]
            if total > 0:
                percents.append(n_right/total)
            else:
                percents.append(0)

        all_percents.append(percents)
    return deltas, all_percents

def main(config=None):
    """Plot event to reference labelled ONT nanopore reads"""
    start = timer()
    if config is None:
        args = parse_args()
        # load model files
        assert os.path.exists(args.config), "Config file does not exist: {}".format(args.config)
        config = load_json(args.config)

    args = create_dot_dict(config)
    threshold = args.threshold
    all_data = []
    names = []
    for experiment in args.plot:
        names.append(experiment.name)
        experiment_data = []
        for sample in experiment.samples:
            tsvs = None
            f5s = None
            if sample.variant_tsvs is not None:
                tsvs = list_dir(sample.variant_tsvs, ext="vc.tsv")
            if sample.embedded_fast5_dir is not None:
                f5s = list_dir(sample.embedded_fast5_dir, ext="fast5")

            data = multiprocess_get_distance_from_guide(sample.embedded_fast5_dir, experiment.sa_number,
                                                                     threshold, sample.label,
                                                                     worker_count=7, debug=False)

            experiment_data.extend([x for x in data if x is not None])
        all_data.append(pd.concat(experiment_data))



    true_deltas = pd.concat([data[data["true_false"]]["guide_delta"] for data in all_data])
    true_starts = pd.concat([data[data["true_false"]]["raw_start"] for data in all_data])

    false_deltas = pd.concat([data[[not x for x in data["true_false"]]]["guide_delta"] for data in all_data])
    false_starts = pd.concat([data[[not x for x in data["true_false"]]]["raw_start"] for data in all_data])

    plot_deviation_vs_time_from_start([true_deltas, false_deltas], [true_starts, false_starts], ["True", "False"],
                                      os.path.join(args.save_fig_dir, "raw_start_vs_alignment_deviation_accuracy.png"))

    plot_deviation_vs_time_from_start([x["guide_delta"] for x in all_data], [x["raw_start"] for x in all_data], names,
                                      os.path.join(args.save_fig_dir, "raw_start_vs_alignment_deviation.png"))

    new_names = []
    new_data = []
    for name, data in zip(names, all_data):
        new_data.append(data[data["true_false"]]["guide_delta"])
        new_names.append(name+"_correct")
        new_data.append(data[[not x for x in data["true_false"]]]["guide_delta"])
        new_names.append(name+"_wrong")

    plot_alignment_deviation(new_data, new_names, bins=np.arange(-10000, 10000, 100),
                             save_fig_path=os.path.join(args.save_fig_dir, "alignment_deviation_hist.png"))

    plot_violin_classication_alignment_deviation(new_data, new_names,
                                                 save_fig_path=os.path.join(args.save_fig_dir,
                                                                            "alignment_deviation_violin.png"))

    plot_classification_accuracy_vs_deviation(all_data, names,
                                              save_fig_path=os.path.join(args.save_fig_dir,
                                                                         "classification_accuracy_vs_deviation.png"))

    stop = timer()
    print("Running Time = {} seconds".format(stop - start), file=sys.stderr)


if __name__ == '__main__':
    main()