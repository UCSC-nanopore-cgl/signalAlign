
from argparse import ArgumentParser
import subprocess
import os
import hashlib
import sys
import shutil

from signalalign.utils.multithread import *
from signalalign.tools.fast5_lookup import RUN_NAME, FAST5_LOCATION, READ_ID, FAST5_ROOT


KEY_ID = "id"
KEY_S3_ROOT = "s3_root"
KEY_READS = "reads"
KEY_SUMMARY = "summary"
KEY_FASTQ = "fastq"

KEY_READ_COUNT = "count"


def parse_args():
    parser = ArgumentParser(description="Build indices from fast5 output stored in S3")

    parser.add_argument('--s3_root', '-i', action='store', dest='s3_root', required=True, type=str,
                        help="s3 location with files")
    parser.add_argument('--workdir', '-w', action='store', dest='workdir', type=str, default="/tmp",
                        help="tmp directory for downloading temporary files")
    parser.add_argument('--index_destination', '-o', action='store', dest='index_destination', required=True, type=str,
                        help="destination for index files")
    parser.add_argument('--threads', '-t', action='store', dest='threads', required=False, type=int, default=1,
                        help="number of threads to use while indexing")

    args = parser.parse_args()
    return args


def log(msg):
    print(msg)


def fast5_promethion_s3_organization_service(work_queue, done_queue, service_name="fast5_organization"):
    # prep
    total_handled = 0
    failure_count = 0
    total_reads = 0

    #catch overall exceptions
    try:
        for f in iter(work_queue.get, 'STOP'):
            # catch exceptions on each element
            try:
                success, count = organize_promethion_s3_fast5s(**f)
                if success:
                    total_reads += count
            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            # increment total handling
            total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))
        done_queue.put("{}:{}".format(KEY_READ_COUNT, total_reads))


def organize_promethion_s3_fast5s(workdir, fast5_data_locations, destination_dir, keep_workdir=False):
    # get run info
    reads = fast5_data_locations[KEY_READS].rstrip("/")
    summary = fast5_data_locations[KEY_SUMMARY]

    id = fast5_data_locations[KEY_ID]
    s3_root = fast5_data_locations[KEY_S3_ROOT]
    summary_s3 = os.path.join(s3_root, summary)

    # evaluate completedness
    run_id = "f5run_{}_{}".format(hashlib.md5(s3_root.encode()).hexdigest(), id)
    completed_index_location = os.path.join(destination_dir, "{}.index.tsv".format(run_id))
    if os.path.isfile(completed_index_location):
        log("Indexing appears to be completed for:")
        log("\ts3_root:     {}".format(s3_root))
        log("\tid:          {}".format(id))
        log("\tdestination: {}".format(completed_index_location))
        return True, 0

    # create workdir
    run_workdir = os.path.join(workdir, run_id)
    os.mkdir(run_workdir)

    # get s3files
    s3cmd_get_cmd = ['s3cmd', 'get', summary_s3]
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(s3cmd_get_cmd, cwd=run_workdir, stdout=devnull)

    # get locations and sanity check
    summary_loc = os.path.join(run_workdir, os.path.basename(summary_s3))
    if not os.path.isfile(summary_loc):
        err = "Error downloading files from {}.  Missing {} ({}).".format(s3_root, summary_loc, summary_s3)
        log(err)
        return False, err

    # todo copy relevant information to textfile
    stored_reads = 0
    fast5_idx, read_id_idx, run_id_idx = None, None, None
    with open(summary_loc, 'r') as input, open(completed_index_location, 'w') as output:
        output.write("##{}:{}\n".format(FAST5_ROOT, s3_root))
        output.write("#{}\t{}\t{}\n".format(FAST5_LOCATION, READ_ID, RUN_NAME))
        for line in input:
            line = line.split("\t")
            if len(line) < 3: continue

            #header
            if fast5_idx is None:
                for column, idx in zip(line, range(len(line))):
                    if column == 'filename': fast5_idx = idx
                    if column == 'read_id': read_id_idx = idx
                    if column == 'run_id': run_id_idx = idx
                if None in [fast5_idx, read_id_idx, run_id_idx]:
                    err = "Invalid header line: {}".format(line)
                    log(err)
                    return False, err
                continue

            # normal line
            output.write("{}\t{}\t{}\n".format(os.path.join(s3_root, reads, line[fast5_idx]),
                                               line[read_id_idx], line[run_id_idx]))
            stored_reads += 1

    # cleanup
    if not keep_workdir:
        shutil.rmtree(run_workdir)

    return True, stored_reads


def s3cmd_ls(location, strip_metadata=True):
    # basically the difference between "ls" and "ls -la"
    if strip_metadata:  postprocess = lambda x: x.strip().split()[-1] if len(x.strip()) > 0 else ""
    else:               postprocess = lambda x: x.strip()

    ls_cmd = ['s3cmd', 'ls', location]
    bytecontents = subprocess.check_output(ls_cmd)
    contents = list(filter(lambda x: len(x) > 0, map(postprocess, bytecontents.decode('ascii', 'ignore').split("\n"))))

    return contents


def get_fast5_datasets(s3_root):
    # get directory
    contents = s3cmd_ls(s3_root, strip_metadata=True)

    # get raw info
    reads = None
    fastqs = list()
    summaries = list()
    other_file_count = 0
    for original_element in contents:
        elem = original_element.replace(s3_root, "")
        if elem.rstrip("/") == "reads": reads = list(map(lambda x: x.replace(s3_root, ""), s3cmd_ls(original_element)))
        elif elem.endswith(".fastq"): fastqs.append(elem)
        elif elem.endswith(".txt") and "summary" in elem: summaries.append(elem)
        else: other_file_count += 1

    # sanity check and report
    assert reads is not None, "Found no 'reads' directory in {}".format(s3_root)
    log("Summary of {}: {} reads, {} fastqs, {} summaries".format(s3_root, len(reads), len(fastqs), len(summaries)))
    assert False not in set(map(lambda x: len(x) >= 2, [reads, fastqs, summaries])), "Expected two or more elements!"

    # organize methods
    def calculate_prefix(a, b):
        for idx, char in zip(range(min(len(a), len(b))), a):
            if b[idx] != char: return a[:idx]
        return a[:idx+1]
    def calculate_suffix(a, b):
        if a[-1] != b[-1]: return ""
        for idx in map(lambda x: -1 * (x + 1), range(min(len(a), len(b)))):
            if b[idx] != a[idx]: return a[idx+1:]
        return a[idx:]
    def calculate_ixes(a, b):
        return calculate_prefix(a, b), calculate_suffix(a, b)

    # get prefix, suffix
    read_prefix, read_suffix = calculate_ixes(reads[0], reads[1])
    fastq_prefix, fastq_suffix = calculate_ixes(fastqs[0], fastqs[1])
    summary_prefix, summary_suffix = calculate_ixes(summaries[0], summaries[1])

    # build identifier map
    id_map = dict()
    def add_to_id_map(elements, key, prefix, suffix):
        for elem in elements:
            id = int(elem.replace(prefix, "").replace(suffix, ""))
            map = {KEY_ID: id, KEY_S3_ROOT: s3_root} if id not in id_map else id_map[id]
            map[key] = elem
            id_map[id] = map
    add_to_id_map(reads, KEY_READS, read_prefix, read_suffix)
    add_to_id_map(fastqs, KEY_FASTQ, fastq_prefix, fastq_suffix)
    add_to_id_map(summaries, KEY_SUMMARY, summary_prefix, summary_suffix)

    # organize
    fast5_dataset = list()
    has_missing_elements = list()
    for elem in id_map.values():
        if len(elem) != 5:
            has_missing_elements.append(elem)
        else:
            fast5_dataset.append(elem)

    # sanity check
    if len(has_missing_elements) > 0:
        log("ERROR data missing for {} fast5 datasets:".format(len(has_missing_elements)))
        for melem in has_missing_elements:
            log("\t{}: {}".format(melem[KEY_ID], melem))
    log("Found all info for {} fast5 datasets".format(len(fast5_dataset)))

    # return
    return fast5_dataset


def main():
    args = parse_args()
    if not args.s3_root.endswith("/"): args.s3_root = "{}/".format(args.s3_root)
    if not os.path.isdir(args.index_destination): os.mkdir(args.index_destination)
    if not os.path.isdir(args.workdir): os.mkdir(args.workdir)

    # loggit
    log("Running S3 FAST5 indexing program:")
    log("\ts3_root:     {}".format(args.s3_root))
    log("\tworkdir:     {}".format(args.workdir))
    log("\tdestination: {}".format(args.index_destination))
    log("")

    # get data from s3
    datasets = get_fast5_datasets(args.s3_root)

    # service
    organization_args = {
        'workdir': args.workdir,
        'destination_dir': args.index_destination
    }
    total, failures, messages = run_service(fast5_promethion_s3_organization_service, datasets, organization_args,
                                            "fast5_data_locations", args.threads)

    # loggit
    log("Indexed {} fast5 datasets with {} failures".format(total, failures))
    total_read_count = 0
    for msg in messages:
        if msg.startswith(KEY_READ_COUNT):
            total_read_count += int(msg.split(":")[1])
    log("Total indexed reads: {}".format(total_read_count))


if __name__ == "__main__":
    main()