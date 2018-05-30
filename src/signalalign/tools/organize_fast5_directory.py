
from argparse import ArgumentParser
import subprocess
import os
import hashlib
import glob
import shutil

from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.multithread import *
from signalalign.tools.fast5_lookup import RUN_NAME, FAST5_LOCATION, READ_ID, FAST5_ROOT


# KEY_ID = "id"
# KEY_S3_ROOT = "s3_root"
# KEY_READS = "reads"
# KEY_SUMMARY = "summary"
# KEY_FASTQ = "fastq"
#
# KEY_READ_COUNT = "count"

FAST5_SRC_LOCATION = "source"


def parse_args():
    parser = ArgumentParser(description="Organize fast5s in directory, optionally add them to s3")

    parser.add_argument('--input_directory', '-i', action='store', dest='input_directory', required=True, type=str,
                        help="file folder with fast5 files to be organized")
    parser.add_argument('--output_directory', '-o', action='store', dest='output_directory', required=True, type=str,
                        help="directory with organized files (files will be moved to here)")
    parser.add_argument('--copy_files', '-c', action='store_true', dest='copy_files', default=False,
                        help="copy files to output_directory (default: move files)")
    parser.add_argument('--output_folder_prefix', '-f', action='store', dest='output_folder_prefix', default="", type=str,
                        help="string to prepend to output folder")
    parser.add_argument('--output_folder_count', '-n', action='store', dest='output_folder_count', default=16, type=int,
                        help="number of output folders")
    parser.add_argument('--s3_upload_directory', '-s', action='store', dest='s3_upload_directory', default=None, type=str,
                        help="if set, files will be uploaded to this s3 directory")
    parser.add_argument('--threads', '-t', action='store', dest='threads', required=False, type=int, default=1,
                        help="number of threads to use while indexing")

    args = parser.parse_args()
    return args


def log(msg):
    print(msg)


def fast5_file_organization_service(work_queue, done_queue, output_dir_count, output_base, output_index_base,
                                    copy_files, service_name="fast5_file_organization"):
    # prep
    total_handled = 0
    failure_count = 0
    total_reads = 0
    name = current_process().name
    index_files = {}

    #catch overall exceptions
    try:

        # each thread and outputdir gets its own index file (start header for all of them)
        for i in range(output_dir_count):
            index_file = "{}{}_{}.tsv".format(output_index_base, i, name)
            write_header = not os.path.isfile(index_file)
            index_files[i] = open(index_file, 'a')
            if write_header:
                index_files[i].write("##{}:{}\n".format(FAST5_ROOT, get_output_directory(output_base, i)))
                index_files[i].write("#{}\t{}\t{}\n".format(FAST5_LOCATION, READ_ID, RUN_NAME))
        assert len(index_files) == output_dir_count, "unexpected count of index files {} (expected {})".format(
            len(index_files), output_dir_count)

        idx = -1
        for f in iter(work_queue.get, 'STOP'):
            # randomly place files in appropriate directory
            idx = (idx + 1) % output_dir_count

            try:
                # sanity check
                assert 0 <= idx < output_dir_count

                # file organization
                destination_dir = get_output_directory(output_base, idx)
                source = f[FAST5_SRC_LOCATION]
                filename = os.path.basename(source)
                destination = os.path.join(destination_dir, filename)
                action = shutil.move
                if copy_files: action = shutil.copy

                # fast5 organization
                read = NanoporeRead(source)
                if not read._initialize_metadata():
                    failure_count += 1
                    continue
                read_id = read.read_label
                run_id = read.run_id
                assert None not in [read_id, run_id], "Missing read or run id for {}".format(source)

                # move or copy the file
                action(source, destination)

                # write the contents to the index
                index_files[idx].write("{}\t{}\t{}\n".format(destination, read_id, run_id))


            except Exception as e:
                # get error and log it
                message = "{}:{}".format(type(e), str(e))
                error = "{} '{}' failed with: {}".format(service_name, current_process().name, message)
                print("[{}] ".format(service_name) + error)
                done_queue.put(error)
                failure_count += 1

            finally:
                # increment total handling
                total_handled += 1

    except Exception as e:
        # get error and log it
        message = "{}:{}".format(type(e), str(e))
        error = "{} '{}' critically failed with: {}".format(service_name, current_process().name, message)
        print("[{}] ".format(service_name) + error)
        done_queue.put(error)

    finally:
        # close all index files
        for index_file in index_files.values():
            if index_file is not None: index_file.close()

        # logging and final reporting
        print("[%s] '%s' completed %d calls with %d failures"
              % (service_name, current_process().name, total_handled, failure_count))
        done_queue.put("{}:{}".format(TOTAL_KEY, total_handled))
        done_queue.put("{}:{}".format(FAILURE_KEY, failure_count))


def get_output_directory(output_base, idx):
    return os.path.abspath("{}{}".format(output_base, idx)) + "/"


def main():
    args = parse_args()
    # sanitize input dir
    assert os.path.isdir(args.input_directory), "input directory {} does not exist"
    args.input_directory = os.path.abspath(args.input_directory)
    if args.input_directory.endswith("/"): args.input_directory = args.input_directory.rstrip("/")
    # sanitize output dir
    if not os.path.isdir(args.output_directory): os.mkdir(args.output_directory)
    args.output_directory = os.path.abspath(args.output_directory)
    if args.output_directory.endswith("/"): args.output_directory = args.output_directory.rstrip("/")
    # sanitize s3 location
    if args.s3_upload_directory is not None and args.s3_upload_directory.endswith("/"):
        args.s3_upload_directory = args.s3_upload_directory.rstrip("/")

    # loggit
    log("Running FAST5 organization program:")
    log("\tinput directory:      {}".format(args.input_directory))
    log("\toutput directory:     {}".format(args.output_directory))
    log("\tcopy files:           {}".format(args.copy_files))
    log("\ts3_upload_directory:  {}".format(args.s3_upload_directory))
    log("\toutput_folder_prefix: {}".format(args.output_folder_prefix))
    log("\toutput_folder_count:  {}".format(args.output_folder_count))
    log("")

    # file locations
    base_output_location = os.path.join(args.output_directory, args.output_folder_prefix)
    file_index_dir = "{}.file_idx".format(args.output_directory)
    file_index_base = "{}/{}".format(file_index_dir, args.output_folder_prefix.replace("/", "."))

    # move and document fast5 files
    input_files = glob.glob(os.path.join(args.input_directory, "*.fast5"))
    if len(input_files) > 0:
        log("Organizing {} files.\n".format(len(input_files)))

        # make necessary output directories
        output_directories = list(map(lambda x: get_output_directory(base_output_location, x),
                                      range(args.output_folder_count)))
        output_directories.append(file_index_dir)
        [os.makedirs(x, exist_ok=True) for x in filter(lambda x: not os.path.isdir(x), output_directories)]

        # service prep
        organization_service_args = {
            'output_dir_count': args.output_folder_count,
            'output_base': base_output_location,
            'output_index_base': file_index_base,
            'copy_files': args.copy_files,
        }

        # run service
        total, failure, messages = run_service(fast5_file_organization_service, input_files, {}, FAST5_SRC_LOCATION,
                                               args.threads, service_arguments=organization_service_args)
    # no files to write
    else:
        log("No input files found")
    log("")


    # upload to s3
    if args.s3_upload_directory is not None:
        # get output folders
        output_folders = glob.glob(get_output_directory(base_output_location, "*"))
        log("Preparing to move {} folders to {}".format(len(output_folders), args.s3_upload_directory))

        # prep for reindexing
        file_index_files = glob.glob("{}*".format(file_index_base))
        s3_index_dir = "{}.s3_idx".format(args.output_directory)
        if not os.path.isdir(s3_index_dir):
            os.mkdir(s3_index_dir)

        # reindex
        log("Reindexing {} files into {}".format(len(file_index_files), s3_index_dir))
        for file_index in file_index_files:
            file_index_basename = os.path.basename(file_index)
            s3_index = os.path.join(s3_index_dir, file_index_basename)

            with open(file_index, 'r') as fin, open(s3_index, 'w') as s3out:
                for line in fin:
                    line = line.replace(args.output_directory, args.s3_upload_directory)
                    s3out.write(line if line.endswith("\n") else (line + "\n"))

        # force the user to do the upload
        log("\nScript to upload files:\n")
        for output_folder in output_folders:
            log("s3cmd put --recursive {}* {}".format(output_folder, output_folder.replace(
                args.output_directory, args.s3_upload_directory)))

    # no s3 uploading
    else:
        log("Not configured to upload files to s3")


    #finished
    log("\nFin.")




if __name__ == "__main__":
    main()