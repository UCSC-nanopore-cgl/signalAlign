#!/usr/bin/env python
"""Train HMMs for alignment of signal data from the MinION
"""

import sys
import os
import urllib.parse
import textwrap
import yaml
import h5py

from argparse import ArgumentParser
from random import shuffle
from shutil import copyfile

from multiprocessing import Process, current_process, Manager

from signalalign import parseFofn, DEFAULT_TRAINMODELS_OPTIONS
from signalalign.signalAlignment import multithread_signal_alignment
from signalalign.hiddenMarkovModel import ContinuousPairHmm, HdpSignalHmm
from signalalign.utils import processReferenceFasta
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.bwaWrapper import buildBwaIndex


class AbstractSamples(object):
    def __init__(self, source, fw_fasta_path, bw_fasta_path):
        self.source = source
        self.fw_fasta_path = fw_fasta_path
        self.bw_fasta_path = bw_fasta_path

    def _parse(self):
        raise NotImplementedError

    def getFiles(self):
        raise NotImplementedError

    def getKey(self):
        return self.source

    def getReferences(self):
        return self.fw_fasta_path, self.bw_fasta_path


class Fast5Directory(AbstractSamples):
    def __init__(self, source, fw_fasta_path, bw_fasta_path):
        AbstractSamples.__init__(self, source, fw_fasta_path, bw_fasta_path)
        self.files = self._parse()
        self.fw_fasta_path = fw_fasta_path
        self.bw_fasta_path = bw_fasta_path

    def _parse(self):
        return [os.path.abspath(os.path.join(self.source, x)) for x in os.listdir(self.source) if x.endswith(".fast5")]

    def getFiles(self):
        return self.files


class FileOfFilenames(AbstractSamples):
    def __init__(self, source, fw_fasta_path, bw_fasta_path):
        AbstractSamples.__init__(self, source, fw_fasta_path, bw_fasta_path)
        self.files = self._parse()
        self.fw_fasta_path = fw_fasta_path
        self.bw_fasta_path = bw_fasta_path

    def _parse(self):
        return parseFofn(self.source)

    def getFiles(self):
        return self.files


def parse_args():
    parser = ArgumentParser(description=__doc__)
    # required arguments
    parser.add_argument('--file_directory', '-d', action='append', default=None,
                        dest='files_dir', required=True, type=str,
                        help="directories with fast5 files to train on. example: ../reads/")
    parser.add_argument('--ref', '-r', action='store', default=None, dest='ref', required=True, type=str,
                        help="location of refrerence sequence in FASTA, example: ../ref.fasta")
    parser.add_argument('--output_location', '-o', action='store', dest='out', default=None,
                        required=True, type=str,
                        help="directory to put the trained model, and use for working directory. example: ./scratch/")
    # optional arguments
    parser.add_argument("--2d", action='store_true', dest="twoD", default=False, help="flag, reads are 2D chemistry.")
    parser.add_argument("--bwt", action='store', dest="bwt", default=None,
                        help="path to BWT files. example: ../ref.fasta")
    parser.add_argument('--stateMachineType', '-smt', action='store', dest='stateMachineType', type=str,
                        default="threeState", required=False,
                        help="StateMachine options: threeState, threeStateHdp")
    parser.add_argument("--file_of_files", "-fofn", action="append", required=False, default=None, dest="fofn",
                        type=str, help="text file with absolute paths of files to use")
    parser.add_argument('--iterations', '-i', action='store', dest='iter', default=10,
                        required=False, type=int, help='number of iterations to perform')
    parser.add_argument('--train_amount', '-a', action='store', dest='amount', default=15000,
                        required=False, type=int,
                        help="limit the total length of sequence to use in training (batch size).")
    parser.add_argument('--in_template_hmm', '-T', action='store', dest='in_T_Hmm',
                        required=True, type=str, help="template model to bootstrap from, find a starting model in the "
                                                      "models directory")
    parser.add_argument('--in_complement_hmm', '-C', action='store', dest='in_C_Hmm',
                        required=True, type=str,
                        help="complement model to bootstrap from, find a starting model in the "
                             "models directory")
    parser.add_argument('--templateHDP', '-tH', action='store', dest='templateHDP', default=None,
                        help="path to template HDP model to use")
    parser.add_argument('--complementHDP', '-cH', action='store', dest='complementHDP', default=None,
                        help="path to complement HDP model to use")
    parser.add_argument('--jobs', '-j', action='store', dest='nb_jobs', required=False, default=4,
                        type=int, help="number of jobs to run concurrently")
    parser.add_argument('--test', action='store_true', default=False, dest='test', help="Used for CI testing")
    parser.add_argument('--ambiguity_positions', '-p', action='store', required=False, default=None,
                        dest='substitution_file', help="Substitution positions")
    parser.add_argument("--motif", action="store", dest="motif_key", default=None)
    parser.add_argument('--ambig_char', '-X', action='append', required=False, default=None, type=str, dest='labels',
                        help="Character to substitute at positions, default is 'X'.")
    parser.add_argument('--diagonalExpansion', '-e', action='store', dest='diag_expansion', type=int,
                        required=False, default=None,
                        help="number of diagonals to expand around each anchor default: 50")
    parser.add_argument('--constraintTrim', '-m', action='store', dest='constraint_trim', type=int,
                        required=False, default=None, help='amount to remove from an anchor constraint')
    parser.add_argument('--debug', action='store_true', dest="DEBUG", default=False)

    args = parser.parse_args()
    return args


def get_2d_length(fast5):
    read = h5py.File(fast5, 'r')
    read_length = 0
    twoD_read_sequence_address = "/Analyses/Basecall_2D_000/BaseCalled_2D/Fastq"
    if not (twoD_read_sequence_address in read):
        print("This read didn't have a 2D read", fast5, end='\n', file=sys.stderr)
        read.close()
        return 0
    else:
        read_length = len(read[twoD_read_sequence_address][()].split()[2])
        read.close()
        return read_length


def get_1d_length(fast5):
    print(fast5)
    read = h5py.File(fast5, "r")
    read_length = 0
    template_fastq_address = "/Analyses/Basecall_1D_000/BaseCalled_template/Fastq"
    if not (template_fastq_address in read):
        print("Read %s has not been basecalled" % fast5)
        read.close()
        return 0
    else:
        read_length = len(read[template_fastq_address][()].split()[2])
        print("read %s has %s bases" % (fast5, read_length))
        read.close()
        return read_length


def cull_training_files(samples, training_amount, twoD):
    print("trainModels - culling training files.\n", end="", file=sys.stderr)
    # training_files = []
    for sample in samples:
        shuffle(sample.getFiles())
        total_amount = 0
        file_count = 0
        get_seq_len_fcn = get_2d_length if twoD else get_1d_length
        # loop over files and add them to training list, break when we have enough bases to complete a batch
        # make a list of tuples [(fast5_path, (plus_ref_seq, minus_ref_seq))]
        for f in sample.getFiles():
            yield f, sample.fw_fasta_path, sample.bw_fasta_path
            # training_files.append((f, sample.fw_fasta_path, sample.bw_fasta_path))
            file_count += 1
            total_amount += get_seq_len_fcn(f)
            if total_amount >= training_amount:
                break
        print("Culled {file_count} training files, for {bases} from {sample}."
              .format(file_count=file_count, bases=total_amount, sample=sample.getKey()), end="\n", file=sys.stderr)

    # shuffle(training_files)
    # return training_files  # [(path_to_fast5, reference_map)...]


def get_model(model_type, model_file):
    assert (model_type in ["threeState", "threeStateHdp"]), "Unsupported StateMachine type"
    # todo clean this up
    assert model_file is not None, "Need to have starting lookup table for {} HMM".format(model_type)
    if model_type == "threeState":
        model = ContinuousPairHmm(model_type=model_type)
        model.load_model(model_file=model_file)
        return model
    if model_type == "threeStateHdp":
        model = HdpSignalHmm(model_type=model_type)
        model.load_model(model_file=model_file)
        return model


def add_and_norm_expectations(path, files, model, hmm_file, update_transitions=False, update_emissions=False):
    if update_emissions is False and update_transitions is False:
        print("[trainModels] NOTICE: Training transitions by default\n", file=sys.stderr)
        update_transitions = True

    model.likelihood = 0
    files_added_successfully = 0
    files_with_problems = 0
    for f in files:
        try:
            success = model.add_expectations_file(path + f)
            os.remove(path + f)
            if success:
                files_added_successfully += 1
            else:
                files_with_problems += 1
        except Exception as e:
            print("Problem adding expectations file {file} got error {e}".format(file=path + f, e=e),
                  file=sys.stderr)
            os.remove(path + f)
            files_with_problems += 1
    model.normalize(update_transitions=update_transitions, update_emissions=update_emissions)
    model.write(hmm_file)
    model.running_likelihoods.append(model.likelihood)
    if type(model) is HdpSignalHmm:
        model.reset_assignments()
    print("[trainModels] NOTICE: Added {success} expectations files successfully, {problem} files had problems\n"
          "".format(success=files_added_successfully, problem=files_with_problems), file=sys.stderr)


def build_hdp(template_hdp_path, complement_hdp_path, template_assignments, complement_assignments, samples,
              burn_in, thinning, verbose=False):
    assert (template_assignments is not None) and (complement_assignments is not None), \
        "trainModels - ERROR: missing assignments"

    if verbose is True:
        verbose_flag = "--verbose "
    else:
        verbose_flag = ""

    command = "./buildHdpUtil {verbose}-v {tHdpP} -w {cHdpP} -E {tExpectations} -W {cExpectations} " \
              "-n {samples} -I {burnIn} -t {thinning}".format(tHdpP=template_hdp_path,
                                                              cHdpP=complement_hdp_path,
                                                              tExpectations=template_assignments,
                                                              cExpectations=complement_assignments,
                                                              samples=samples, burnIn=burn_in,
                                                              thinning=thinning,
                                                              verbose=verbose_flag)
    print("[trainModels] Running command:{}".format(command), file=sys.stderr)
    os.system(command)  # todo try checkoutput
    print("trainModels - built HDP.", file=sys.stderr)
    return


def generateConfig(config_path):
    if os.path.exists(config_path):
        raise RuntimeError
    config_content = textwrap.dedent("""\
                # SignalAlign model training config file
                output_dir: signalAlign_unittest/
                samples: [{
                    fast5_dir: ../tests/minion_test_reads/canonical_ecoli_R9/,
                    fofn:,
                    positions_file:,
                    motif:,
                    label:,
                }]
                bwa_reference: ../tests/test_sequences/E.coli_k12.fasta
                bwt:
                stateMachineType: threeState
                in_T_Hmm: ../models/testModelR9_5mer_acegot_template.model
                in_C_Hmm: ../models/testModelR9_5mer_acegot_complement.model
                templateHdp:
                complementHdp:
                iterations: 3
                training_bases: 10000
                job_count: 4
                diagonal_expansion:
                constraint_trim:
                twoD: true
                alignment_file:
                DEBUG:
                test: true
                """)
    fH = open(config_path, "w")
    fH.write(config_content)
    fH.flush()
    fH.close()


def process_sample(sample, reference, working_folder):
    """Process a set of Fast5 files. Creates edited reference sequences if needed

    :param samples: dictionary of a sample of fast5's
            fast5_dir:  path to fast5s,
            fofn: file with paths to fast5s
            positions_file: changes to nucleotides by position,
            edited_fw_reference: forward reference if different than canonical
            edited_bw_reference: backward reference if needed

    :param reference: reference sequence if needed to process
    :param working_folder: FolderHandler object with a folder already created.
    """

    options = dict(**DEFAULT_TRAINMODELS_OPTIONS)
    options.update(sample)
    if options["fast5_dir"] is None and options["fofn"] is None:
        raise RuntimeError("Need to provide path to .fast5 files or file with filenames (fofn)")

    if options["edited_fw_reference"] is None:
        assert os.path.isfile(reference), "Must specify a bwa_reference in order to create signalAlignments. {}" \
                                          "".format(reference)
        fw_fasta_path, bw_fasta_path = processReferenceFasta(fasta=reference,
                                                             work_folder=working_folder,
                                                             motif_key=options["motif"],
                                                             sub_char=options["label"],
                                                             positions_file=options["positions_file"])

    else:
        fw_fasta_path = options["edited_fw_reference"]
        bw_fasta_path = options["edited_bw_reference"]

    if options["fast5_dir"] is not None:
        if options["fofn"] is not None:
            print("WARNING Only using files is directory %s ignoring fofn %s"
                  % (options["files_dir"], options["fofn"]))
        sample = Fast5Directory(options["fast5_dir"], fw_fasta_path, bw_fasta_path)
    else:
        sample = FileOfFilenames(options["fofn"], fw_fasta_path, bw_fasta_path)
    return sample


def trainHMMTransitions(config):
    """Train HMM transitions.

    :param config: config dictionary with required arguments

    config required keys
    _____________
    output_dir : path to output directory
    samples: list of dictionary of samples
            fast5_dir:  path to fast5s,
            fofn: file with paths to fast5s
            positions_file: changes to nucleotides by position,
            edited_fw_reference: forward reference if different than canonical
            edited_bw_reference: backward reference if needed

    reference: reference to index if bwt does not exist
    bwa_reference : original reference to create guide alignments
    in_T_Hmm: path to template HMM model file
    in_C_Hmm: path to complement HMM model file
    templateHdp: path to template HDP model file
    complementHdp: path to complement HDP model file
    twoD: boolean option for 2D reads
    training_bases: number of bases to use for each sample during training
    job_count: number of processes to use when running SignalAlign
    iterations: number of iterations over the entire updating pipeline
    stateMachineType: type of stateMachine ("threeStateHdp" or "threeStateHmm")
    diagonal_expansion: alignment algorithm param to expand how far a path can get from guide alignment
    constraint_trim: alignment algorithm param for how much to trim the guide alignment anchors
    """
    # make directory to put the files we're using
    output_dir = config["output_dir"]
    bwa_reference = config["bwa_reference"]
    alignment_file = config["alignment_file"]
    training_amount = config["training_bases"]
    workers = config["job_count"]
    iterations = config["iterations"]
    twoD = config['twoD']
    state_machine_type = config["stateMachineType"]
    template_model_path = config["in_T_Hmm"]
    complement_model_path = config["in_C_Hmm"]
    original_template_hdp_path = config["templateHdp"]
    original_complement_hdp_path = config["complementHdp"]
    diagonal_expansion = config["diagonal_expansion"]
    constraint_trim = config["constraint_trim"]
    test = config["test"]

    working_folder = FolderHandler()
    working_folder_path = working_folder.open_folder(os.path.join(output_dir, "temp_trainModels"))
    print(working_folder_path)

    # find model files
    # make some paths to files to hold the HMMs
    complement_model = None
    complement_hmm = None

    if twoD:
        assert os.path.exists(complement_model_path), \
            "Missing complement model %s" % (complement_model_path)
        complement_model_path = os.path.abspath(complement_model_path)
        complement_model = get_model(state_machine_type, complement_model_path)
        complement_hmm = working_folder.add_file_path("complement_trained.hmm")
        copyfile(complement_model_path, complement_hmm)
        assert os.path.exists(complement_hmm), "Problem copying default model to {}".format(complement_hmm)

    assert os.path.exists(template_model_path), \
        "Missing template model %s" % (template_model_path)
    template_model_path = os.path.abspath(template_model_path)
    template_model = get_model(state_machine_type, template_model_path)
    template_hmm = working_folder.add_file_path("template_trained.hmm")
    copyfile(os.path.abspath(template_model_path), template_hmm)
    assert os.path.exists(template_hmm), "Problem copying default model to {}".format(template_hmm)
    # determine if we train HMM emissions or HDP
    update_template_hmm_emissions = False
    update_complement_hmm_emissions = False
    # get the input HDP, if we're using it
    if state_machine_type == "threeStateHdp":
        assert os.path.exists(original_template_hdp_path), "Templace HDP path not found {}".format(original_template_hdp_path)
        original_template_hdp_path = os.path.abspath(original_template_hdp_path)
        template_hdp = working_folder.add_file_path("%s" % os.path.basename(original_template_hdp_path))
        copyfile(original_template_hdp_path, template_hdp)
        if twoD:
            assert os.path.exists(original_complement_hdp_path), "Templace HDP path not found {}".format(original_complement_hdp_path)
            original_complement_hdp_path = os.path.abspath(original_complement_hdp_path)

            complement_hdp = working_folder.add_file_path("%s" % os.path.basename(original_complement_hdp_path))
            copyfile(original_complement_hdp_path, complement_hdp)
        else:
            update_complement_hmm_emissions = True
            complement_hdp = None
    else:
        update_template_hmm_emissions = True
        update_complement_hmm_emissions = True
        template_hdp = None
        complement_hdp = None

    # start iterating
    if bwa_reference:
        print("BWA_REFERENCE")
        bwa_reference = os.path.abspath(bwa_reference)
        print("os.path.isfile(bwa_reference)", os.path.isfile(bwa_reference))
    print("BWA_REFERENCE222", bwa_reference)

    i = 0
    samples = [process_sample(s, bwa_reference, working_folder) for s in config["samples"]]

    while i < iterations:
        # first cull a set of files to get expectations on
        # training_files = cull_training_files(samples=samples,
        #                                      training_amount=config["training_bases"],
        #                                      twoD=config["twoD"])

        for sample in samples:
            shuffle(sample.getFiles())
            total_amount = 0
            file_count = 0
            get_seq_len_fcn = get_2d_length if twoD else get_1d_length
            # loop over files and add them to training list, break when we have enough bases to complete a batch
            # collect paths to fast5 files
            list_of_fast5s = []
            for f in sample.getFiles():
                # yield f, sample.fw_fasta_path, sample.bw_fasta_path
                list_of_fast5s.append(f)
                # training_files.append((f, sample.fw_fasta_path, sample.bw_fasta_path))
                file_count += 1
                total_amount += get_seq_len_fcn(f)
                if total_amount >= training_amount:
                    break
                print("[trainModels_HMM] Culled {file_count} training files, for {bases} from {sample}."
                      .format(file_count=file_count, bases=total_amount, sample=sample.getKey()), end="\n", file=sys.stderr)
            if alignment_file:
                alignment_file = os.path.abspath(alignment_file)

            alignment_args = {
                "backward_reference": sample.bw_fasta_path,
                "forward_reference": sample.fw_fasta_path,
                "destination": working_folder_path,
                "stateMachineType": state_machine_type,
                "in_templateHmm": template_hmm,
                "in_complementHmm": complement_hmm,
                "in_templateHdp": template_hdp,
                "in_complementHdp": complement_hdp,
                "threshold": 0.01,
                "diagonal_expansion": diagonal_expansion,
                "constraint_trim": constraint_trim,
                "target_regions": None,
                "degenerate": None,
                "twoD_chemistry": twoD,
                "alignment_file": alignment_file,
                "bwa_reference": bwa_reference,
                'track_memory_usage': False,
                'get_expectations': True
            }
            multithread_signal_alignment(alignment_args, list_of_fast5s, workers)


        # load then normalize the expectations
        template_expectations_files = [x for x in os.listdir(working_folder_path)
                                       if x.endswith(".template.expectations")]

        complement_expectations_files = [x for x in os.listdir(working_folder_path)
                                         if x.endswith(".complement.expectations")]

        if len(template_expectations_files) > 0:
            add_and_norm_expectations(path=working_folder_path,
                                      files=template_expectations_files,
                                      model=template_model,
                                      hmm_file=template_hmm,
                                      update_transitions=True,
                                      update_emissions=False)
        if twoD and len(complement_expectations_files) > 0:
            add_and_norm_expectations(path=working_folder_path,
                                      files=complement_expectations_files,
                                      model=complement_model,
                                      hmm_file=complement_hmm,
                                      update_transitions=True,
                                      update_emissions=False)

        # log the running likelihood
        if len(template_model.running_likelihoods) > 0 and \
                (twoD and len(complement_model.running_likelihoods)) > 0:
            print("[trainModels_HMM] {i}| {t_likelihood}\t{c_likelihood}".format(t_likelihood=template_model.running_likelihoods[-1],
                                                               c_likelihood=complement_model.running_likelihoods[-1],
                                                               i=i))
            if test and (len(template_model.running_likelihoods) >= 2) and \
                    (config["twoD"] and len(complement_model.running_likelihoods) >= 2):
                assert (template_model.running_likelihoods[-2] < template_model.running_likelihoods[-1]) and \
                       (complement_model.running_likelihoods[-2] < complement_model.running_likelihoods[-1]), \
                    "Testing: Likelihood error, went up"
        elif len(template_model.running_likelihoods) > 0:
            print("[trainModels_HMM] {i}| {t_likelihood}".format(t_likelihood=template_model.running_likelihoods[-1], i=i))
            if test and (len(template_model.running_likelihoods) >= 2):
                assert (template_model.running_likelihoods[-2] < template_model.running_likelihoods[-1]), "Testing: Likelihood error, went up"

        i += 1

    # if we're using HDP, trim the final Hmm (remove assignments)

    print("trainModels - finished training routine", file=sys.stdout)
    print("trainModels - finished training routine", file=sys.stderr)

    return [template_hmm, complement_hmm, template_hdp, complement_hdp]


def main():
    def parse_args():
        parser = ArgumentParser()
        subparsers = parser.add_subparsers(dest="command")

        # parsers for running the full pipeline
        run_parser = subparsers.add_parser("run", help="runs full workflow ")
        run_parser.add_argument('--config', default='trainModels-config.yaml', type=str,
                                help='Path to the (filled in) config file, generated with "generate".')
        subparsers.add_parser("generate", help="generates a config file for your run, do this first")
        return parser.parse_args()

    args = parse_args()
    if args.command == "generate":
        try:
            config_path = os.path.join(os.getcwd(), "trainModels-config.yaml")
            generateConfig(config_path)
        except RuntimeError:
            print("Using existing config file {}".format(config_path))
            pass
    elif args.command == "run":
        if not os.path.exists(args.config):
            print("{config} not found run generate-config".format(config=args.config))
            exit(1)
        # Parse config
        config = {x.replace('-', '_'): y for x, y in yaml.load(open(args.config).read()).items()}
        trainHMMTransitions(config)


if __name__ == "__main__":
    sys.exit(main())
