#!/usr/bin/env python3

import os
import sys
import csv
import numpy as np
import pysam
# from sonLib.bioio import fastaWrite
from signalalign import defaultModelFromVersion
from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.bwaWrapper import generateGuideAlignment
from signalalign.utils.fileHandlers import FolderHandler
from signalalign.utils.sequenceTools import fastaWrite
from signalalign.mea_algorithm import mea_alignment_from_signal_align, match_events_with_signalalign


class SignalAlignment(object):
    def __init__(self,
                 in_fast5,
                 destination,
                 stateMachineType,
                 bwa_index,
                 in_templateHmm,
                 in_complementHmm,
                 in_templateHdp,
                 in_complementHdp,
                 threshold,
                 diagonal_expansion,
                 constraint_trim,
                 degenerate,
                 twoD_chemistry,
                 forward_reference,
                 backward_reference=None,
                 target_regions=None,
                 output_format="full",
                 embed=False,
                 event_table=False):
        self.in_fast5 = in_fast5  # fast5 file to align
        self.destination = destination  # place where the alignments go, should already exist
        self.stateMachineType = stateMachineType  # flag for signalMachine
        self.bwa_index = bwa_index  # index of reference sequence
        self.threshold = threshold  # min posterior probability to keep
        self.diagonal_expansion = diagonal_expansion  # alignment algorithm param
        self.constraint_trim = constraint_trim  # alignment algorithm param
        self.output_format = output_format  # smaller output files
        self.degenerate = degenerate  # set of nucleotides for degenerate characters
        self.twoD_chemistry = twoD_chemistry  # flag for 2D sequencing runs
        self.temp_folder = FolderHandler()  # object for holding temporary files (non-toil)
        self.read_name = self.in_fast5.split("/")[-1][:-6]  # get the name without the '.fast5'
        self.target_regions = target_regions
        self.output_formats = {"full": 0, "variantCaller": 1, "assignments": 2}
        self.embed = embed  # embed the output into the fast5 file
        self.event_table = event_table  # specify which event table to use to generate alignments
        self.backward_reference = backward_reference  # fasta path to backward reference if modified bases are used
        self.forward_reference = forward_reference  # fasta path to forward reference


        if (in_templateHmm is not None) and os.path.isfile(in_templateHmm):
            self.in_templateHmm = in_templateHmm
        else:
            self.in_templateHmm = None
        if (in_complementHmm is not None) and os.path.isfile(in_complementHmm):
            self.in_complementHmm = in_complementHmm
        else:
            self.in_complementHmm = None

        # similarly for HDPs
        if (in_templateHdp is not None) and os.path.isfile(in_templateHdp):
            self.in_templateHdp = in_templateHdp
        else:
            self.in_templateHdp = None
        if (in_complementHdp is not None) and os.path.isfile(in_complementHdp):
            self.in_complementHdp = in_complementHdp
        else:
            self.in_complementHdp = None

    def run(self, get_expectations=False):
        print("[SignalAlignment.run]INFO: Starting on {read}".format(read=self.in_fast5), file=sys.stderr)
        if get_expectations:
            assert self.in_templateHmm is not None and self.in_complementHmm is not None, \
                "Need HMM files for model training"
        # file checks
        if os.path.isfile(self.in_fast5) is False:
            print("[SignalAlignment.run]ERROR: Did not find .fast5 at{file}".format(file=self.in_fast5))
            return False

        self.openTempFolder("tempFiles_%s" % self.read_name)
        npRead_ = self.addTempFilePath("temp_%s.npRead" % self.read_name)
        # TODO is this totally fucked for RNA because of 3'-5' mapping?
        npRead = NanoporeRead(fast_five_file=self.in_fast5, twoD=self.twoD_chemistry, event_table=self.event_table)
        fH = open(npRead_, "w")
        ok = npRead.Write(parent_job=None, out_file=fH, initialize=True)
        fH.close()
        if not ok:
            self.failStop("[SignalAlignment.run]File: %s did not pass initial checks" % self.read_name, npRead)
            return False

        read_label = npRead.read_label  # use this to identify the read throughout
        read_fasta_ = self.addTempFilePath("temp_seq_%s.fa" % read_label)
        temp_samfile_ = self.addTempFilePath("temp_sam_file_%s.sam" % read_label)
        cigar_file_ = self.addTempFilePath("temp_cigar_%s.txt" % read_label)
        if self.twoD_chemistry:
            ok, version, pop1_complement = self.prepare_twod(nanopore_read=npRead, twod_read_path=read_fasta_)
        else:
            ok, version, _ = self.prepare_oned(nanopore_read=npRead, oned_read_path=read_fasta_)
            pop1_complement = None
        # add an indicator for the model being used
        if self.stateMachineType == "threeState":
            model_label = ".sm"
            stateMachineType_flag = ""
        elif self.stateMachineType == "threeStateHdp":
            model_label = ".sm3Hdp"
            stateMachineType_flag = "--sm3Hdp "
            if self.twoD_chemistry:
                assert (self.in_templateHdp is not None) and (self.in_complementHdp is not None), "Need to provide HDPs"
            else:
                assert self.in_templateHdp is not None, "Need to provide Template HDP"
        else:  # make invalid stateMachine control?
            model_label = ".sm"
            stateMachineType_flag = ""

        guide_alignment = generateGuideAlignment(bwa_index=self.bwa_index, query=read_fasta_,
                                                 temp_sam_path=temp_samfile_,
                                                 target_regions=self.target_regions)
        # ok = guide_alignment.validate(list(self.reference_map.keys()))
        ok = guide_alignment.validate()

        if not ok:
            self.failStop("[SignalAlignment.run]ERROR getting guide alignment", npRead)
            return False

        cig_handle = open(cigar_file_, "w")
        cig_handle.write(guide_alignment.cigar + "\n")
        cig_handle.close()

        # next section makes the output file name with the format: /directory/for/files/file.model.orientation.tsv
        posteriors_file_path = ''
        # forward strand
        if guide_alignment.strand == "+":
            if self.output_format == "full":
                posteriors_file_path = self.destination + read_label + model_label + ".forward.tsv"
            elif self.output_format == "variantCaller":
                posteriors_file_path = self.destination + read_label + model_label + ".tsv"
            else:
                posteriors_file_path = self.destination + read_label + model_label + ".assignments"

        # backward strand
        if guide_alignment.strand == "-":
            if self.output_format == "full":
                posteriors_file_path = self.destination + read_label + model_label + ".backward.tsv"
            elif self.output_format == "variantCaller":
                posteriors_file_path = self.destination + read_label + model_label + ".tsv"
            else:
                posteriors_file_path = self.destination + read_label + model_label + ".assignments"

        # Alignment/Expectations routine
        path_to_signalAlign = "./signalMachine"

        # flags

        # input (match) models
        if self.in_templateHmm is None:
            self.in_templateHmm = defaultModelFromVersion(strand="template", version=version)
        if self.twoD_chemistry:
            if self.in_complementHmm is None:
                self.in_complementHmm = defaultModelFromVersion(strand="complement", version=version,
                                                                pop1_complement=pop1_complement)

        assert self.in_templateHmm is not None
        if self.twoD_chemistry:
            if self.in_complementHmm is None:
                self.failStop("[SignalAlignment.run]ERROR Need to have complement HMM for 2D analysis", npRead)
                return False

        template_model_flag = "-T {} ".format(self.in_templateHmm)
        if self.twoD_chemistry:
            complement_model_flag = "-C {} ".format(self.in_complementHmm)
        else:
            complement_model_flag = ""

        print("[SignalALignment.run]NOTICE: template model {t} complement model {c}"
              "".format(t=self.in_templateHmm, c=self.in_complementHmm), file=sys.stderr)

        # reference sequences
        assert os.path.isfile(self.forward_reference)
        forward_ref_flag = "-f {f_ref} ".format(f_ref=self.forward_reference)
        if self.backward_reference:
            assert os.path.isfile(self.backward_reference)
            backward_ref_flag = "-b {b_ref} ".format(b_ref=self.backward_reference)
        else:
            backward_ref_flag = ""

        # input HDPs
        if (self.in_templateHdp is not None) or (self.in_complementHdp is not None):
            hdp_flags = "-v {tHdp_loc} ".format(tHdp_loc=self.in_templateHdp)
            if self.twoD_chemistry and self.in_complementHdp is not None:
                hdp_flags += "-w {cHdp_loc} ".format(cHdp_loc=self.in_complementHdp)
        else:
            hdp_flags = ""

        # threshold
        if self.threshold is not None:
            threshold_flag = "-D {threshold} ".format(threshold=self.threshold)
        else:
            threshold_flag = ""

        # diagonal expansion
        if self.diagonal_expansion is not None:
            diag_expansion_flag = "-x {expansion} ".format(expansion=self.diagonal_expansion)
        else:
            diag_expansion_flag = ""

        # constraint trim
        if self.constraint_trim is not None:
            trim_flag = "-m {trim} ".format(trim=self.constraint_trim)
        else:
            trim_flag = ""

        # output format
        if self.output_format not in list(self.output_formats.keys()):
            self.failStop("[SignalAlignment.run]ERROR illegal output format selected %s" % self.output_format)
            return False
        out_fmt = "-s {fmt} ".format(fmt=self.output_formats[self.output_format])

        # degenerate nucleotide information
        if self.degenerate is not None:
            degenerate_flag = "-o {} ".format(self.degenerate)
        else:
            degenerate_flag = ""

        if self.twoD_chemistry:
            twoD_flag = "--twoD"
        else:
            twoD_flag = ""
        # commands
        if get_expectations:
            template_expectations_file_path = self.destination + read_label + ".template.expectations"
            complement_expectations_file_path = self.destination + read_label + ".complement.expectations"

            command = \
                "{vA} {td} {degen}{sparse}{model} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} {hdp}-L {readLabel} -p {cigarFile} " \
                "-t {templateExpectations} -c {complementExpectations} -n {seq_name} {f_ref_fa} {b_ref_fa}" \
                    .format(vA=path_to_signalAlign, model=stateMachineType_flag,
                            cigarFile=cigar_file_,
                            npRead=npRead_, readLabel=read_label, td=twoD_flag,
                            templateExpectations=template_expectations_file_path, hdp=hdp_flags,
                            complementExpectations=complement_expectations_file_path, t_model=template_model_flag,
                            c_model=complement_model_flag, thresh=threshold_flag, expansion=diag_expansion_flag,
                            trim=trim_flag, degen=degenerate_flag, sparse=out_fmt, seq_name=guide_alignment.reference_name,
                            f_ref_fa=forward_ref_flag, b_ref_fa=backward_ref_flag)
        else:
            command = \
                "{vA} {td} {degen}{sparse}{model} -q {npRead} " \
                "{t_model}{c_model}{thresh}{expansion}{trim} -p {cigarFile} " \
                "-u {posteriors} {hdp}-L {readLabel} -n {seq_name} {f_ref_fa} {b_ref_fa}" \
                    .format(vA=path_to_signalAlign, model=stateMachineType_flag, sparse=out_fmt,
                            cigarFile=cigar_file_,
                            readLabel=read_label, npRead=npRead_, td=twoD_flag,
                            t_model=template_model_flag, c_model=complement_model_flag,
                            posteriors=posteriors_file_path, thresh=threshold_flag, expansion=diag_expansion_flag,
                            trim=trim_flag, hdp=hdp_flags, degen=degenerate_flag, seq_name=guide_alignment.reference_name,
                            f_ref_fa=forward_ref_flag, b_ref_fa=backward_ref_flag)


        # run
        print("signalAlign - running command: ", command, end="\n", file=sys.stderr)
        os.system(command)
        if self.embed:
            print("signalAlign - embedding into Fast5 ", file=sys.stderr)

            data = self.read_in_signal_align_tsv(posteriors_file_path, file_type=self.output_format)
            npRead = NanoporeRead(fast_five_file=self.in_fast5, twoD=self.twoD_chemistry, event_table=self.event_table)
            npRead.Initialize(None)
            signal_align_path = npRead.get_latest_basecall_edition("/Analyses/SignalAlign_00{}", new=False)
            assert signal_align_path, "There is no path in Fast5 file: {}".format("/Analyses/SignalAlign_00{}")
            output_path = npRead._join_path(signal_align_path, self.output_format)
            npRead.write_data(data, output_path)

            # Todo add attributes to signalalign output
            if self.output_format == "full":
                print("signalAlign - writing maximum expected alignment ", file=sys.stderr)
                alignment = mea_alignment_from_signal_align(None, events=data)
                mae_path = npRead._join_path(signal_align_path, "MEA_alignment_labels")
                events = npRead.get_template_events()
                if events:
                    if guide_alignment.strand == "-":
                        minus = True
                    else:
                        minus = False
                    labels = match_events_with_signalalign(sa_events=alignment,
                                                           event_detections=np.asanyarray(npRead.template_events),
                                                           minus=minus,
                                                           rna=npRead.is_read_rna())
                    npRead.write_data(labels, mae_path)
                    sam_string = str()
                    with open(temp_samfile_, 'r') as test:
                        for line in test:
                            sam_string += line
                    sam_path = npRead._join_path(signal_align_path, "sam")
                    # print(sam_string)
                    npRead.write_data(data=sam_string, location=sam_path, compression=None)

        # self.temp_folder.remove_folder()
        return True

    def prepare_oned(self, nanopore_read, oned_read_path):
        try:
            read_file = open(oned_read_path, "w")
            fastaWrite(fileHandleOrFile=read_file,
                       name=nanopore_read.read_label,
                       seq=nanopore_read.template_read)

            version = nanopore_read.version
            read_file.close()

            return True, version, False
        except Exception as e:
            return False, None, False

    def prepare_twod(self, nanopore_read, twod_read_path):
        # check for table to make 'assembled' 2D alignment table fasta with
        if nanopore_read.has2D_alignment_table is False:
            nanopore_read.close()
            return False, None, False
        fasta_handle = open(twod_read_path, "w")
        fastaWrite(fileHandleOrFile=fasta_handle,
                   name=nanopore_read.read_label,
                   seq=nanopore_read.alignment_table_sequence)
        if nanopore_read.complement_model_id == "complement_median68pA_pop1.model":
            pop1_complement = True
        else:
            pop1_complement = False
        version = nanopore_read.version
        fasta_handle.close()
        nanopore_read.close()
        return True, version, pop1_complement

    def openTempFolder(self, temp_dir):
        self.temp_folder.open_folder("%s%s" % (self.destination, temp_dir))

    def addTempFilePath(self, path_to_add):
        return self.temp_folder.add_file_path(path_to_add)

    def failStop(self, message, nanopore_read=None):
        self.temp_folder.remove_folder()
        if nanopore_read is not None:
            nanopore_read.close()
        print(message, file=sys.stderr)

    def read_in_signal_align_tsv(self, tsv_path, file_type):
        """Read in tsv file"""
        assert file_type in ("full", "assignments", "variantCaller")
        with open(tsv_path, 'r') as tsvin:
            if file_type == "full":
                dtype = [('contig', 'S10'), ('reference_index', int),
                         ('reference_kmer', 'S5'), ('read_file', 'S57'),
                         ('strand', 'S1'), ('event_index', int),
                         ('event_mean', float), ('event_noise', float),
                         ('event_duration', float), ('aligned_kmer', 'S5'),
                         ('scaled_mean_current', float), ('scaled_noise', float),
                         ('posterior_probability', float), ('descaled_event_mean', float),
                         ('ont_model_mean', float), ('path_kmer', 'S5')]
            elif file_type == "assignments":
                dtype = [('k-mer', 'S10'), ('read_file', 'S57'),
                         ('descaled_event_mean', float), ('posterior_probability', float)]

            else:
                dtype = [('event_index', int), ('reference_position', int),
                         ('base', 'S6'), ('posterior_probability', float), ('strand', 'S1'),
                         ('forward_mapped', int), ('read_file', 'S57')]

            event_table = np.loadtxt(tsvin, dtype=dtype)

            def remove_field_name(a, name):
                names = list(a.dtype.names)
                if name in names:
                    names.remove(name)
                b = a[names]
                return b

            event_table = remove_field_name(event_table, "read_file")

        return event_table
