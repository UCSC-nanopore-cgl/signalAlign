from __future__ import print_function
import sys
import os
import re
import numpy as np
from itertools import islice
from signalalign.fast5 import Fast5
from signalalign.utils.sequenceTools import get_full_nucleotide_read_from_alignment
from signalalign.event_detection import load_from_raw, load_from_raw2
from py3helpers.utils import check_numpy_table
from Bio import SeqIO

try:
    from StringIO import StringIO
except ImportError:
    from io import StringIO

TEMPLATE_BASECALL_KEY = Fast5.__default_basecall_1d_analysis__  # "/Analyses/Basecall_1D_00{}"
TWOD_BASECALL_KEY = Fast5.__default_basecall_2d_analysis__  # "/Analyses/Basecall_2D_00{}"
TWOD_BASECALL_KEY_0 = os.path.join(Fast5.__base_analysis__, TWOD_BASECALL_KEY + "_000")  # "/Analyses/Basecall_2D_000"
METADATA_PATH_KEY = Fast5.__tracking_id_path__  # "/UniqueGlobalKey/tracking_id"
READS_KEY = Fast5.__raw_path__  # "/Raw/Reads/"
VERSION_KEY = ("version", "dragonet version", "nanotensor version", "signalAlign version")
SUPPORTED_1D_VERSIONS = ("1.0.1", "1.2.1", "1.2.4", "1.23.0", "1.22.4", "2.1.0", "0.2.0", "0.1.7", "2.3.1", "1.22.2")
SUPPORTED_2D_VERSIONS = ("1.15.0", "1.19.0", "1.20.0", "1.22.2", "1.22.4", "1.23.0")


# promethion read_name: self.fast5['PreviousReadInfo'].attrs['previous_read_id'].decode()


class NanoporeRead(object):
    def __init__(self, fast_five_file, twoD=False, event_table='', initialize=False, path_to_bin="./",
                 alignment_file=None, model_file_location=None, perform_kmer_event_alignment=None,
                 enforce_supported_versions=True, filter_reads=False, aligned_segment=None, rna=None):
        # load the fast5
        self.filename = fast_five_file  # fast5 file path
        self.fastFive = None  # fast5 object
        self.is_open = False  # bool, is the member .fast5 open?
        self.open()  # attempt to open (will set is_open)
        self.read_label = ""  # read label, template read by default
        self.run_id = ""  # run id describes which reads were sequenced together
        self.alignment_table_sequence = ""  # the sequence made by assembling the alignment table
        self.template_events = []  # template event sequence
        self.complement_events = []  # complement event sequence
        self.template_read = ""  # template strand base (from fastq) sequence
        self.complement_read = ""  # complement strand base (from fastq) sequence
        self.template_strand_event_map = []  # map of events to kmers in the 1D template read
        self.complement_strand_event_map = []  # map of events to kmers in the 1D complement read
        self.template_event_map = []  # map of template events to kmers in 2D read
        self.complement_event_map = []  # map of complement events to kmers in 2D read
        self.stay_prob = 0  # TODO, do I need this?
        self.template_model_name = ""  # legacy, for reads base-called by a specific model (template)
        self.complement_model_name = ""  # legacy, for reads base-called by a specific model (complement)
        self.template_scale = 1  # initial values for scaling parameters
        self.template_shift = 1  #
        self.template_drift = 0  # Template Parameters
        self.template_var = 1  #
        self.template_scale_sd = 1  #
        self.template_var_sd = 1  # --------------------------------------
        self.complement_scale = 1  #
        self.complement_shift = 1  # Complement Parameters
        self.complement_drift = 0  #
        self.complement_var = 1  #
        self.complement_scale_sd = 1  #
        self.complement_var_sd = 1  #
        self.event_table = event_table  # if we look for alternative events table
        self.path_to_bin = path_to_bin  # path to bin
        self.alignment_file = alignment_file
        self.fastq_sequence_address = None
        self.model_file_location = model_file_location
        self.initialize_success = None  # set if initialize was attempted
        self.filter_reads = filter_reads
        # perform_kmer_event_alignment: True - always perform, False - never perform, None - perform if required
        self.perform_kmer_event_alignment = perform_kmer_event_alignment
        self.twoD = twoD  # 2D read flag, necessary right now, and the client should know
        # if set, unsupported versions will cause failure
        self.enforce_supported_versions = enforce_supported_versions
        self.aligned_segment = aligned_segment  # pysam aligned_segment object

        if type(self) == NanoporeRead:
            if twoD:
                raise Exception("The 'twoD' initialization flag is deprecated for NanoporeRead.  "
                                "Use 'NanoporeRead2D' object instead.")
            self.twoD = False
        elif type(self) == NanoporeRead2D:
            self.twoD = True

        # determination of RNA read
        self.rna = False
        if self.is_read_rna() or rna:
            self.rna = True
            assert self.twoD is False, "Cannot perform 2D analysis when using RNA data"
        print("[NanoporeRead:open] is this an rna read?: ", self.rna)
        # initialize if appropriate
        if initialize:
            self.Initialize()

    def open(self):
        if self.is_open: return True
        try:
            self.fastFive = Fast5(self.filename, 'r+')
            self.is_open = True
            return True
        except Exception as e:
            self.close()
            self.logError("[NanoporeRead:open] ERROR opening {filename}, {e}".format(filename=self.filename, e=e))
            return False

    def close(self):
        self.is_open = False
        if self.fastFive:
            self.fastFive.close()

    def get_latest_basecall_edition(self, address, new=False):
        if address.startswith(Fast5.__base_analysis__):
            address = address.replace(Fast5.__base_analysis__, "").lstrip("/")

        try:
            if new:
                return self.fastFive.get_analysis_new(address)
            else:
                return self.fastFive.get_analysis_latest(address)
        except (ValueError, IndexError) as e:
            print("[NanoporeRead:get_latest_basecall_edition] could not find {} in {}".format(address, self.filename),
                  file=sys.stderr)
            return False

    def Initialize(self):
        if not self.open():
            return False

        ok = self._initialize_metadata()
        ok &= self._initialize()
        self.initialize_success = ok

        if not ok:
            self.close()

        return ok

    def _initialize_metadata(self):
        if not self.open(): return False

        path_locations = [METADATA_PATH_KEY, READS_KEY]
        missing_paths = list(filter(lambda x: x not in self.fastFive, path_locations))

        if len(missing_paths) > 0:
            self.logError("[NanoporeRead:Initialize] ERROR keys are missing from {}: {}"
                          .format(self.filename, missing_paths))
            self.close()
            return False

        self.run_id = self.bytes_to_string(self.fastFive[METADATA_PATH_KEY].attrs.get('run_id'))
        reads = list(self.fastFive[READS_KEY])
        if len(reads) != 1:
            self.logError("[NanoporeRead:Initialize] ERROR expected 1 read in {} ({}): got {}"
                          .format(self.filename, READS_KEY, reads))
            self.close()
            return False
        read_loc = os.path.join(READS_KEY, reads[0])

        self.read_label = self.bytes_to_string(self.fastFive[read_loc].attrs.get('read_id'))

        return self.run_id is not None and self.read_label is not None

    def _initialize(self):
        """Routine setup 1D NanoporeReads, returns false if basecalled with upsupported
        version or is not base-called
        """
        if not self.open():
            return False

        # parameter clarity
        perform_kmer_event_aln_if_required = self.perform_kmer_event_alignment is None
        perform_kmer_event_aln_always = self.perform_kmer_event_alignment == True

        # are we required to perform kmer event realignment?
        if perform_kmer_event_aln_always:
            oned_root_address = self.generate_new_event_table()
            if not oned_root_address:
                self.logError("[NanoporeRead:_initialize] required kmer event alignment failed for {}".format(
                    self.filename))
                self.close()
                return False
        else:
            # try to find current basecall address
            oned_root_address = self.get_latest_basecall_edition(
                self.event_table if self.event_table else TEMPLATE_BASECALL_KEY)
            # if not found, generate
            if not oned_root_address and perform_kmer_event_aln_if_required:
                oned_root_address = self.generate_new_event_table()

        # Some RNA reads have incorrectly formatted basecall tables so we check and then force load_from_raw
        if oned_root_address and self.rna and not self.has_valid_event_table_format(oned_root_address):
            self.logError("[NanoporeRead:_initialize] WARN invalid event table format for RNA read")
            if perform_kmer_event_aln_if_required:
                oned_root_address = self.generate_new_event_table()
        # sanity check
        if not oned_root_address:
            self.logError("[NanoporeRead:_initialize] ERROR could not find 1D root address in {}"
                          .format(self.filename))
            self.close()
            return False

        print("[NanoporeRead._initialize] oned_root_address {}".format(oned_root_address))
        # get basecall version
        # if not any(x in self.fastFive[oned_root_address].attrs.keys() for x in VERSION_KEY):
        #     self.logError("[NanoporeRead:_initialize] ERROR %s missing version" % self.filename)
        #     self.close()
        #     return False
        # if "version" in self.fastFive[oned_root_address].attrs.keys():
        #     self.version = self.bytes_to_string(self.fastFive[oned_root_address].attrs["version"])
        # elif "dragonet version" in self.fastFive[oned_root_address].attrs.keys():
        #     self.version = self.bytes_to_string(self.fastFive[oned_root_address].attrs["dragonet version"])
        # elif "nanotensor version" in self.fastFive[oned_root_address].attrs.keys():
        #     self.version = self.fastFive[oned_root_address].attrs["nanotensor version"]
        # else:
        #     self.version = self.fastFive[oned_root_address].attrs["signalAlign version"]

        # if self.version not in SUPPORTED_1D_VERSIONS:
        #     self.logError("[NanoporeRead:_initialize] ERROR %s unsupported version %s " % (self.filename, self.version))
        #     if self.enforce_supported_versions:
        #         self.close()
        #         return False
        #     else:
        #         self.logError(
        #             "[NanoporeRead:_initialize] unexpected behavior may be due to unexpected nanopore read version")

        self.template_event_table_address = os.path.join(oned_root_address, "BaseCalled_template/Events")
        self.template_model_address = os.path.join(oned_root_address, "BaseCalled_template/Model")
        self.template_model_id = None

        self.fastq_sequence_address = os.path.join(oned_root_address, "BaseCalled_template/Fastq")
        if self.fastq_sequence_address not in self.fastFive:
            self.logError("[NanoporeRead:_initialize]ERROR %s missing fastq" % self.filename)
            self.close()
            return False

        self.fastq = self.bytes_to_string(self.fastFive[self.fastq_sequence_address][()])
        if self.filter_reads is not None:
            record = SeqIO.read(StringIO(self.fastq), "fastq")
            average_quality = np.mean(record.letter_annotations["phred_quality"])
            if average_quality < self.filter_reads:
                self.logError("[NanoporeRead:_initialize]ERROR {} Fastq quality is below minimum threshold "
                              "{} < 7".format(self.filename, average_quality))
                self.close()
                return False

        self.template_read = self.fastq.split('\n')[1]
        if self.rna:
            # reverse and replace "U"
            self.template_read = self.template_read.replace("U", "T")[::-1]

        self.kmer_length = -1 if len(self.fastFive[self.template_event_table_address]) == 0 else \
            len(self.bytes_to_string(self.fastFive[self.template_event_table_address][0]['model_state']))
        self.template_read_length = len(self.template_read)
        if self.template_read_length <= 0 or self.kmer_length <= 0:
            self.logError("[NanoporeRead:_initialize] ERROR illegal read parameters "
                          "template_read_length: %s, kmer_length: %s"
                          % (self.template_read_length, self.kmer_length))
            self.close()
            return False

        return True

    def generate_new_event_table(self):
        if self.aligned_segment is None:
            nucleotide_sequence, nucleotide_qualities, _, _, self.aligned_segment = \
                get_full_nucleotide_read_from_alignment(self.alignment_file, self.read_label)

        oned_root_address = load_from_raw2(self, self.aligned_segment, self.model_file_location, self.path_to_bin,
                                           analysis_identifier=self.event_table if self.event_table else None,
                                           write_failed_alignments=True,
                                           rna=self.rna)
        if oned_root_address:
            self.logError(
                "[NanoporeRead:generate_new_event_table] INFO generated event table at {}".format(oned_root_address))
        else:
            self.logError("[NanoporeRead:generate_new_event_table] ERROR failed to generate event table")
        return oned_root_address

    def has_valid_event_table_format(self, oned_root_address):
        """Check if the 'start' and 'length' values are in the time scale, NOT the index scale
        :param oned_root_address: Basecalled analysis path
        :return: boolean if start and length are correct format
        """
        template_event_table_address = os.path.join(oned_root_address, "BaseCalled_template/Events")
        if template_event_table_address in self.fastFive:
            template_events = np.asarray(self.fastFive[template_event_table_address])
        check_numpy_table(template_events, req_fields=('start', 'length'))
        if template_events["start"].dtype is np.dtype('uint64'):
            return False
        else:
            return True

    @staticmethod
    def make_event_map(events, kmer_length):
        event_map = [0]
        previous_prob = 0
        for i, line in islice(enumerate(events), 1, None):
            move = line['move']
            this_prob = line['p_model_state']
            if move == 1:
                event_map.append(i)
            if move > 1:
                for skip in range(move - 1):
                    event_map.append(i - 1)
                event_map.append(i)
            if move == 0:
                if this_prob > previous_prob:
                    event_map[-1] = i
            previous_prob = this_prob
        final_event_index = [event_map[-1]]
        padding = final_event_index * (kmer_length - 1)
        event_map = event_map + padding
        return event_map

    def init_event_map(self):
        """Maps the events from the template and complement strands to their base called kmers. The map
        generated by this function is called the "strand_event_map" because it only works for mapping the
        strand read (1D read) to to it's events. Uses the same fields as 'get_twoD_event_map' below.
        """

        self.template_strand_event_map = self.make_event_map(self.template_events, self.kmer_length)
        assert len(self.template_strand_event_map) == len(self.template_read), \
            "Read and event map lengths do not match {} != {}".format(len(self.template_read),
                                                                      len(self.template_strand_event_map))

        return True

    @staticmethod
    def sequence_from_events(events):
        """Get new read from event table"""
        bases = []
        for i, event in enumerate(events):
            if i == 0:
                bases.extend([chr(x) for x in event['model_state']])

            else:
                if event['move'] > 0:
                    bases.append(bytes.decode
                                 (event['model_state'][-event['move']:]))
        sequence = ''.join(bases)
        return sequence

    def _join_path(self, *args):
        return '/'.join(args)

    def get_template_events(self):
        if self.template_event_table_address in self.fastFive:
            self.template_events = self.fastFive[self.template_event_table_address]
            return self.template_events
        return False

    def get_template_read(self, initalize_bypass=True):
        # if we have it, return it
        if self.template_read is not None and len(self.template_read) > 0:
            return self.template_read

        # if we haven't (successfully) initialized, try to get it in the same way as init
        if not self.initialize_success and initalize_bypass:
            # get configured root address
            oned_root_address = False
            if self.event_table:
                oned_root_address = self.get_latest_basecall_edition(self.event_table)
            if not oned_root_address:
                oned_root_address = self.get_latest_basecall_edition(TEMPLATE_BASECALL_KEY)
                if not oned_root_address: return False

            # try to get fastq address
            self.fastq_sequence_address = os.path.join(oned_root_address, "BaseCalled_template/Fastq")
            if self.fastq_sequence_address not in self.fastFive:
                return False

            # set and return the read (initialize will overwrite this if run)
            self.template_read = self.bytes_to_string(self.fastFive[self.fastq_sequence_address][()]).split('\n')[1]
            self.template_read_length = len(self.template_read)
            if self.rna:
                self.template_read = self.template_read.replace("U", "T")[::-1]
            return self.template_read

        return False

    def get_model_adjustments(self):
        if self.template_model_address in self.fastFive:
            self.has_template_model = True
            self.template_scale = self.fastFive[self.template_model_address].attrs["scale"]
            self.template_shift = self.fastFive[self.template_model_address].attrs["shift"]
            self.template_drift = self.fastFive[self.template_model_address].attrs["drift"]
            self.template_var = self.fastFive[self.template_model_address].attrs["var"]
            self.template_scale_sd = self.fastFive[self.template_model_address].attrs["scale_sd"]
            self.template_var_sd = self.fastFive[self.template_model_address].attrs["var_sd"]

        if self.template_model_address not in self.fastFive:
            self.has_template_model = False
        return

    def get_model_id(self, address):
        if address in self.fastFive:
            model_name = bytes.decode(self.fastFive[address].attrs["model_file"])
            model_name = model_name.split('/')[-1]
            return model_name
        else:
            return None

    def write_data(self, data, location, compression=True):
        """Write numpy data to fast5 file"""
        try:
            del self.fastFive[location]
        except KeyError:
            pass
        self.fastFive.create_dataset(location, data=data, compression=compression)

    def assert_events_and_event_map(self):
        template_events_check = self.get_template_events()
        oneD_event_map_check = self.init_event_map()

        ok = False not in [template_events_check, oneD_event_map_check]
        return ok

    def Write(self, out_file, initialize=False):
        if initialize:
            ok = self.Initialize()
            if not ok:
                self.close()
                return False

        ok = self.assert_events_and_event_map()
        if not ok:
            self.close()
            return False

        # get model params
        self.get_model_adjustments()

        # Make the npRead
        # line 1 parameters
        print(len(self.alignment_table_sequence), end=' ', file=out_file)  # 0alignment read length
        print(len(self.template_events), end=' ', file=out_file)  # 1nb of template events
        print(len(self.complement_events), end=' ', file=out_file)  # 2nb of complement events
        print(len(self.template_read), end=' ', file=out_file)  # 3length of template read
        print(len(self.complement_read), end=' ', file=out_file)  # 4length of complement read
        print(self.template_scale, end=' ', file=out_file)  # 5template scale
        print(self.template_shift, end=' ', file=out_file)  # 6template shift
        print(self.template_var, end=' ', file=out_file)  # 7template var
        print(self.template_scale_sd, end=' ', file=out_file)  # 8template scale_sd
        print(self.template_var_sd, end=' ', file=out_file)  # 9template var_sd
        print(self.template_drift, end=' ', file=out_file)  # 0template_drift
        print(self.complement_scale, end=' ', file=out_file)  # 1complement scale
        print(self.complement_shift, end=' ', file=out_file)  # 2complement shift
        print(self.complement_var, end=' ', file=out_file)  # 3complement var
        print(self.complement_scale_sd, end=' ', file=out_file)  # 4complement scale_sd
        print(self.complement_var_sd, end=' ', file=out_file)  # 5complement var_sd
        print(self.complement_drift, end=' ', file=out_file)  # 6complement_drift
        print((1 if self.twoD else 0), end='\n', file=out_file)  # has 2D

        # line 2 alignment table sequence
        print(self.alignment_table_sequence, end='\n', file=out_file)

        # line 3 template read
        print(self.template_read, end='\n', file=out_file)

        # line 4 template strand map
        for _ in self.template_strand_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 5 complement read
        print(self.complement_read, end='\n', file=out_file)

        # line 6 complement strand map
        for _ in self.complement_strand_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 7 template 2D event map
        for _ in self.template_event_map:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 8 template events
        template_start_time = self.template_events[0]['start']
        for mean, stdev, length, start in self.template_events['mean', 'stdv', 'length', 'start']:
            print(mean, stdev, length, (start - template_start_time), sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 9 complement 2D event map
        for _ in self.complement_event_map[::-1]:
            print(_, end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 10 complement events
        if self.twoD:
            complement_start_time = self.complement_events[0]['start']
            for mean, stdev, length, start in self.complement_events['mean', 'stdv', 'length', 'start']:
                print(mean, stdev, length, (start - complement_start_time), sep=' ', end=' ', file=out_file)
        else:
            pass
        print("", end="\n", file=out_file)

        # line 11 model_state (template)
        for _ in self.template_events['model_state']:
            print(bytes.decode(_), sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 12 p(model) (template)
        for _ in self.template_events['p_model_state']:
            print(_, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 13 model_state (complement)
        if self.twoD:
            for _ in self.complement_events['model_state']:
                print(bytes.decode(_), sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        # line 14 p(model) (complement)
        if self.twoD:
            for _ in self.complement_events['p_model_state']:
                print(_, sep=' ', end=' ', file=out_file)
        print("", end="\n", file=out_file)

        return True

    def is_read_rna(self):
        """
        Determine if a read is RNA or DNA
        source: https://github.com/nanoporetech/tombo/blob/master/tombo/tombo_helper.py
        """
        # check both experiment type and kit slots for "rna"
        exp_type, exp_kit = None, None
        try:
            exp_type = bytes.decode(self.fastFive['UniqueGlobalKey/context_tags'].attrs[
                                        'experiment_type'])
            # remove the word internal since it contains rna.
            exp_type = exp_type.replace('internal', '')
        except:
            pass
        try:
            exp_kit = bytes.decode(self.fastFive['UniqueGlobalKey/context_tags'].attrs[
                                       'experiment_kit'])
            # remove the word internal since it contains rna.
            exp_kit = exp_kit.replace('internal', '')
        except:
            pass

        if exp_type is None and exp_kit is None:
            rna = False
        else:
            rna = (
                    (exp_type is not None and re.search('rna', exp_type) is not None) or
                    (exp_kit is not None and re.search('rna', exp_kit) is not None))

        return rna

    @staticmethod
    def logError(message, parent_job=None):
        if parent_job is None:
            print(message, file=sys.stderr)
        else:
            parent_job.fileStore.logToMaster(message)

    @staticmethod
    def bytes_to_string(string):
        """Check string. If bytes, convert to string and return string

        :param string: string or bytes
        """
        if string is None or type(string) == str:
            return string
        elif 'bytes' in str(type(string)):
            return string.decode()
        else:
            raise AssertionError("String needs to be bytes or string ")


class NanoporeRead2D(NanoporeRead):

    def _initialize(self):
        """
        Separate initialization routine for 2D reads
        :return:
        """

        if not self.open(): return False

        self.has2D = False
        self.has2D_alignment_table = False

        if TWOD_BASECALL_KEY_0 not in self.fastFive:
            self.logError("[NanoporeRead::initialize_twoD] Didn't find twoD address, looked here {}"
                          .format(TWOD_BASECALL_KEY_0))
            self.close()
            return False

        twoD_address = self.get_latest_basecall_edition(TWOD_BASECALL_KEY)
        if twoD_address not in self.fastFive:
            self.logError("[NanoporeRead::initialize_twoD] Didn't find twoD address, looked here %s " % twoD_address)
            self.close()
            return False

        self.version = bytes.decode(self.fastFive[twoD_address].attrs["dragonet version"])

        if self.version not in SUPPORTED_2D_VERSIONS:
            self.logError("[NanoporeRead::initialize_twoD] Unsupported Version {} ({} supported)".format(
                self.version, SUPPORTED_2D_VERSIONS))
            if self.enforce_supported_versions:
                self.close()
                return False
            else:
                self.logError(
                    "[NanoporeRead:_initialize] unexpected behavior may be due to unexpected nanopore read version")

        if self.version == "1.15.0":
            oneD_address = self.get_latest_basecall_edition(TWOD_BASECALL_KEY)
        else:
            oneD_address = self.get_latest_basecall_edition(TEMPLATE_BASECALL_KEY)

        twoD_alignment_table_address = twoD_address + "/BaseCalled_2D/Alignment"
        if twoD_alignment_table_address in self.fastFive:
            self.twoD_alignment_table = self.fastFive[twoD_alignment_table_address]
            if len(self.twoD_alignment_table) > 0:
                self.has2D_alignment_table = True
            self.kmer_length = len(self.twoD_alignment_table[0][2])

        twoD_read_sequence_address = twoD_address + "/BaseCalled_2D/Fastq"
        if twoD_read_sequence_address in self.fastFive:
            self.has2D = True
            self.twoD_read_sequence = bytes.decode(self.fastFive[twoD_read_sequence_address][()].split()[2])

        # initialize version-specific paths
        if self.version == "1.15.0":
            self.template_event_table_address = twoD_address + '/BaseCalled_template/Events'
            self.template_model_address = twoD_address + "/BaseCalled_template/Model"
            self.template_model_id = self.get_model_id(twoD_address + "/Summary/basecall_1d_template")
            self.template_read = bytes.decode(self.fastFive[twoD_address + "/BaseCalled_template/Fastq"][()].split()[2])

            self.complement_event_table_address = twoD_address + '/BaseCalled_complement/Events'
            self.complement_model_address = twoD_address + "/BaseCalled_complement/Model"
            self.complement_model_id = self.get_model_id(twoD_address + "/Summary/basecall_1d_complement")
            self.complement_read = bytes.decode(
                self.fastFive[twoD_address + "/BaseCalled_complement/Fastq"][()].split()[2])
            return True

        elif self.version == "1.19.0" or self.version == "1.20.0":
            self.template_event_table_address = oneD_address + '/BaseCalled_template/Events'
            self.template_model_address = oneD_address + "/BaseCalled_template/Model"
            self.template_model_id = self.get_model_id(oneD_address + "/Summary/basecall_1d_template")
            self.template_read = bytes.decode(self.fastFive[oneD_address + "/BaseCalled_template/Fastq"][()].split()[2])

            self.complement_event_table_address = oneD_address + '/BaseCalled_complement/Events'
            self.complement_model_address = oneD_address + "/BaseCalled_complement/Model"
            self.complement_model_id = self.get_model_id(oneD_address + "/Summary/basecall_1d_complement")
            self.complement_read = bytes.decode(
                self.fastFive[oneD_address + "/BaseCalled_complement/Fastq"][()].split()[2])
            return True

        elif self.version == "1.22.2" or self.version == "1.22.4" or self.version == "1.23.0":
            self.template_event_table_address = oneD_address + '/BaseCalled_template/Events'
            self.template_model_address = ""
            self.template_model_id = None
            self.template_read = bytes.decode(self.fastFive[oneD_address + "/BaseCalled_template/Fastq"][()].split()[2])

            self.complement_event_table_address = oneD_address + '/BaseCalled_complement/Events'
            self.complement_model_address = ""
            self.complement_model_id = None
            self.complement_read = bytes.decode(
                self.fastFive[oneD_address + "/BaseCalled_complement/Fastq"][()].split()[2])
            return True
        else:
            self.logError("Unsupported Version (1.15.0, 1.19.0, 1.20.0, 1.22.2, 1.22.4 supported)")
            return False

    def assemble_2d_sequence_from_table(self):
        """The 2D read sequence contains kmers that may not map to a template or complement event, which can make
        mapping difficult downstream. This function makes a sequence from the 2D alignment table, which is usually
        pretty similar to the 2D read, except it is guaranteed to have an event map to every position.

        returns: sequence made from alignment table
        """

        def find_kmer_overlap(k_i, k_j):
            """ finds the overlap between two non-identical kmers.
            k_i: one kmer
            k_j: another kmer
            returns: The number of positions not matching
            """
            for i in range(1, len(k_i)):
                sk_i = k_i[i:]
                sk_j = k_j[:-i]
                if sk_i == sk_j:
                    return i
            return len(k_i)

        self.alignment_table_sequence = ''
        self.alignment_table_sequence = self.twoD_alignment_table[0][2]
        p_kmer = self.twoD_alignment_table[0][2]

        # iterate through the k-mers in the alignment table
        for t, c, kmer in self.twoD_alignment_table:
            # if we're at a new 6-mer
            if kmer != p_kmer:
                # find overlap, could move up to len(k-mer) - 1 bases
                i = find_kmer_overlap(p_kmer, kmer)
                # append the suffix of the new 6-mer to the sequence
                self.alignment_table_sequence += kmer[-i:]
                # update
                p_kmer = kmer
            else:
                continue
        self.alignment_table_sequence = bytes.decode(self.alignment_table_sequence)
        return

    def get_twoD_event_map(self):
        """Maps the kmers in the alignment table sequence read to events in the template and complement strand reads
        """

        def kmer_iterator(dna, k):
            for i in range(len(dna)):
                kmer = dna[i:(i + k)]
                if len(kmer) == k:
                    yield kmer

        # initialize
        alignment_row = 0
        prev_alignment_kmer = ''
        nb_template_gaps = 0
        previous_complement_event = None
        previous_template_event = None

        # twoD_init = self.initialize_twoD()
        # if twoD_init is False:
        #    return False

        if not self.has2D_alignment_table:
            print("{file} doesn't have 2D alignment table".format(file=self.filename))
            return False

        self.assemble_2d_sequence_from_table()

        # go thought the kmers in the read sequence and match up the events
        for i, seq_kmer in enumerate(kmer_iterator(self.alignment_table_sequence, self.kmer_length)):
            # assign the current row's kmer
            current_alignment_kmer = bytes.decode(self.twoD_alignment_table[alignment_row][2])

            # in the situation where there is a repeat kmer in the alignment then
            # we want to pick the best event to kmer alignment, TODO implement this
            # right now we just use the first alignment
            while current_alignment_kmer == prev_alignment_kmer:
                alignment_row += 1
                current_alignment_kmer = bytes.decode(self.twoD_alignment_table[alignment_row][2])
            # a match
            if seq_kmer == current_alignment_kmer:
                template_event = self.twoD_alignment_table[alignment_row][0]
                complement_event = self.twoD_alignment_table[alignment_row][1]

                # handle template event
                # if there is a gap, count it and don't add anything to the map
                if template_event == -1:
                    nb_template_gaps += 1

                # if there is an aligned event
                if template_event != -1:
                    # if it is an aligned event and there are no gaps, add it to the map
                    if nb_template_gaps == 0:
                        self.template_event_map.append(template_event)
                        # update
                        previous_template_event = template_event
                    # if there were gaps in the alignment we have to add 'best guess'
                    # event alignments to the map which is the current aligned event
                    if nb_template_gaps > 0:
                        self.template_event_map += [template_event] * (nb_template_gaps + 1)
                        # reset template gaps
                        nb_template_gaps = 0
                        # update
                        previous_template_event = template_event

                # handle complement event
                # if there is a gap, add the last aligned complement event to the map
                if complement_event == -1:
                    self.complement_event_map.append(previous_complement_event)

                # if there is an aligned complement event add it to the map
                if complement_event != -1:
                    self.complement_event_map.append(complement_event)
                    # update the most recent aligned complement event
                    previous_complement_event = complement_event

                # update previous alignment kmer and increment alignment row
                prev_alignment_kmer = current_alignment_kmer
                alignment_row += 1
                continue

            # not a match, meaning that this kmer in the read sequence is not
            # in the event alignment but we need to assign an event to it so
            # we use the heuristic that we use the alignment of the most
            # recent aligned events to this base
            if seq_kmer != current_alignment_kmer:
                self.template_event_map.append(previous_template_event)
                self.complement_event_map.append(previous_complement_event)
                continue

        # fill in the final events for the partial last kmer
        for _ in range(self.kmer_length - 1):
            self.template_event_map += [previous_template_event] * (nb_template_gaps + 1)
            self.complement_event_map.append(previous_complement_event)
            nb_template_gaps = 0

        # check that we have mapped all of the bases in the 2D read
        assert (len(self.template_event_map) == len(self.alignment_table_sequence))
        assert (len(self.complement_event_map) == len(self.alignment_table_sequence))
        return True

    def init_event_map(self):
        """Maps the events from the template and complement strands to their base called kmers. The map
        generated by this function is called the "strand_event_map" because it only works for mapping the
        strand read (1D read) to to it's events. Uses the same fields as 'get_twoD_event_map' below.
        """

        super(NanoporeRead2D, self).init_event_map()

        assert self.get_complement_events(), "Complement event table not present at {} in {}".format(
            self.complement_event_table_address, self.filename)

        self.complement_strand_event_map = self.make_event_map(self.complement_events, self.kmer_length)
        assert len(self.complement_strand_event_map) == len(self.complement_read), \
            "Complement read and event map lengths do not match {} != {}".format(len(self.complement_read),
                                                                                 len(self.complement_strand_event_map))
        return True

    def get_complement_events(self):
        if self.complement_event_table_address in self.fastFive:
            self.complement_events = self.fastFive[self.complement_event_table_address]
            return self.complement_events
        return False

    def get_model_adjustments(self):

        super(NanoporeRead2D, self).get_model_adjustments()

        if self.complement_model_address in self.fastFive:
            self.has_complement_model = True
            self.complement_scale = self.fastFive[self.complement_model_address].attrs["scale"]
            self.complement_shift = self.fastFive[self.complement_model_address].attrs["shift"]
            self.complement_drift = self.fastFive[self.complement_model_address].attrs["drift"]
            self.complement_var = self.fastFive[self.complement_model_address].attrs["var"]
            self.complement_scale_sd = self.fastFive[self.complement_model_address].attrs["scale_sd"]
            self.complement_var_sd = self.fastFive[self.complement_model_address].attrs["var_sd"]

        if self.complement_model_address not in self.fastFive:
            self.has_complement_model = False
        return

    def assert_events_and_event_map(self):

        oneD_check = super(NanoporeRead2D, self).assert_events_and_event_map()
        twoD_map_check = self.get_twoD_event_map()
        complement_events_check = self.get_complement_events()

        ok = False not in [oneD_check, twoD_map_check, complement_events_check]
        return ok
