#!/usr/bin/env python
"""Fast5 python handler"""
########################################################################
# File: fast5.py
#  executable: fast5.py
#
# Authors: Most of the initial work for this class was completed by the people at Nanonet
#  https://github.com/nanoporetech/nanonet.
#
# Further additions done by Andrew Bailey
#
# History: 12/19/17
########################################################################

import os
import sys
from glob import glob
import subprocess
import shutil
import re
import h5py
import numpy as np
import numpy.lib.recfunctions as nprf
from copy import deepcopy
from py3helpers.utils import check_numpy_table
from py3helpers.seq_tools import check_fastq_line


def short_names(fname):
    filename_short = os.path.splitext(os.path.basename(fname))[0]
    short_name_match = re.search(re.compile(r'ch\d+_file\d+'), filename_short)
    name_short = filename_short
    if short_name_match:
        name_short = short_name_match.group()
    return filename_short, name_short


class Fast5(h5py.File):
    """Class for grabbing data from single read fast5 files. Many attributes/
    groups are assumed to exist currently (we're concerned mainly with reading).
    Needs some development to make robust and for writing.

    """
    __base_analysis__ = '/Analyses'
    __event_detect_name__ = 'EventDetection'
    __default_event_path__ = 'Reads'
    __raw_path__ = '/Raw/Reads'
    __raw_name_old__ = 'RawData'
    __raw_path_old__ = '{}/{}/'.format(__base_analysis__, __raw_name_old__)
    __raw_signal_path_old__ = '{}/Signal'.format(__raw_path_old__)
    __raw_meta_path_old__ = '{}/Meta'.format(__raw_path_old__)
    __channel_meta_path__ = '/UniqueGlobalKey/channel_id'
    __tracking_id_path__ = 'UniqueGlobalKey/tracking_id'
    __context_tags_path__ = 'UniqueGlobalKey/context_tags'

    __default_basecall_2d_analysis__ = 'Basecall_2D'
    __default_basecall_1d_analysis__ = 'Basecall_1D'

    __default_seq_section__ = '2D'
    __default_basecall_fastq__ = 'BaseCalled_{}/Fastq'
    __default_basecall_1d_events__ = 'BaseCalled_{}/Events'
    __default_basecall_1d_model__ = 'BaseCalled_{}/Model'
    __default_basecall_1d_summary__ = 'Summary/basecall_1d_{}'

    __default_alignment_analysis__ = 'Alignment'

    __default_hairpin_split_analysis__ = 'Hairpin_Split'
    __template_section__ = 'template'
    __complement_section__ = 'complement'
    __default_section__ = __template_section__

    __default_mapping_analysis__ = 'Squiggle_Map'
    __default_mapping_events__ = 'SquiggleMapped_{}/Events'
    __default_mapping_model__ = 'SquiggleMapped_{}/Model'
    __default_mapping_summary__ = 'Summary/squiggle_map_{}'

    __default_substep_mapping_analysis__ = 'Substate_Map'
    __default_substep_mapping_events__ = '/Events'

    __default_basecall_mapping_analysis__ = 'AlignToRef'
    __default_basecall_mapping_events__ = 'CurrentSpaceMapped_{}/Events/'
    __default_basecall_mapping_summary__ = '/Summary/current_space_map_{}/'  # under AlignToRef analysis
    __default_basecall_alignment_summary__ = '/Summary/genome_mapping_{}/'  # under Alignment analysis

    #todo fix the form of these
    __default_corrected_genome__ = '/Analyses/RawGenomeCorrected_000/BaseCalled_template'  # nanoraw
    __default_signalalign_events__ = '/Analyses/SignalAlign_00{}'  # signalalign events
    __default_resegment_basecall__ = '/Analyses/ReSegmentBasecall_00{}'
    __default_eventalign_events__ = '/Analyses/EventAlign_00{}'

    __default_event_table_fields__ = ('start', 'length', 'mean', 'stdv')

    __default_engine_state_path__ = '/EngineStates/'
    __temp_fields__ = ('heatsink', 'asic')

    def __init__(self, fname, read='r'):
        super(Fast5, self).__init__(fname, read)

        # Attach channel_meta as attributes, slightly redundant
        for k, v in self[self.__channel_meta_path__].attrs.items():
            setattr(self, k, v)
        # Backward compat.
        self.sample_rate = self.sampling_rate

        self.filename_short, self.name_short = short_names(self.filename)

    @classmethod
    def New(cls, fname, read='a', tracking_id={}, context_tags={}, channel_id={}):
        """Construct a fresh single-read file, with meta data written to
        standard locations.

        """
        # TODO: tracking_id and channel_id checks, do we care for these?
        #       simply copies the data "verbatim" (it doesn't copy the group
        #       directly, but reconstructs it in a parallel layout). channel_id
        #       is enough for most purposes I (cjw) think, if we take filenames
        #       to be a key to filter by.
        req_fields = ['channel_number', 'offset', 'range', 'digitisation', 'sampling_rate']
        if not set(req_fields).issubset(set(channel_id.keys())):
            raise KeyError(
                'channel_id does not contain required fields: {},\ngot {}.'.format(req_fields, channel_id.keys())
            )

        # Start a new file, populate it with meta
        with h5py.File(fname, 'w') as h:
            h.attrs['file_version'] = 1.0
            for data, location in zip(
                    [tracking_id, context_tags],
                    [cls.__tracking_id_path__, cls.__context_tags_path__]
            ):
                # cjw: no idea why these must be str, just following ossetra
                cls.__add_attrs(h, data, location, convert=str)
            # These aren't forced to be str
            cls.__add_attrs(h, channel_id, cls.__channel_meta_path__)

        # return instance from new file
        return cls(fname, read, update=False)

    def _add_attrs(self, data, location, convert=None):
        """Convenience method for adding attrs to a possibly new group.
        :param data: dict of attrs to add
        :param location: hdf path
        :param convert: function to apply to all dictionary values
        """
        self.__add_attrs(self, data, location, convert=convert)

    @staticmethod
    def __add_attrs(self, data, location, convert=None):
        """Implementation of _add_attrs as staticmethod. This allows
        functionality to be used in .New() constructor but is otherwise nasty!
        """
        if location not in self:
            self.create_group(location)
        attrs = self[location].attrs
        for k, v in data.items():
            if convert is not None:
                attrs[k] = convert(v)
            else:
                attrs[k] = v

    def _add_string_dataset(self, data, location):
        assert type(data) == str, 'Need to supply a string'
        self.create_dataset(location, data=data)

    def _add_numpy_table(self, data, location):
        self.create_dataset(location, data=data, compression=True)

    def _add_event_table(self, data, location):
        if not isinstance(data, np.ndarray):
            raise TypeError('Table is not a ndarray.')

        req_fields = ['mean', 'stdv', 'start', 'length']
        if not set(req_fields).issubset(data.dtype.names):
            raise KeyError(
                'Array does not contain fields for event array: {}, got {}.'.format(
                    req_fields, data.dtype.names
                )
            )
        self._add_numpy_table(data, location)

    def _join_path(self, *args):
        return '/'.join(args)

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


    @property
    def writable(self):
        """Can we write to the file."""
        if self.mode is 'r':
            return False
        else:
            return True

    def assert_writable(self):
        assert self.writable, "File not writable, opened with {}.".format(self.mode)

    @property
    def channel_meta(self):
        """Channel meta information as python dict"""
        return dict(self[self.__channel_meta_path__].attrs)

    @property
    def tracking_id(self):
        """Tracking id meta information as python dict"""
        return dict(self[self.__tracking_id_path__].attrs)

    @property
    def raw_attributes(self):
        """Attributes for a read, assumes one read in file"""
        return dict(self.get_read(group=True, raw=True).attrs)

    @property
    def event_attributes(self):
        """Attributes for a read, assumes one read in file"""
        return dict(self.get_read(group=True).attrs)

    def summary(self, rename=True, delete=True, scale=True):
        """A read summary, assumes one read in file"""
        to_rename = zip(
            ('start_mux', 'abasic_found', 'duration', 'median_before'),
            ('mux', 'abasic', 'strand_duration', 'pore_before')
        )
        to_delete = ('read_number', 'scaling_used')

        data = deepcopy(self.event_attributes)
        data['filename'] = os.path.basename(self.filename)
        data['run_id'] = self.tracking_id['run_id']
        data['channel'] = self.channel_meta['channel_number']
        if scale:
            data['duration'] /= self.channel_meta['sampling_rate']
            data['start_time'] /= self.channel_meta['sampling_rate']

        if rename:
            for i,j in to_rename:
                try:
                    data[j] = data[i]
                    del data[i]
                except KeyError:
                    pass
        if delete:
            for i in to_delete:
                try:
                    del data[i]
                except KeyError:
                    pass

        for key in data:
            if isinstance(data[key], float):
                data[key] = np.round(data[key], 4)

        return data

    def strip_analyses(self, keep=('{}_000'.format(__event_detect_name__), __raw_path__)):
        """Remove all analyses from file

        :param keep: whitelist of analysis groups to keep

        """
        analyses = self[self.__base_analysis__]
        for name in analyses.keys():
            if name not in keep:
                del analyses[name]

    def repack(self):
        """Run h5repack on the current file. Returns a fresh object."""
        path = os.path.abspath(self.filename)
        path_tmp = '{}.tmp'.format(path)
        mode = self.mode
        self.close()
        subprocess.check_call(['h5repack', path, path_tmp], stderr=subprocess.STDOUT, shell=False)
        shutil.move(path_tmp, path)
        return Fast5(path, mode)

    def create_copy(self, copy_path):
        """Run h5repack on the current file. Returns a fresh object."""
        path = os.path.abspath(self.filename)
        path_tmp = copy_path
        mode = self.mode
        self.close()
        subprocess.check_call(['h5repack', path, path_tmp], stderr=subprocess.STDOUT, shell=False)
        return Fast5(path_tmp, mode)

    ###
    # Extracting read event data

    def get_reads(self, group=False, raw=False, read_numbers=None, scale=True):
        """Iterator across event data for all reads in file

        :param group: return hdf group rather than event data
        """
        if not raw:
            event_group = self.get_analysis_latest(self.__event_detect_name__)
            event_path = self._join_path(event_group, self.__default_event_path__)
            reads = self[event_path]
        else:
            try:
                reads = self[self.__raw_path__]
            except:
                raise KeyError('No raw data available in file {}.'.format(self.filename))

        if read_numbers is None:
            it = reads.keys()
        else:
            it = (k for k in reads.keys()
                  if reads[k].attrs['read_number'] in read_numbers)

        if group == 'all':
            for read in it:
                yield reads[read], read
        elif group:
            for read in it:
                yield reads[read]
        else:
            for read in it:
                if not raw:
                    yield self._get_read_data(reads[read])
                else:
                    yield self._get_read_data_raw(reads[read], scale=scale)

    def get_read(self, group=False, raw=False, read_number=None, scale=False):
        """Like get_reads, but only the first read in the file

        :param group: return hdf group rather than event/raw data
        """
        if read_number is None:
            return self.get_reads(group, raw, scale=scale).__next__()
        else:
            return self.get_reads(group, raw, read_numbers=[read_number], scale=scale).__next__()

    def get_corrected_events(self):
        """Returns corrected events table along with the start relative to raw data"""
        try:
            reads = self[self.__default_corrected_genome__]
            events = reads['Events']
            attributes = dict(events.attrs.items())
            corr_start_rel_to_raw = attributes['read_start_rel_to_raw']
        except KeyError:
            raise KeyError('Read does not contain required fields: {}'.format(self.__default_corrected_genome__))
        return np.asarray(events), corr_start_rel_to_raw

    #todo fix path creation
    def get_signalalign_events(self, mea=False, sam=False):
        """Get signal align events, sam or mea alignment"""
        assert (not mea or not sam), "Both mea and sam cannot be set to True"
        try:
            path = self.check_path(self.__default_signalalign_events__, latest=True)
            reads = self[path]
            if mea:
                events = np.asarray(reads['MEA_alignment_labels'])
            elif sam:
                events = str(np.asarray(reads['sam']))
            else:
                events = np.asarray(reads['full'])

        except KeyError:
            raise KeyError('Read does not contain required fields: {}'.format(path))
        return events

    #todo fix path creation
    def get_eventalign_events(self, section=__default_section__):
        """Get signal align events, sam or mea alignment"""
        assert section in [self.__template_section__, self.__complement_section__], \
            "Section must be template or complement: {}".format(section)
        try:
            path = self.check_path(self.__default_eventalign_events__, latest=True)
            reads = self[path]
            events = np.asarray(reads['BaseCalled_{}/Events'.format(section)])

        except KeyError:
            raise KeyError('Read does not contain required fields: {}'.format(path))
        return events

    #todo fix path creation
    def get_resegment_basecall(self, number=None):
        """Get most recent resegmented basecall events table

        :param number: manually select the resegment version

        """
        try:
            if number:
                path = self.__default_resegment_basecall__.format(number)
            else:
                path = self.check_path(self.__default_resegment_basecall__, latest=True)
            reads = self[path]
            events = reads['BaseCalled_template/Events']
        except KeyError:
            raise KeyError('Read does not contain required fields: {}'.format(path))
        return np.asarray(events)


    def _get_read_data(self, read, indices=None):
        """Private accessor to read event data"""
        # We choose the following to always be floats
        float_fields = ('start', 'length', 'mean', 'stdv')

        events = read['Events']

        # We assume that if start is an int or uint the data is in samples
        #    else it is in seconds already.
        needs_scaling = False
        if events['start'].dtype.kind in ['i', 'u']:
            needs_scaling = True

        dtype = np.dtype([(d[0], 'float') if d[0] in float_fields else d
                          for d in events.dtype.descr
                          ])
        data = None
        with events.astype(dtype):
            if indices is None:
                data = events[()]
            else:
                try:
                    data = events[indices[0]:indices[1]]
                except:
                    raise ValueError(
                        'Cannot retrieve events using {} as indices'.format(indices)
                    )

        # File spec mentions a read.attrs['scaling_used'] attribute,
        #    its not clear what this is. We'll ignore it and hope for
        #    the best.
        if needs_scaling:
            data['start'] /= self.sample_rate
            data['length'] /= self.sample_rate
        return data

    def _get_read_data_raw(self, read, indices=None, scale=True):
        """Private accessor to read raw data"""
        raw = read['Signal']
        dtype = float if scale else int

        data = None
        with raw.astype(dtype):
            if indices is None:
                data = raw[()]
            else:
                try:
                    data = raw[indices[0]:indices[1]]
                except:
                    raise ValueError(
                        'Cannot retrieve events using {} as indices'.format(indices)
                    )

        # Scale data to pA
        if scale:
            meta = self.channel_meta
            raw_unit = meta['range'] / meta['digitisation']
            data = (data + meta['offset']) * raw_unit
        return data

    #todo fix path creation
    def set_read(self, data, meta, scale=True):
        """Write event data to file

        :param data: event data
        :param meta: meta data to attach to read
        :param read_number: per-channel read counter
        """
        req_fields = [
            'start_time', 'duration', 'read_number',
            'start_mux', 'read_id', 'scaling_used'
        ]
        if not set(req_fields).issubset(meta.keys()):
            raise KeyError(
                'Read meta does not contain required fields: {}, got {}'.format(
                    req_fields, meta.keys()
                )
            )
        self.test_event_table(data)
        path = self._join_path(
            self.__event_path__, 'Read_{}'.format(meta['read_number'])
        )
        self._add_attrs(meta, path)

        # (see _get_read_data()). If the data is not an int or uint
        #    we assume it is in seconds and scale appropriately
        if scale:
            data['start'] *= self.sample_rate
            data['length'] *= self.sample_rate

        self._add_event_table(data, self._join_path(path, 'Events'))

    def set_fastq(self, destination_root, data, section=__default_section__, overwrite=False):
        """Write new fasta file to file

        :param destination_root: root directory; data will be stored in {destination_root}/Basecalled_{section}/Fastq
        :param data: fastq file
        :param section: name of basecall analysis default (template)
        """
        check_fastq_line(data)

        # get location and sanity check
        path = self._join_path(destination_root, self.__default_basecall_fastq__.format(section))
        if path in self:
            if overwrite:
                self.delete(path, ignore=True)
            else:
                raise Exception("Destination {} already exists in {}".format(path, self.filename))

        # save
        self._add_string_dataset(data, path)

    def set_event_table(self, destination_root, data, meta, section=__default_section__, scale=False, overwrite=False):
        """Write new event data to file

        :param destination_root: root directory; data will be stored in {destination_root}/Basecalled_{section}/Events
        :param data: event data
        :param meta: meta data to attach to read
        :param section: name of basecall analysis default (template)
        :param scale: scale the the start and length by the sample rate
        :param overwrite: overwrite most recent path
        """

        self.assert_writable()
        self.test_event_table(data)

        # modification to data
        if meta:
            #todo add attrs to dest_root or dest_events?
            self._add_attrs(meta, destination_root)
        if scale:
            data['start'] *= self.sample_rate
            data['length'] *= self.sample_rate

        # get location and sanity check
        destination_events = self._join_path(destination_root, self.__default_basecall_1d_events__.format(section))
        if destination_events in self:
            if overwrite:
                self.delete(destination_events, ignore=True)
            else:
                raise Exception("Destination {} already exists in {}".format(destination_events, self.filename))

        # save
        self._add_event_table(data, destination_events)

    #todo fix path creation
    def set_eventalign_table(self, template, complement, meta, overwrite=False):
        """Write eventalign table to fast5 file

        :param template: template dataset
        :param complement: complement dataset
        :param meta: meta data to attach to read
        :param overwrite: overwrite most recent path
        """
        assert template or complement, "Must set template and/or complement dataset"
        self.assert_writable()
        path = "EventAlign_00{}"
        path = self._join_path(self.__base_analysis__, path)
        path = self.check_path(path, latest=overwrite)
        if overwrite:
            self.delete(path, ignore=True)
        if meta:
            self._add_attrs(meta, path)
        if not isinstance(template, np.ndarray):
            raise TypeError('Table is not a ndarray.')
        if not isinstance(complement, np.ndarray):
            raise TypeError('Table is not a ndarray.')

        if template is not None:
            self._add_numpy_table(template, self._join_path(path, "BaseCalled_{}".format("template"), 'Events'))
        if complement is not None:
            self._add_numpy_table(complement, self._join_path(path, "BaseCalled_{}".format("complement"), 'Events'))

        return True

    #todo change to look like get_analysis_latest
    def check_path(self, path, latest=False):
        """Check if path exists, if it does increment numbering

        :param path: path to fast5 object. Needs to have a field where string.format can work! """
        highest = 0
        while highest < 20:
            if path.format(highest) in self:
                highest += 1
                continue
            else:
                if latest and highest > 0:
                    return path.format(highest-1)  # the last base-called version we saw
                else:
                    return path.format(highest)  # the new base-called version

    def get_read_stats(self):
        """Combines stats based on events with output of .summary, assumes a
        one read file.

        """
        data = deepcopy(self.summary())
        read = self.get_read()
        sorted_means = np.sort(read['mean'])
        n_events = len(sorted_means)
        n10 = int(0.1*n_events)
        n90 = int(0.9*n_events)
        data['range_current'] = sorted_means[n90] - sorted_means[n10]
        data['median_current'] = sorted_means[int(0.5*n_events)] # could be better
        data['num_events'] = n_events
        data['median_sd'] = np.median(read['stdv'])
        return data

    @staticmethod
    def test_event_table(data, req_fields=__default_event_table_fields__):
        """Wrapper function to test if event tables have required fields
        :param data: numpy array
        :param req_fields: required fields for event table """
        return check_numpy_table(data, req_fields)

    ###
    # Raw Data

    def set_raw(self, raw, meta=None, read_number=None):
        """Set the raw data in file.

        :param raw: raw data to add
        :param read_number: read number (as usually given in filename and
            contained within HDF paths, viz. Reads/Read_<>/). If not given
            attempts will be made to guess the number (assumes single read
            per file).
        """
        # Attempt to guess read_number
        if read_number is None:
            if sum(1 for _ in self.get_reads()) == 1:
                read_number = self.get_read(group=True).attrs['read_number']
            else:
                raise RuntimeError("'read_number' not given and cannot guess.")

        # Attempt to guess meta
        if meta is None:
            try:
                meta = dict(self.get_read(group=True, read_number=read_number).attrs)
            except KeyError:
                raise RuntimeError("'meta' not given and cannot guess.")

        # Clean up keys as per spec. Note: TANG-281 found that 'read_id' is not
        #   always present on the read event data, such that if we have copied
        #   meta from there, we won't have 'read_id'. The following:
        #      https://wiki/display/OFAN/Single-read+fast5+file+format
        #   notes that prior to MinKNOW version 49.2 'read_id' was not present.
        #   Why do we have a specification?
        req_keys = ['start_time', 'duration', 'read_number', 'start_mux'] #'read_id'
        meta = {k:v for k,v in meta.iteritems() if k in req_keys}
        if len(meta.keys()) != len(req_keys):
            raise KeyError(
                'Raw meta data must contain keys: {}.'.format(req_keys)
            )
        # Check meta is same as that for event data, if any
        try:
            event_meta = dict(self.get_read(group=True, read_number=read_number).attrs)
        except:
            pass
        else:
            if sum(meta[k] != event_meta[k] for k in meta.keys()) > 0:
                raise ValueError(
                    "Attempted to set raw meta data as {} "
                    "but event meta is {}".format(meta, event_meta)
                )

        # Good to go!
        read_path = self._join_path(self.__raw_path__, 'Read_{}'.format(read_number))
        data_path = self._join_path(read_path, 'Signal')
        self._add_attrs(meta, read_path)
        self[data_path] = raw

        ###
    # Analysis path resolution

    def get_analysis_latest(self, name):
        """Get group of latest (present) analysis with a given base path.

        :param name: Get the (full) path of newest analysis with a given base
            name.
        """
        try:
            return self._join_path(
                self.__base_analysis__,
                sorted(filter(
                    lambda x: name in x, self[self.__base_analysis__].keys()
                ))[-1]
            )
        except (IndexError, KeyError):
            raise IndexError('No analyses with name {} present.'.format(name))

    def get_analysis_new(self, name):
        """Get group path for new analysis with a given base name.

        :param name: desired analysis name
        """

        # Formatted as 'base/name_000'
        try:
            latest = self.get_analysis_latest(name)
            root, counter = latest.rsplit('_', 1)
            counter = int(counter) + 1
        except IndexError:
            # Nothing present
            root = self._join_path(
                self.__base_analysis__, name
            )
            counter = 0
        return '{}_{:03d}'.format(root, counter)

    # The remaining are methods to read and write data as chimaera produces
    #    It is necessarily all a bit nasty, but should provide a more
    #    consistent interface to the files. Paths are defaulted

    ###
    # Temperature etc.

    def get_engine_state(self, state, time=None):
        """Retrieve engine state from {}, either across the whole read
        (default) or at a given time.

        :param state: name of engine state
        :param time: time (in seconds) at which to retrieve temperature

        """
        location = self._join_path(
            self.__default_engine_state_path__, state
        )
        states = self[location][()]
        if time is None:
            return states
        else:
            i = np.searchsorted(states['time'], time) - 1
            return states[state][i]

    def get_temperature(self, time=None, field=__temp_fields__[0]):
        """Retrieve temperature data from {}, either across the whole read
        (default) or at a given time.

        :param time: time at which to get temperature
        :param field: one of {}

        """
        if field not in self.__temp_fields__:
            raise RuntimeError("'field' argument must be one of {}.".format(self.__temp_fields__))

        return self.get_engine_state('minion_{}_temperature'.format(field), time)

    def set_engine_state(self, data):
        """Set the engine state data.

        :param data: a 1D-array containing two fields, the first of which
            must be named 'time'. The name of the second field will be used
            to name the engine state and be used in the dataset path.
        """
        fields = data.dtype.names
        if fields[0] != 'time':
            raise ValueError("First field of engine state data must be 'time'.")
        if len(fields) != 2:
            raise ValueError("Engine state data must contain exactly two fields.")

        state = fields[1]
        location = self._join_path(
            self.__default_engine_state_path__, state
        )
        self[location] = data


    ###
    # Template/complement splitting data
    __split_summary_location__ = '/Summary/split_hairpin'

    def set_split_data(self, data, analysis=__default_hairpin_split_analysis__):
        """Write a dict containing split point data.

        :param data: `dict`-like object containing attrs to add
        :param analysis: Base analysis name (under {})

        .. warning::
            Not checking currently for required fields.
        """

        location = self._join_path(
            self.get_analysis_new(analysis), self.__split_summary_location__
        )
        self._add_attrs(data, location)

    def get_split_data(self, analysis=__default_hairpin_split_analysis__):
        """Get template-complement segmentation data.

        :param analysis: Base analysis name (under {})
        """

        location = self._join_path(
            self.get_analysis_latest(analysis), self.__split_summary_location__
        )
        try:
            return dict(self[location].attrs)
        except:
            raise ValueError(
                'Could not retrieve template-complement split point data from attributes of {}'.format(location)
            )

    def get_section_indices(self, analysis=__default_hairpin_split_analysis__):
        """Get two tuples indicating the event indices for the template and
        complement boundaries.

        :param analysis: Base analysis path (under {})
        """

        # TODO: if the below fails, calculating the values on the fly would be
        #       a useful feature. Which brings to mind could we do such a thing
        #       in all cases of missing data? Probably not reasonble.
        attrs = self.get_split_data(analysis)
        try:
            return (
                (attrs['start_index_temp'], attrs['end_index_temp']),
                (attrs['start_index_comp'], attrs['end_index_comp'])
            )
        except:
            raise ValueError('Could not retrieve template-complement segmentation data.')

    def get_section_events(self, section, analysis=__default_hairpin_split_analysis__):
        """Get the template event data.

        :param analysis: Base analysis path (under {})
        """

        indices = self.get_section_indices(analysis)
        read = self.get_read(group=True)
        events = None
        if section == 'template':
            events = self._get_read_data(read, indices[0])
        elif section == 'complement':
            events = self._get_read_data(read, indices[1])
        else:
            raise ValueError(
                '"section" parameter for fetching events must be "template" or "complement".'
            )
        return events

    ###
    # 1D Basecalling data


    def has_basecall_data(self, section=__default_section__, analysis=__default_basecall_1d_analysis__):
        """
        Determines whether events table exists

        :param section: String to use in paths, e.g. 'template' or 'complement'.
        :param analysis: Base analysis name (under {})
        """

        try:
            self.get_basecall_data(section=section, analysis=analysis)
            return True
        except:
            return False


    def get_basecall_data(self, section=__default_section__, analysis=__default_basecall_1d_analysis__):
        """Read the annotated basecall_1D events from the fast5 file.

        :param section: String to use in paths, e.g. 'template' or 'complement'.
        :param analysis: Base analysis name (under {})
        """

        base = self.get_analysis_latest(analysis)
        events_path = self._join_path(base, self.__default_basecall_1d_events__.format(section))
        try:
            return self[events_path][()]
        except:
            raise ValueError('Could not retrieve basecall_1D data from {}'.format(events_path))


    def get_alignment_attrs(self, section=__default_section__, analysis=__default_alignment_analysis__):
        """Read the annotated alignment meta data from the fast5 file.

        :param section: String to use in paths, e.g. 'template' or 'complement'.
        :param analysis: Base analysis name (under {})

        """

        attrs = None
        base = self.get_analysis_latest(analysis)
        attr_path = self._join_path(base,
                                    self.__default_basecall_alignment_summary__.format(section))
        try:
            attrs = dict(self[attr_path].attrs)
        except:
            raise ValueError('Could not retrieve alignment attributes from {}'.format(attr_path))

        return attrs

    ###
    # Mapping data

    def get_mapping_data(self, section=__default_section__, analysis=__default_mapping_analysis__):
        """Read the annotated mapping events from the fast5 file.

        .. note::
            The seq_pos column for the events table returned from basecall_mapping is
            adjusted to be the genome position (consistent with squiggle_mapping)

        :param section: String to use in paths, e.g. 'template' or 'complement'.
        :param analysis: Base analysis name (under {}). For basecall mapping
            use analysis = 'AlignToRef'.
        """

        events = None
        if analysis == self.__default_mapping_analysis__:
            # squiggle_mapping
            base = self.get_analysis_latest(analysis)
            event_path = self._join_path(base, self.__default_mapping_events__.format(section))
            try:
                events = self[event_path][()]
            except:
                raise ValueError('Could not retrieve squiggle_mapping data from {}'.format(event_path))
            attrs = self.get_mapping_attrs(section=section)

        elif analysis == self.__default_substep_mapping_analysis__:
            # substep mapping
            base = self.get_analysis_latest(analysis)
            event_path = self._join_path(base, self.__default_substep_mapping_events__.format(section))
            try:
                events = self[event_path][()]
            except:
                raise ValueError('Could not retrieve substep_mapping data from {}'.format(event_path))
            attrs=None

        else:
            # basecall_mapping
            base = self.get_analysis_latest(analysis)
            event_path = self._join_path(base, self.__default_basecall_mapping_events__.format(section))
            try:
                events = self[event_path][()]
            except:
                raise ValueError('Could not retrieve basecall_mapping data from {}'.format(event_path))

            # Modify seq_pos to be the actual genome position (consistent with squiggle_map)
            attrs = self.get_mapping_attrs(section=section, analysis=self.__default_alignment_analysis__)
            if attrs['direction'] == '+':
                events['seq_pos'] = events['seq_pos'] + attrs['ref_start']
            else:
                events['seq_pos'] = attrs['ref_stop'] - events['seq_pos']

        # add transition field
        if attrs:
            move = np.ediff1d(events['seq_pos'], to_begin=0)
            if attrs['direction'] == '-':
                move *= -1
            events = nprf.append_fields(events, 'move', move)

        return events


    def get_mapping_attrs(self, section=__default_section__, analysis=__default_mapping_analysis__):
        """Read the annotated mapping meta data from the fast5 file.
        Names which are inconsistent between squiggle_mapping and basecall_mapping are added to
        basecall_mapping (thus duplicating the attributes in basecall mapping).

        :param section: String to use in paths, e.g. 'template' or 'complement'.
        :param analysis: Base analysis name (under {})
                         For basecall mapping use analysis = 'Alignment'
        """

        attrs = None
        if analysis == self.__default_mapping_analysis__:
            # squiggle_mapping
            base = self.get_analysis_latest(analysis)
            attr_path = self._join_path(base, self.__default_mapping_summary__.format(section))
            try:
                attrs = dict(self[attr_path].attrs)
            except:
                raise ValueError('Could not retrieve squiggle_mapping meta data from {}'.format(attr_path))
        else:
            # basecall_mapping

            # AligToRef attributes (set AlignToRef first so that Alignment attrs are not overwritten)
            base = self.get_analysis_latest(self.__default_basecall_mapping_analysis__)
            attr_path = self._join_path(base, self.__default_basecall_mapping_summary__.format(section))
            try:
                attrs = dict(self[attr_path].attrs)
            except:
                raise ValueError('Could not retrieve basecall_mapping meta data from {}'.format(attr_path))

            # Rename some of the fields
            rename = [
                ('genome_start', 'ref_start'),
                ('genome_end', 'ref_stop'),
            ]
            for old, new in rename:
                attrs[new] = attrs.pop(old)

            # Alignment attributes
            base = self.get_analysis_latest(analysis)
            attr_path = self._join_path(
                base, self.__default_basecall_alignment_summary__.format(section))
            try:
                genome = self[attr_path].attrs.get('genome')
            except:
                raise ValueError('Could not retrieve basecall_mapping genome field from {}'.format(attr_path))
            try:
                attrs['reference'] = (self.get_reference_fasta(section = section)).split('\n')[1]
            except:
                raise ValueError('Could not retrieve basecall_mapping fasta from Alignment analysis')

            # Add attributes with keys consistent with Squiggle_map
            rc = '_rc'
            is_rc = genome.endswith(rc)
            attrs['ref_name'] = genome[:-len(rc)] if is_rc else genome
            attrs['direction'] = '-' if is_rc else '+'

        # Trim any other fields, the allowed are those produced by
        #   squiggle_mapping. We allow strand_score but do not require
        #   it since our writer does not require it.
        required = [
            'direction', 'ref_start', 'ref_stop', 'ref_name',
            'num_skips', 'num_stays', 'reference'
        ]
        additional = ['strand_score', 'shift', 'scale', 'drift', 'var', 'scale_sd', 'var_sd']
        keep = required + additional
        assert set(required).issubset(set(attrs)), 'Required mapping attributes not found'
        for key in (set(attrs) - set(keep)):
            del(attrs[key])

        return attrs


    def get_any_mapping_data(self, section=__default_section__, attrs_only=False):
        """Convenience method for extracting whatever mapping data might be
        present, favouring squiggle_mapping output over basecall_mapping.

        :param section: (Probably) one of '2D', 'template', or 'complement'
        :param attrs_only: Use attrs_only=True to return mapping attributes without events

        :returns: the tuple (events, attrs) or attrs only
        """
        events = None
        attrs = None
        try:
            if not attrs_only:
                events = self.get_mapping_data(section=section)
            attrs = self.get_mapping_attrs(section=section)
        except Exception as e:
            try:
                if not attrs_only:
                    events = self.get_mapping_data(section=section,
                                                   analysis=self.__default_basecall_mapping_analysis__)
                attrs = self.get_mapping_attrs(section=section,
                                               analysis=self.__default_alignment_analysis__)
            except Exception as e:
                raise ValueError(
                    "Cannot find any mapping data at paths I know about in {}. "
                    "Consider using get_mapping_data() with analysis argument."
                        .format(self.filename)
                )
        if not attrs_only:
            return events, attrs
        else:
            return attrs

    ###
    # Sequence data

    def get_fastq(self, analysis=__default_basecall_2d_analysis__, section=__default_seq_section__, custom=None):
        """Get the fastq (sequence) data.

        :param analysis: Base analysis name (under {})
        :param section: (Probably) one of '2D', 'template', or 'complement'
        :param custom: Custom hdf path overriding all of the above.
        """

        err_msg = 'Could not retrieve sequence data from {}'

        if custom is not None:
            location = custom
        else:
            location = self._join_path(
                self.get_analysis_latest(analysis), self.__default_basecall_fastq__.format(section)
            )
        try:
            return self.bytes_to_string(self[location][()])
        except:
            # Did we get given section != 2D and no analysis, that's
            #    more than likely incorrect. Try alternative analysis
            if section != self.__default_seq_section__ and analysis == self.__default_basecall_2d_analysis__:
                location = self._join_path(
                    self.get_analysis_latest(self.__default_basecall_1d_analysis__),
                    self.__default_basecall_fastq__.format(section)
                )
                try:
                    return self.bytes_to_string(self[location][()])
                except:
                    raise ValueError(err_msg.format(location))
            else:
                raise ValueError(err_msg.format(location))

    def get_sam(self, analysis=__default_alignment_analysis__, section=__default_seq_section__, custom=None):
        """Get SAM (alignment) data.

        :param analysis: Base analysis name (under {})
        :param section: (Probably) one of '2D', 'template', or 'complement'
        :param custom: Custom hdf path overriding all of the above.
        """

        if custom is not None:
            location = custom
        else:
            location = self._join_path(
                self.get_analysis_latest(analysis), 'Aligned_{}'.format(section), 'SAM'
            )
        try:
            return self[location][()]
        except:
            raise ValueError('Could not retrieve SAM data from {}'.format(location))

    def get_reference_fasta(self, analysis=__default_alignment_analysis__, section=__default_seq_section__, custom=None):
        """Get fasta sequence of known DNA fragment for the read.

        :param analysis: Base analysis name (under {})
        :param section: (Probably) one of '2D', 'template', or 'complement'
        :param custom: Custom hdf path overriding all of the above.
        """

        if custom is not None:
            location = custom
        else:
            location = self._join_path(
                self.get_analysis_latest(analysis), 'Aligned_{}'.format(section), 'Fasta'
            )
        try:
            return self[location][()]
        except:
            raise ValueError('Could not retrieve sequence data from {}'.format(location))

    def delete(self, section, ignore=False):
        """Delete a section of the H5file"""
        self.assert_writable()
        try:
            del self[section]
        except KeyError:
            if ignore:
                pass
            else:
                raise KeyError("{} not found in Fast5 file".format(section))

    def is_read_rna(self):
        """
        Determine if a read is RNA or DNA
        source: https://github.com/nanoporetech/tombo/blob/master/tombo/tombo_helper.py
        """
        # check both experiment type and kit slots for "rna"
        exp_type, exp_kit = None, None
        try:
            exp_type = self.bytes_to_string(self['UniqueGlobalKey/context_tags'].attrs[
                                        'experiment_type'])
            # remove the word internal since it contains rna.
            exp_type = exp_type.replace('internal', '')
        except:
            pass
        try:
            exp_kit = self.bytes_to_string(self['UniqueGlobalKey/context_tags'].attrs[
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


def iterate_fast5(path, strand_list=None, paths=False, mode='r', limit=None, files_group_pattern=None, sort_by_size=None):
    """Iterate over directory or list of .fast5 files.

    :param path: Directory in which single read fast5 are located or filename.
    :param strand_list: list of filenames to iterate, will be combined with path.
    :param paths: yield file paths instead of fast5 objects.
    :param mode: mode for opening files.
    :param limit: limit number of files to consider.
    :param files_group_pattern: yield file paths in groups of specified pattern
    :param sort_by_size: 'desc' - from largest to smallest, 'asc' - opposite
    """
    if strand_list is None:
        #  Could make glob more specific to filename pattern expected
        if os.path.isdir(path):
            files = glob(os.path.join(path, '*.fast5'))
        else:
            files = [path]
    elif os.path.isfile(strand_list):
        names = np.genfromtxt(
            strand_list, delimiter='\t', dtype=None, names=True
        )['filename']
        files = [os.path.join(path, x) for x in names]
    else:
        files = [os.path.join(path, x) for x in strand_list]

    if sort_by_size is not None:
        reverse = True if sort_by_size == 'desc' else False
        files.sort(reverse=reverse, key=lambda x: os.path.getsize(x))

    for f in files[:limit]:
        if not os.path.exists(f):
            sys.stderr.write('File {} does not exist, skipping\n'.format(f))
            continue
        if not paths:
            fh = Fast5(f, read=mode)
            yield fh
            fh.close()
        else:
            yield os.path.abspath(f)


def main():
    fast5_file = "/Users/andrewbailey/CLionProjects/nanopore-RNN/nanotensor/tests/test_files/minion-reads/canonical/miten_PC_20160820_FNFAD20259_MN17223_sequencing_run_AMS_158_R9_WGA_Ecoli_08_20_16_43623_ch100_read214_strand.fast5"
    f5fh = Fast5(fast5_file, read='r+')
    print(f5fh.raw_attributes)
    print(f5fh.get_read_stats())

if __name__ == '__main__':
    main()
