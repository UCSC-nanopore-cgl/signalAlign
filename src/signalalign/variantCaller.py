#!/usr/bin/env python3
"""Generate non-canonical nucleotide probability predictions using signal align output
"""

import os

import numpy as np
import pandas as pd
from py3helpers.utils import list_dir, merge_lists
from signalalign.nanoporeRead import NanoporeRead
from signalalign.signalAlignment import SignalAlignment
from signalalign.train.trainModels import read_in_alignment_file
from signalalign.utils.sequenceTools import CustomAmbiguityPositions


class MarginalizeVariants(object):

    def __init__(self, variant_data, variants, read_name):
        """Marginalize over all posterior probabilities to give a per position read probability
        :param variants: bases to track probabilities
        :param variant_data: variant data
        """
        self.read_name = read_name
        self.variant_data = variant_data
        self.variants = sorted(variants)
        self.columns = merge_lists([['read_name', 'contig', 'position', 'strand', 'forward_mapped'],
                                    list(self.variants)])
        self.contig = NanoporeRead.bytes_to_string(self.variant_data["contig"][0])
        self.position_probs = pd.DataFrame()
        self.has_data = False
        self.per_read_calls = pd.DataFrame()
        self.per_read_columns = merge_lists([['read_name', 'contig', 'strand', "forward_mapped",
                                              "n_sites"], list(self.variants)])

    def get_data(self):
        """Calculate the normalized probability of variant for each nucleotide and across the read"""
        # final location of per position data and per read data
        data = []
        per_read_data = []
        for read_strand in (b"t", b"c"):
            read_strand_specifc_data = self.variant_data[self.variant_data["strand"] == read_strand]
            read_strand = read_strand.decode("utf-8")
            if len(read_strand_specifc_data) == 0:
                continue
            for forward_mapped in set(self.variant_data["forward_mapped"]):
                mapping_strand = "-"
                if forward_mapped == b"forward":
                    mapping_strand = "+"
                strand_specifc_data = read_strand_specifc_data[read_strand_specifc_data["forward_mapped"] ==
                                                               forward_mapped]
                if len(strand_specifc_data) == 0:
                    continue
                # get positions on strand
                positions = set(strand_specifc_data["reference_position"])
                n_positions = len(positions)
                strand_read_nuc_data = [0] * len(self.variants)

                # marginalize probabilities for each position
                for pos in positions:
                    pos_data = strand_specifc_data[strand_specifc_data["reference_position"] == pos]
                    total_prob = 0
                    position_nuc_dict = {x: 0.0 for x in self.variants}
                    # Get total probability for each nucleotide
                    for nuc in set(pos_data["base"]):
                        nuc_data = pos_data[pos_data["base"] == nuc]
                        nuc_prob = sum(nuc_data["posterior_probability"])
                        total_prob += nuc_prob
                        position_nuc_dict[NanoporeRead.bytes_to_string(nuc)] = nuc_prob
                    # normalize probabilities over each position
                    nuc_data = [0] * len(self.variants)
                    for nuc in position_nuc_dict.keys():
                        index = self.variants.index(nuc)
                        nuc_data[index] = position_nuc_dict[nuc] / total_prob
                        strand_read_nuc_data[index] += nuc_data[index]

                    data.append(merge_lists([[self.read_name, self.contig, pos, read_strand, mapping_strand],
                                             nuc_data]))
                if n_positions > 0:
                    per_read_data.append(merge_lists([[self.read_name, self.contig, read_strand, mapping_strand,
                                                       n_positions],
                                                      [prob / n_positions for prob in strand_read_nuc_data]]))

            self.position_probs = pd.DataFrame(data, columns=self.columns)
            self.per_read_calls = pd.DataFrame(per_read_data, columns=self.per_read_columns)
            self.has_data = True

        return self.position_probs


class MarginalizeFullVariants(object):

    def __init__(self, full_data, variants, read_name, forward_mapped):
        """Marginalize over all posterior probabilities to give a per position read probability
        :param variants: bases to track probabilities
        :param full_data: path to full tsv file

                             ['contig', 'reference_index',
                              'reference_kmer', 'read_file',
                              'strand', 'event_index',
                              'event_mean', 'event_noise',
                              'event_duration', 'aligned_kmer',
                              'scaled_mean_current', 'scaled_noise',
                              'posterior_probability', 'descaled_event_mean',
                              'ont_model_mean', 'path_kmer']
        """
        self.read_name = read_name
        self.full_data = full_data
        self.variant_data = self.full_data[["X" in kmer for kmer in self.full_data["reference_kmer"]]]
        self.variants = sorted(variants)
        self.forward_mapped = forward_mapped
        self.columns = merge_lists([['read_name', 'contig', 'position', 'strand', 'forward_mapped'],
                                    list(self.variants)])
        self.contig = NanoporeRead.bytes_to_string(self.full_data["contig"][0])
        self.position_probs = pd.DataFrame()
        self.has_data = False
        self.per_read_calls = pd.DataFrame()
        self.per_read_columns = merge_lists([['read_name', 'contig', 'strand', "forward_mapped", "n_sites"],
                                             list(self.variants)])

    def get_data(self):
        """Calculate the normalized probability of variant for each nucleotide and across the read"""
        # final location of per position data and per read data
        data = []
        per_read_data = []
        if self.forward_mapped:
            mapping_strands = ["+", "-"]
        else:
            mapping_strands = ["-", "+"]

        if len(self.variant_data) > 0:
            kmer_len_1 = len(self.variant_data["reference_kmer"].iloc[0]) - 1
            mapping_index = 0
            for read_strand in ("t", "c"):
                read_strand_specifc_data = self.variant_data[self.variant_data["strand"] == read_strand]
                # read_strand = read_strand.decode("utf-8")
                if len(read_strand_specifc_data) == 0:
                    continue
                # get positions on strand
                positions = sorted(set(read_strand_specifc_data["reference_index"]))

                if mapping_strands[mapping_index] == "-":
                    positions = positions[::-1]

                strand_read_nuc_data = [0] * len(self.variants)

                # marginalize probabilities for each position
                n_positions = 0
                for pos in positions:
                    pos_data = read_strand_specifc_data[read_strand_specifc_data["reference_index"] == pos]
                    if pos_data["aligned_kmer"].iloc[0][kmer_len_1] != "X":
                        continue
                    n_positions += 1
                    total_prob = 0
                    position_nuc_dict = {x: 0.0 for x in self.variants}
                    # Get total probability for each nucleotide
                    for nuc in self.variants:
                        # kmer_len_1 = pos_data["reference_kmer"].iloc[0].find("X")
                        # print(pos_data["reference_kmer"].iloc[0])
                        nuc_data = pos_data[[nuc == kmer[kmer_len_1] for kmer in pos_data["path_kmer"]]]
                        nuc_prob = sum(nuc_data["posterior_probability"])
                        total_prob += nuc_prob
                        position_nuc_dict[NanoporeRead.bytes_to_string(nuc)] = nuc_prob
                    # normalize probabilities over each position
                    nuc_data = [0] * len(self.variants)
                    for index, nuc in enumerate(self.variants):
                        assert total_prob > 0, "Check 'variants' parameter. There seems to be no kmers with those " \
                                               "variant characters"
                        nuc_data[index] = position_nuc_dict[nuc] / total_prob
                        strand_read_nuc_data[index] += nuc_data[index]
                    data.append(merge_lists([[self.read_name, self.contig, pos, read_strand,
                                              mapping_strands[mapping_index]], nuc_data]))
                if n_positions > 0:
                    per_read_data.append(merge_lists([[self.read_name, self.contig, read_strand,
                                                       mapping_strands[mapping_index], n_positions],
                                                      [prob / n_positions for prob in strand_read_nuc_data]]))
                mapping_index += 1
            self.position_probs = pd.DataFrame(data, columns=self.columns)
            self.per_read_calls = pd.DataFrame(per_read_data, columns=self.per_read_columns)
            self.has_data = True

        else:
            self.has_data = False

        return self.position_probs


class AggregateOverReads(object):

    def __init__(self, variant_tsv_dir, variants="ATGC", verbose=False):
        """Marginalize over all posterior probabilities to give a per position read probability
        :param variant_tsv_dir: directory of variantCaller output from signalAlign
        :param variants: bases to track probabilities
        """
        self.variant_tsv_dir = variant_tsv_dir
        self.variants = sorted(variants)
        self.columns = merge_lists([['contig', 'position', 'strand', 'forward_mapped'], list(self.variants)])
        self.variant_tsvs = list_dir(self.variant_tsv_dir, ext=".vc.tsv")
        self.aggregate_position_probs = pd.DataFrame()
        self.per_position_data = pd.DataFrame()
        self.verbose = verbose
        self.per_read_data = pd.DataFrame()
        self.has_data = self._aggregate_all_variantcalls()

    def _aggregate_all_variantcalls(self):
        """Aggregate all the variant calls"""
        for v_tsv in self.variant_tsvs:
            if os.stat(v_tsv).st_size == 0:
                continue
            read_name = os.path.basename(v_tsv)
            variant_data = SignalAlignment.read_in_signal_align_tsv(v_tsv, "variantCaller")
            mv_h = MarginalizeVariants(variant_data, variants=self.variants, read_name=read_name)
            mv_h.get_data()
            if self.verbose:
                print(v_tsv)
            self.per_position_data = self.per_position_data.append(mv_h.position_probs, ignore_index=True)
            self.per_read_data = self.per_read_data.append(mv_h.per_read_calls, ignore_index=True)

        return True

    def marginalize_over_all_reads(self):
        """Calculate the per position posterior probability"""
        assert self.has_data, "AggregateOverReads does not have data. Make sure you initialized correctly"
        self.aggregate_position_probs = pd.concat([pd.DataFrame([i], columns=self.columns)
                                                   for i in self._normalize_all_data(self.per_position_data)],
                                                  ignore_index=True)
        return self.aggregate_position_probs

    def _normalize_all_data(self, all_data):
        """Helper function to normalize all probability data"""
        for strand in set(all_data["strand"]):
            strand_data = all_data[all_data["strand"] == strand]
            for contig in set(strand_data["contig"]):
                contig_data = strand_data[strand_data["contig"] == contig]
                for mapped_strand in set(contig_data["forward_mapped"]):
                    strand_mapped_data = contig_data[contig_data["forward_mapped"] == mapped_strand]
                    for position in set(strand_mapped_data["position"]):
                        position_data = strand_mapped_data[strand_mapped_data["position"] == position]
                        sum_total = sum(sum(position_data.loc[:, base]) for base in self.variants)
                        normalized_probs = [np.round(sum(position_data.loc[:, base]) / sum_total, 6) for base
                                            in self.variants]
                        yield merge_lists([[contig, position, strand, mapped_strand], normalized_probs])

    def write_data(self, out_path):
        """Write out aggregate_position_probs to tsv file"""
        self.aggregate_position_probs.to_csv(out_path, sep='\t', index=False)

    def generate_labels(self, labelled_positions, predicted_data):
        """Generate labels for predictions given labelled positions.
        Note: This will drop sites that do not have labels in 'labelled_positions'
        """
        for char in self.variants:
            predicted_data.loc[:, char+"_label"] = pd.Series(0, index=predicted_data.index)

        for i in range(len(predicted_data)):
            contig = predicted_data.loc[i]["contig"]
            forward_mapped = predicted_data.loc[i]["forward_mapped"]
            position = predicted_data.loc[i]["position"]
            true_char = get_true_character(labelled_positions, contig, forward_mapped, position)
            if true_char is None:
                print("No variant found in labelled data at chr:{} pos:{} "
                      "forward_mapped:{}: Check positions file".format(contig, position, forward_mapped))
                predicted_data = predicted_data.drop([i])
            else:
                predicted_data.loc[i, true_char+"_label"] = 1

        return predicted_data

    def generate_labels2(self, predicted_data, true_char):
        """Generate labels for predictions given labelled positions"""
        for char in self.variants:
            predicted_data.loc[:, char+"_label"] = pd.Series(0, index=predicted_data.index)

        for i in range(len(predicted_data)):
            predicted_data.loc[i, true_char+"_label"] = 1

        return predicted_data


class AggregateOverReadsFull(object):

    def __init__(self, sa_full_tsv_dir, variants="ATGC", verbose=False):
        """Marginalize over all posterior probabilities to give a per position read probability
        :param sa_full_tsv_dir: directory of full output from signalAlign
        :param variants: bases to track probabilities
        """
        self.sa_full_tsv_dir = sa_full_tsv_dir
        self.variants = sorted(variants)
        self.columns = merge_lists([['contig', 'position', 'strand', 'forward_mapped'], list(self.variants)])
        self.forward_tsvs = list_dir(self.sa_full_tsv_dir, ext=".forward.tsv")
        self.backward_tsvs = list_dir(self.sa_full_tsv_dir, ext=".backward.tsv")
        self.verbose = verbose

        self.aggregate_position_probs = pd.DataFrame()
        self.per_position_data = pd.DataFrame()
        self.per_read_data = pd.DataFrame()
        self.has_data = self._aggregate_all_variantcalls()

    def _aggregate_all_variantcalls(self):
        """Aggregate all the variant calls"""
        for v_tsv in self.forward_tsvs:
            if os.stat(v_tsv).st_size == 0:
                continue
            read_name = os.path.basename(v_tsv)
            variant_data = read_in_alignment_file(v_tsv)
            mv_h = MarginalizeFullVariants(variant_data, variants=self.variants, read_name=read_name,
                                           forward_mapped=True)
            mv_h.get_data()
            if self.verbose:
                print(v_tsv)
            if mv_h.has_data:
                self.per_position_data = self.per_position_data.append(mv_h.position_probs, ignore_index=True)
                self.per_read_data = self.per_read_data.append(mv_h.per_read_calls, ignore_index=True)

        for v_tsv in self.backward_tsvs:
            if os.stat(v_tsv).st_size == 0:
                continue
            read_name = os.path.basename(v_tsv)
            variant_data = read_in_alignment_file(v_tsv)
            mv_h = MarginalizeFullVariants(variant_data, variants=self.variants, read_name=read_name,
                                           forward_mapped=False)
            mv_h.get_data()
            if self.verbose:
                print(v_tsv)
            if mv_h.has_data:
                self.per_position_data = self.per_position_data.append(mv_h.position_probs, ignore_index=True)
                self.per_read_data = self.per_read_data.append(mv_h.per_read_calls, ignore_index=True)

        return True

    def marginalize_over_all_reads(self):
        """Calculate the per position posterior probability"""
        assert self.has_data, "AggregateOverReads does not have data. Make sure you initialized correctly"
        self.aggregate_position_probs = pd.concat([pd.DataFrame([i], columns=self.columns)
                                                   for i in self._normalize_all_data(self.per_position_data)],
                                                  ignore_index=True)
        return self.aggregate_position_probs

    def _normalize_all_data(self, all_data):
        """Helper function to normalize all probability data"""
        for strand in set(all_data["strand"]):
            strand_data = all_data[all_data["strand"] == strand]
            for contig in set(strand_data["contig"]):
                contig_data = strand_data[strand_data["contig"] == contig]
                for mapped_strand in set(contig_data["forward_mapped"]):
                    strand_mapped_data = contig_data[contig_data["forward_mapped"] == mapped_strand]
                    for position in set(strand_mapped_data["position"]):
                        position_data = strand_mapped_data[strand_mapped_data["position"] == position]
                        sum_total = sum(sum(position_data.loc[:, base]) for base in self.variants)
                        normalized_probs = [np.round(sum(position_data.loc[:, base]) / sum_total, 6) for base in
                                            self.variants]
                        yield merge_lists([[contig, position, strand, mapped_strand], normalized_probs])

    def write_data(self, out_path):
        """Write out aggregate_position_probs to tsv file"""
        self.aggregate_position_probs.to_csv(out_path, sep='\t', index=False)

    def generate_labels(self, labelled_positions, predicted_data):
        """Generate labels for predictions given labelled positions.
        Note: This will drop sites that do not have labels in 'labelled_positions'
        """
        for char in self.variants:
            predicted_data.loc[:, char+"_label"] = pd.Series(0, index=predicted_data.index)

        for i in range(len(predicted_data)):
            contig = predicted_data.loc[i]["contig"]
            forward_mapped = predicted_data.loc[i]["forward_mapped"]
            position = predicted_data.loc[i]["position"]
            true_char = get_true_character(labelled_positions, contig, forward_mapped, position)
            if true_char is None:
                print("No variant found in labelled data at chr:{} pos:{} forward_mapped:{}: "
                      "Check positions file".format(contig, position, forward_mapped))
                predicted_data = predicted_data.drop([i])
            else:
                predicted_data.loc[i, true_char+"_label"] = 1

        return predicted_data

    def generate_labels2(self, predicted_data, true_char):
        """Generate labels for predictions given labelled positions"""
        for char in self.variants:
            if char == true_char:
                predicted_data.loc[:, char+"_label"] = pd.Series(1, index=predicted_data.index)
            else:
                predicted_data.loc[:, char+"_label"] = pd.Series(0, index=predicted_data.index)

        return predicted_data


def get_true_character(true_positions_data, contig, forward_mapped, position):
    """Get true character from an positions/ambiguity file"""
    true_char = true_positions_data.loc[(true_positions_data['contig'] == contig) &
                                        (true_positions_data['strand'] == forward_mapped) &
                                        (true_positions_data['position'] == position)]["change_to"]
    if len(true_char) > 0:
        return true_char.values[0]
    else:
        return None


def create_labels_from_positions_file(positions_file_path, variants):
    """Read in positions file and add classification labels
    :param positions_file_path: path to a tsv positions file with the following columns
                        contig  position    forward_mapped  change_from change_to

    :param variants: string of characters to expect from change_to
    :return: data with new columns
    """
    data = CustomAmbiguityPositions.parseAmbiguityFile(positions_file_path)
    for char in variants:
        data.loc[:, char] = pd.Series(data["change_to"] == char, index=data.index)
    return data
