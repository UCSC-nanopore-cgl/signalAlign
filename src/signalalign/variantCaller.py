#!/usr/bin/env python3
"""Generate non-canonical nucleotide probability predictions using signal align output
"""

import sys
import os
import pandas as pd
import numpy as np
from py3helpers.utils import list_dir, merge_lists
from signalalign.signalAlignment import SignalAlignment
from signalalign.nanoporeRead import NanoporeRead
from signalalign.utils.sequenceTools import CustomAmbiguityPositions


class MarginalizeVariants(object):

    def __init__(self, variant_tsv, variants):
        """Marginalize over all posterior probabilities to give a per position read probability
        :param variants: bases to track probabilities
        :param variant_tsv: path to variantCaller tsv file
        """
        self.variant_tsv = variant_tsv
        self.read_name = os.path.basename(variant_tsv)
        assert os.path.exists(self.variant_tsv), "Variant tsv path does not exist: {}".format(variant_tsv)
        self.variant_data = SignalAlignment.read_in_signal_align_tsv(self.variant_tsv, "variantCaller")
        self.variants = sorted(variants)
        self.columns = merge_lists([['read_name', 'contig', 'position', 'strand'], list(self.variants)])
        self.contig = NanoporeRead.bytes_to_string(self.variant_data["contig"][0])
        self.position_probs = pd.DataFrame()
        self.has_data = False
        self.per_read_calls = pd.DataFrame()
        self.per_read_columns = merge_lists([['read_name', 'contig', 'strand'], list(self.variants)])

    def get_data(self):
        """Calculate the normalized probability of variant for each nucleotide and across the read"""
        # final location of per position data and per read data
        data = []
        per_read_data = []
        for mapping_strand in set(self.variant_data["forward_mapped"]):
            strand = "-"
            if mapping_strand == b"forward":
                strand = "+"

            stand_specifc_data = self.variant_data[self.variant_data["forward_mapped"] == mapping_strand]
            # get positions on strand
            positions = set(stand_specifc_data["reference_position"])
            n_positions = len(positions)
            strand_read_nuc_data = [0] * len(self.variants)

            # marginalize probabilities for each position
            for pos in positions:
                pos_data = stand_specifc_data[stand_specifc_data["reference_position"] == pos]
                total_prob = 0
                position_nuc_dict = {x: 0.0 for x in self.variants}
                # Get total probability for each nucleotide
                for nuc in set(pos_data["base"]):
                    nuc_data = pos_data[pos_data["base"] == nuc]
                    nuc_prob = max(nuc_data["posterior_probability"])
                    total_prob += nuc_prob
                    position_nuc_dict[NanoporeRead.bytes_to_string(nuc)] = nuc_prob
                # normalize probabilities over each position
                nuc_data = [0] * len(self.variants)
                for nuc in position_nuc_dict.keys():
                    index = self.variants.index(nuc)
                    nuc_data[index] = position_nuc_dict[nuc] / total_prob
                    strand_read_nuc_data[index] += nuc_data[index]
                data.append(merge_lists([[self.read_name, self.contig, pos, strand], nuc_data]))

            per_read_data.append(merge_lists([[self.read_name, self.contig, strand],
                                              [prob / n_positions for prob in strand_read_nuc_data]]))

        self.position_probs = pd.DataFrame(data, columns=self.columns)
        self.per_read_calls = pd.DataFrame(per_read_data, columns=self.per_read_columns)
        self.has_data = True

        return self.position_probs


class AggregateOverReads(object):

    def __init__(self, variant_tsv_dir, variants="ATGC"):
        """Marginalize over all posterior probabilities to give a per position read probability
        :param variant_tsv_dir: directory of variantCaller output from signalAlign
        :param variants: bases to track probabilities
        """
        self.variant_tsv_dir = variant_tsv_dir
        self.variants = sorted(variants)
        self.columns = merge_lists([['contig', 'position', 'strand'], list(self.variants)])
        self.variant_tsvs = list_dir(self.variant_tsv_dir, ext="tsv")
        self.aggregate_position_probs = pd.DataFrame()
        self.per_position_data = pd.DataFrame()
        self.per_read_data = pd.DataFrame()
        self.has_data = self._aggregate_all_variantcalls()

    def _aggregate_all_variantcalls(self):
        """Aggregate all the variant calls"""
        for v_tsv in self.variant_tsvs:
            if os.stat(v_tsv).st_size == 0:
                continue
            mv_h = MarginalizeVariants(v_tsv, variants=self.variants)
            mv_h.get_data()
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
        for contig in set(all_data["contig"]):
            contig_data = all_data[all_data["contig"] == contig]
            for strand in set(contig_data["strand"]):
                strand_data = contig_data[contig_data["strand"] == strand]
                for position in set(strand_data["position"]):
                    position_data = strand_data[strand_data["position"] == position]
                    sum_total = sum(sum(position_data.loc[:, base]) for base in self.variants)
                    normalized_probs = [sum(position_data.loc[:, base]) / sum_total for base in self.variants]
                    yield merge_lists([[contig, position, strand], normalized_probs])

    def write_data(self, out_path):
        """Write out aggregate_position_probs to tsv file"""
        self.aggregate_position_probs.to_csv(out_path, sep='\t', index=False)

    def generate_labels(self, labelled_positions, predicted_data):
        """Generate labels for predictions given labelled positions"""
        for char in self.variants:
            predicted_data.loc[:, char+"_label"] = pd.Series(0, index=predicted_data.index)

        for i in range(len(predicted_data)):
            contig = predicted_data.loc[i]["contig"]
            strand = predicted_data.loc[i]["strand"]
            position = predicted_data.loc[i]["position"]
            true_char = get_true_character(labelled_positions, contig, strand, position)
            if true_char is None:
                print("No variant found in labelled data at chr:{} pos:{} strand:{}: Check positions file".format(contig, position, strand))
            else:
                predicted_data.loc[i, true_char+"_label"] = 1

        return predicted_data


def get_true_character(true_positions_data, contig, strand, position):
    """Get true character from an positions/ambiguity file"""
    true_char = true_positions_data.loc[(true_positions_data['contig'] == contig) &
                                        (true_positions_data['strand'] == strand) &
                                        (true_positions_data['position'] == position)]["change_to"]
    if len(true_char) > 0:
        return true_char.values[0]
    else:
        return None


def create_labels_from_positions_file(positions_file_path, variants):
    """Read in positions file and add classification labels
    :param positions_file_path: path to a tsv positions file with the following columns
                        contig  position    strand  change_from change_to

    :param variants: string of characters to expect from change_to
    :return: data with new columns
    """
    data = CustomAmbiguityPositions.parseAmbiguityFile(positions_file_path)
    for char in variants:
        data.loc[:, char] = pd.Series(data["change_to"] == char, index=data.index)
    return data



def main():
    variant_caller_tsvs = "/Users/andrewbailey/data/m6a_rna_data/tempFiles_alignment/m6a"
    all_calls = list_dir(variant_caller_tsvs, ext="tsv")
    variants = "AF"
    labeled_positions_file = "/Users/andrewbailey/data/m6a_rna_data/positions_file.tsv"
    # margin_a = MarginalizeVariants(all_calls[0], variants="AF")
    # print(margin_a.get_data())
    data = create_labels_from_positions_file(labeled_positions_file, variants=variants)
    aor_h = AggregateOverReads(variant_caller_tsvs, variants)
    data2 = aor_h.marginalize_over_all_reads()
    data3 = aor_h.aggregate_all_variantcalls()

    data4 = aor_h.generate_labels(data, data2)
    # labels = data[list(variants)]
    # probs1 = data2[list(variants)]
    # probs2 = data3[list(variants)]
    for x in data4.iterrows():
        print(x)

if __name__ == "__main__":
    sys.exit(main())
