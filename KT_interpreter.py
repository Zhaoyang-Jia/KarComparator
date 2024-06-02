from Report_Genes import *
# from reportlab.pdfgen import canvas
# from reportlab.lib import colors
# from reportlab.lib.pagesizes import letter
# from reportlab.lib.styles import getSampleStyleSheet
# from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image
from read_OMKar_output import *
from datetime import datetime
import subprocess
import os


class Aligned_Haplotype:
    id: int  # for debugging
    chrom: str
    alignment_score: int
    mt_aligned: [str]
    wt_aligned: [str]
    block_indices: {int: (int, int)}  # [block_id: (a, b)] where event block are [a, b)
    mt_blocks: {int: (str,)}  # insertion discordant blocks
    wt_blocks: {int: (str,)}  # deletion discordant blocks
    concordant_blocks: {int: (str,)}
    discordant_block_assignment: {int: str}  # block index: event name, exactly one event can be assigned to one block; only for discordant blocks
    size_dict: {str: int}  # seg_name: bp size

    def __init__(self, mt_aligned, wt_aligned, alignment_score, size_dict, id, chrom, mt_hap):
        self.id = id
        self.chrom = chrom
        self.mt_aligned = mt_aligned
        self.wt_aligned = wt_aligned
        self.alignment_score = alignment_score
        self.size_dict = size_dict
        self.mt_hap = mt_hap

        ## get event_blocks
        event_block = []
        block_types = []  # [str], a block is either ins/del
        # we know len(wt_hap) == len(mt_hap) as they are global alignments
        hap_len = len(self.mt_aligned)
        seg_idx = 0
        while seg_idx < hap_len:
            mt_seg = self.mt_aligned[seg_idx]
            wt_seg = self.wt_aligned[seg_idx]
            if mt_seg == wt_seg:
                seg_idx += 1
                continue
            else:
                # this will be the first seg of (potentially) a section, that is indel
                if wt_seg == '-':
                    # insertion
                    # block ends whenever continuous section OR insertion section ends
                    continuous_end = continuous_extension(self.mt_aligned, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if self.wt_aligned[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('ins')
                elif mt_seg == '-':
                    # deletion
                    continuous_end = continuous_extension(self.wt_aligned, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if self.mt_aligned[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('del')
                else:
                    raise ValueError('mismatch not allowed')
                event_block.append((seg_idx, block_end))
                seg_idx = block_end

        ## extract concordant blocks, mt/wt blocks
        discordant_blocks = {event_block_boundaries: block_types[event_block_idx] for event_block_idx, event_block_boundaries in enumerate(event_block)}
        # all gaps in discordant blocks are concordant blocks
        self.block_indices = {}
        self.wt_blocks = {}
        self.mt_blocks = {}
        self.concordant_blocks = {}
        self.discordant_block_assignment = {}
        block_id = 0
        p_endpoint = 0
        for discordant_block in discordant_blocks:
            c_startpoint = int(discordant_block[0])
            if p_endpoint < c_startpoint:
                # fill gap with concordant blocks
                # between two discordant blocks, the region belongs to at most 1 concordant block because the wt is always continuous
                # (i.e. will not split into 2 blocks)
                self.block_indices[block_id] = (p_endpoint, c_startpoint)
                self.concordant_blocks[block_id] = wt_aligned[p_endpoint: c_startpoint]
                block_id += 1
            self.block_indices[block_id] = discordant_block
            self.discordant_block_assignment[block_id] = ''
            if discordant_blocks[discordant_block] == 'del':
                self.wt_blocks[block_id] = wt_aligned[discordant_block[0]: discordant_block[1]]
            elif discordant_blocks[discordant_block] == 'ins':
                self.mt_blocks[block_id] = mt_aligned[discordant_block[0]: discordant_block[1]]
            block_id += 1
            p_endpoint = discordant_block[1]
        # add the last concordant block, if present
        if p_endpoint < len(self.wt_aligned):
            self.block_indices[block_id] = (p_endpoint, len(self.wt_aligned))
            self.concordant_blocks[block_id] = wt_aligned[p_endpoint: len(self.wt_aligned)]

    def __str__(self):
        return "ID<{}>".format(self.id)

    def get_block_segs_and_block_type(self, block_idx):
        """
        given block_idx, get block segs and block type
        :param block_idx:
        :return:
        """
        if block_idx in self.wt_blocks:
            return self.wt_blocks[block_idx], self.wt_blocks
        elif block_idx in self.mt_blocks:
            return self.mt_blocks[block_idx], self.mt_blocks
        elif block_idx in self.concordant_blocks:
            return self.concordant_blocks[block_idx], self.concordant_blocks
        else:
            raise RuntimeError('block idx not found in the Haplotype')

    def get_closest_mt_block(self, block_idx, search_direction):
        """
        real block only exists in mt-block and conc-block
        :param block_idx:
        :param search_direction:
        :return:
        """
        if search_direction == 'left':
            nearest_idx_found = -1
            return_indexed_seg = ''
            for block_idx_in_mt, indexed_seg in self.mt_blocks.items():
                if block_idx_in_mt < block_idx:
                    # on the leftside
                    if block_idx_in_mt > nearest_idx_found:
                        nearest_idx_found = block_idx_in_mt
                        return_indexed_seg = indexed_seg
            for block_idx_in_conc, indexed_seg in self.concordant_blocks.items():
                if block_idx_in_conc < block_idx:
                    # on the leftside
                    if block_idx_in_conc > nearest_idx_found:
                        nearest_idx_found = block_idx_in_conc
                        return_indexed_seg = indexed_seg
            if nearest_idx_found == -1:
                return_indexed_seg = 'p-ter'
            else:
                return_indexed_seg = return_indexed_seg[-1]
        elif search_direction == 'right':
            nearest_idx_found = 9999
            return_indexed_seg = ''
            for block_idx_in_mt, indexed_seg in self.mt_blocks.items():
                if block_idx_in_mt > block_idx:
                    # on the rightside
                    if block_idx_in_mt < nearest_idx_found:
                        nearest_idx_found = block_idx_in_mt
                        return_indexed_seg = indexed_seg
            for block_idx_in_conc, indexed_seg in self.concordant_blocks.items():
                if block_idx_in_conc > block_idx:
                    # on the rightside
                    if block_idx_in_conc < nearest_idx_found:
                        nearest_idx_found = block_idx_in_conc
                        return_indexed_seg = indexed_seg
            if nearest_idx_found == 9999:
                return_indexed_seg = 'q-ter'
            else:
                return_indexed_seg = return_indexed_seg[0]
        else:
            raise ValueError('direction illegal')
        return return_indexed_seg

    def report_SV(self, event_blocks, event_types):
        """
        :param event_blocks: cumulative dict, where value is a list of event_block, this helps to group paired-events that were on different chromosomes
        :param event_types: cumulative dict, where value is a list of event_type, same reason as above
        :return:
        """

        def get_block_boundaries(block_idx):
            if block_idx == 0:
                # first block
                left = 'p-ter'
            else:
                # previous_block, _ = self.get_block_segs_and_block_type(block_idx - 1)
                left = self.get_closest_mt_block(block_idx, 'left')
                # left = previous_block[-1]
            if block_idx == len(self.block_indices) - 1:
                # last block
                right = 'q-ter'
            else:
                # next_block, _ = self.get_block_segs_and_block_type(block_idx + 1)
                right = self.get_closest_mt_block(block_idx, 'right')
                # right = next_block[0]
            return left, right

        for mt_block_idx, mt_block in self.mt_blocks.items():
            mt_block_str = ','.join(list(mt_block))
            event_id = int(self.discordant_block_assignment[mt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[mt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(mt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.{}.mt({}).{}.{}'.format(self.id, mt_block_idx, mt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.{}.mt({}).{}.{}'.format(self.id, mt_block_idx, mt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]
        for wt_block_idx, wt_block in self.wt_blocks.items():
            wt_block_str = ','.join(list(wt_block))
            event_id = int(self.discordant_block_assignment[wt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[wt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(wt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.{}.wt({}).{}.{}'.format(self.id, wt_block_idx, wt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.{}.wt({}).{}.{}'.format(self.id, wt_block_idx, wt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]

        # sort each event_blocks by chr, then by block_idx
        for event_id, event_block_list in event_blocks.items():
            event_blocks[event_id] = sorted(event_block_list, key=lambda x: (int(x.split('.')[0]), int(x.split('.')[1])), reverse=False)

        return event_blocks, event_types

    def search_seed(self, query_section, seed_type, d, eps):
        """
        :param query_section:
        :param seed_type: 'ins' OR 'del', which type of seed to be searched of
        :param d:
        :param eps:
        :return: block_idx if found, -1 if not
        """
        if seed_type == 'del':
            query_block = self.wt_blocks
        elif seed_type == 'ins':
            query_block = self.mt_blocks
        else:
            raise ValueError('invalid seed type')
        for c_block_idx, c_block in query_block.items():
            if len(self.discordant_block_assignment[c_block_idx]) != 0:
                # already has assignment
                # TODO: can improve this by bipartite matching, instead of greedy matching
                continue
            seed_start, seed_end = is_seeded(c_block, query_section, self.size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                return c_block_idx
        return -1


def lcs(list1, list2, size_dict):
    """
    longest common subsequence finding
    :param list1:
    :param list2:
    :param size_dict:
    :return:
    """
    indel_penalty_per_nt = -1
    alignment_1 = []
    alignment_2 = []
    scoring_matrix = [[0 for i in range(0, len(list2) + 1)] for j in range(0, len(list1) + 1)]
    backtrack_matrix = [["" for i in range(0, len(list2) + 1)] for j in range(0, len(list1) + 1)]

    # initialize starting grid
    scoring_matrix[0][0] = 0
    for row_index in range(1, len(list1) + 1):
        current_segment = list1[row_index - 1]
        len_current_segment = size_dict[current_segment[:-1]]  # remove sign
        scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0] + len_current_segment * indel_penalty_per_nt
        backtrack_matrix[row_index][0] = "down"

    for col_index in range(1, len(list2) + 1):
        current_segment = list2[col_index - 1]
        len_current_segment = size_dict[current_segment[:-1]]  # remove sign
        scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1] + len_current_segment * indel_penalty_per_nt
        backtrack_matrix[0][col_index] = "rigt"

    # DP recursion
    for row_index in range(1, len(list1) + 1):
        for col_index in range(1, len(list2) + 1):
            len_down_segment = size_dict[list1[row_index - 1][:-1]]
            down_value = scoring_matrix[row_index - 1][col_index] + len_down_segment * indel_penalty_per_nt

            len_right_segment = size_dict[list2[col_index - 1][:-1]]
            right_value = scoring_matrix[row_index][col_index - 1] + len_right_segment * indel_penalty_per_nt

            if list1[row_index - 1] != list2[col_index - 1]:
                # mismatch: not allowed
                diagonal_value = float('-inf')
            else:
                # match
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1]

            if diagonal_value >= down_value and diagonal_value >= right_value:
                scoring_matrix[row_index][col_index] = diagonal_value
                backtrack_matrix[row_index][col_index] = "diag"
            elif down_value >= right_value:
                scoring_matrix[row_index][col_index] = down_value
                backtrack_matrix[row_index][col_index] = "down"
            else:
                scoring_matrix[row_index][col_index] = right_value
                backtrack_matrix[row_index][col_index] = "rigt"

    # backtracking
    final_score = scoring_matrix[len(list1)][len(list2)]
    current_row = len(list1)
    current_col = len(list2)

    while True:
        if current_row == 0 and current_col == 0:
            break
        if backtrack_matrix[current_row][current_col] == "diag":
            alignment_1.insert(0, list1[current_row - 1])
            alignment_2.insert(0, list2[current_col - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, list1[current_row - 1])
            alignment_2.insert(0, "-")
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "rigt":
            alignment_1.insert(0, "-")
            alignment_2.insert(0, list2[current_col - 1])
            current_col -= 1
        else:
            raise RuntimeError("error in backtrack matrix")

    return final_score, alignment_1, alignment_2


def interpret_haplotypes(mt_hap_list: [[str]], wt_hap_list: [[str]], chrom_identities: [str], segment_size_dict: {str: int}, d=4000000, eps=4000000):
    """
    Input a haplotype and a WT, report SVs; hap from mt list must correspond to hap from wt list
    :param segment_size_dict: {str: int} segment mapped to size of the segment (eg. 11: 12,345)
    :param mt_hap_list:
    :param wt_hap_list:
    :param chrom_identities:
    :return:
    """
    event_id = 0  # static variable, keep different events separate, keep paired events under same id (eg. balanced trans)
    hap_id = 0  # used for naming haplotypes, preserves the order they come-in in mt_hap/wt_hap lists
    aligned_haplotypes = []
    for idx, mt_hap in enumerate(mt_hap_list):
        wt_hap = wt_hap_list[idx]
        c_chrom = chrom_identities[idx]
        a1, a2, a3 = lcs(mt_hap, wt_hap, segment_size_dict)
        aligned_haplotypes.append(Aligned_Haplotype(a2, a3, a1, segment_size_dict, hap_id, c_chrom, mt_hap))
        hap_id += 1

    def genomewide_seed_search(query_section, query_block_idx, query_hap_idx, query_type):
        """
        :param query_section: (segs) to be searched for
        :param query_block_idx: needed for writing discordant_block_assignment if seed found
        :param query_hap_idx: needed to prioritize searching inter-chr first, then intra-chr if inter-chr not found
        :param query_type:
        :return:
        """
        seed_found = False
        for aligned_hap_idx1, aligned_hap1 in enumerate(aligned_haplotypes):
            if aligned_hap_idx1 == query_hap_idx:
                # first get inter-chr event
                continue
            seeded_block_idx = aligned_hap1.search_seed(query_section, query_type, d=d, eps=eps)
            if seeded_block_idx != -1:
                aligned_hap1.discordant_block_assignment[seeded_block_idx] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[query_block_idx] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
                break
        if not seed_found:
            # at last, check intra-chr event
            seeded_block_idx = aligned_haplotypes[query_hap_idx].search_seed(query_section, query_type, d=d, eps=eps)
            if seeded_block_idx != -1:
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[seeded_block_idx] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[query_block_idx] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
        if seed_found:
            return True
        else:
            return False

    ### Order of resolving: balanced-translocation, inv, dup-inv, tandem-dup, del, unbalanced-translocation
    ## resolve all balanced translocations
    for aligned_hap_idx, aligned_hap in enumerate(aligned_haplotypes):
        for c_wt_block_idx, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_idx]) > 0:
                # already has assignment
                continue
            aligned_hap.discordant_block_assignment[
                c_wt_block_idx] = 'under investigation'  # block it from being matched in intra-chr event; maybe not necessary
            balanced_translocation_found = genomewide_seed_search(c_wt_block, c_wt_block_idx, aligned_hap_idx, 'ins')
            if balanced_translocation_found:
                event_id += 1
            else:
                # reset block assignment for query
                aligned_hap.discordant_block_assignment[c_wt_block_idx] = ''
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                # already has assignment
                continue
            aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'under investigation'  # block it from being matched in intra-chr event
            balanced_translocation_found = genomewide_seed_search(c_mt_block, c_mt_block_idx, aligned_hap_idx, 'del')
            if balanced_translocation_found:
                event_id += 1
            else:
                # reset block assignment for query
                aligned_hap.discordant_block_assignment[c_mt_block_idx] = ''

    ## inversion: mt{- k-} wt{k+ -} OR mt{k- -} wt{- k+}
    for aligned_hap in aligned_haplotypes:
        # case: mt{k- -} wt{- k+}, ins of inverted-seg + del of uninverted-seg
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if c_mt_block[0][-1] == "+":
                # not inverted
                continue
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            if c_mt_block_idx + 1 >= len(aligned_hap.block_indices) or \
                    (c_mt_block_idx + 1) not in aligned_hap.wt_blocks or \
                    len(aligned_hap.discordant_block_assignment[c_mt_block_idx + 1]) > 0:
                # this is the last block, next block is not a del-block, OR next del-block has assignment
                continue
            uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
            next_del_block = aligned_hap.wt_blocks[c_mt_block_idx + 1]
            seed_start, seed_end = is_seeded(next_del_block, uninverted_segs, segment_size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'inversion,{}'.format(event_id)
                aligned_hap.discordant_block_assignment[c_mt_block_idx + 1] = 'inversion,{}'.format(event_id)
                event_id += 1
        # case: mt{- k-} wt{k+ -}, del of uninverted-seg + ins of inverted-seg
        for c_wt_block_idx, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_idx]) > 0:
                continue
            if c_wt_block_idx + 1 >= len(aligned_hap.block_indices) or \
                    (c_wt_block_idx + 1) not in aligned_hap.mt_blocks or \
                    len(aligned_hap.discordant_block_assignment[c_wt_block_idx + 1]) > 0:
                # this is the last block, next block is not a ins-block, OR next ins-block has assignment
                continue
            inverted_segs = [seg[:-1] + '-' for seg in c_wt_block[::-1]]
            next_ins_block = aligned_hap.mt_blocks[c_wt_block_idx + 1]
            seed_start, seed_end = is_seeded(next_ins_block, inverted_segs, segment_size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                aligned_hap.discordant_block_assignment[c_wt_block_idx] = 'inversion,{}'.format(event_id)
                aligned_hap.discordant_block_assignment[c_wt_block_idx + 1] = 'inversion,{}'.format(event_id)
                event_id += 1

    ## duplication inversion
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            # left-dup-inv: ins of inverted-seg + concordant block w/ uninverted-seg as seed
            if c_mt_block[0][-1] == "+":
                # not inverted
                continue
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
            # this is not the last block, AND next block is concordant
            if c_mt_block_idx + 1 < len(aligned_hap.block_indices) and \
                    (c_mt_block_idx + 1) in aligned_hap.concordant_blocks:
                next_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx + 1]
                seed_start, seed_end = is_seeded(next_concordant_block, uninverted_segs, segment_size_dict, indel_direction='left', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'left_duplication_inversion,{}'.format(event_id)
                    event_id += 1
                    continue
            # right-dup-inv: concordant block w/ uninverted-seg as seed + ins of inverted-seg
            # this is not the first block, AND previous block is concordant
            if c_mt_block_idx - 1 >= 0 and \
                    (c_mt_block_idx - 1) in aligned_hap.concordant_blocks:
                previous_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx - 1]
                seed_start, seed_end = is_seeded(previous_concordant_block, uninverted_segs, segment_size_dict, indel_direction='right', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'right_duplication_inversion,{}'.format(event_id)
                    event_id += 1

    ## tandem-dup: mt{k+ k+}, wt{- k+} OR wt{k+ -}
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            block_segs = list(c_mt_block)
            # case: mt{k+ k+}, wt{- k+}, ins seg + concordant block w/ ins-seg as seed
            # this is not the last block, AND next block is concordant
            if c_mt_block_idx + 1 < len(aligned_hap.block_indices) and \
                    (c_mt_block_idx + 1) in aligned_hap.concordant_blocks:
                next_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx + 1]
                seed_start, seed_end = is_seeded(next_concordant_block, block_segs, segment_size_dict, indel_direction='left', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1
                    continue
            # case: mt{k+ k+}, wt{k+ -}, concordant block w/ ins-seg as seed + ins seg
            # this is not the first block, AND previous block is concordant
            if c_mt_block_idx - 1 >= 0 and \
                    (c_mt_block_idx - 1) in aligned_hap.concordant_blocks:
                previous_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx - 1]
                seed_start, seed_end = is_seeded(previous_concordant_block, block_segs, segment_size_dict, indel_direction='right', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1

    ## deletion: all remaining wt_blocks are deletions
    for aligned_hap in aligned_haplotypes:
        for c_wt_block_idx, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_idx]) > 0:
                continue
            else:
                aligned_hap.discordant_block_assignment[c_wt_block_idx] = 'deletion,{}'.format(event_id)
                event_id += 1

    ## unbalanced translocation: all remaining mt_blocks are insertions (i.e. unbalanced translocations)
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            else:
                aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'insertion,{}'.format(event_id)
                event_id += 1

    ## congregate events for report
    event_blocks = {}
    event_types = {}
    for aligned_hap in aligned_haplotypes:
        event_blocks, event_types = aligned_hap.report_SV(event_blocks, event_types)
    sorted_event_id = sorted(list(event_blocks.keys()))

    # if for the same event ID, we have different event types, we have an issue, otherwise, name the event as the singular name
    conglomerated_event_types = {}
    for c_id, event_type_list in event_types.items():
        c_event_type = event_type_list[0]
        for event_type in event_type_list:
            if event_type != c_event_type:
                raise ValueError('same ID, multiple event types')
        conglomerated_event_types[c_id] = c_event_type

    # consolidate neighboring blocks for events with paired blocks
    for event_id, event_type in event_types.items():
        event_type = event_type[0]
        event_block_list = event_blocks[event_id]
        if event_type == 'inversion':
            if len(event_block_list) != 2:
                raise RuntimeError('inversion needs to have 2 consecutive discordant blocks')
            block1_segs = event_block_list[0].split('.')[2].split('(')[1].split(')')[0].split(',')
            block2_segs = event_block_list[1].split('.')[2].split('(')[1].split(')')[0].split(',')
            block1_size = section_size(block1_segs, segment_size_dict)
            block2_size = section_size(block2_segs, segment_size_dict)
            bp1 = event_block_list[0].split('.')[3]
            bp2 = event_block_list[1].split('.')[4]
            if block1_segs[0][-1] == '-':
                # earlier is ins of inverted-segs
                blocks1_uninverted_segs = [seg[:-1] + '+' for seg in block1_segs[::-1]]
                block2_uninverted_segs = block2_segs
            elif block2_segs[0][-1] == '-':
                # later is ins of inverted-segs
                blocks1_uninverted_segs = block1_segs
                block2_uninverted_segs = [seg[:-1] + '+' for seg in block2_segs[::-1]]
            else:
                raise RuntimeError('inversion not found')
            if block1_size >= block2_size:
                c_block_segs = blocks1_uninverted_segs
            else:
                c_block_segs = block2_uninverted_segs
            c_block_str = '({})'.format(','.join(c_block_segs))
            new_block_list = ['{}.{}.{}.{}.{}'.format(event_block_list[0].split('.')[0],
                                                      event_block_list[0].split('.')[1] + ',' + event_block_list[1].split('.')[1],
                                                      c_block_str + ',' + event_block_list[0].split('.')[2] + ',' + event_block_list[1].split('.')[2],
                                                      bp1,
                                                      bp2)]
            event_blocks[event_id] = new_block_list
        elif event_type == 'balanced_translocation':
            if conglomerated_event_types[event_id] != 'balanced_translocation':
                # event already labeled by previous balanced translocation chaining, skip
                continue
            # origin_path_idx = int(event_block_list[0].split('.')[0])
            # origin_block_idx = int(event_block_list[0].split('.')[1])
            c_path_idx = int(event_block_list[1].split('.')[0])
            c_block_idx = int(event_block_list[1].split('.')[1])
            origin_path_idx = c_path_idx
            origin_block_idx = c_block_idx
            event_id_visited = [event_id]
            loop_closed = False
            while True:
                ## find the next paired balanced translocations
                next_block_idx = -1
                next_block_event_id = -1
                # search left
                left_block_idx = c_block_idx - 1
                if left_block_idx in aligned_haplotypes[c_path_idx].discordant_block_assignment:
                    # note: if left_block_idx < 0, won't be in discordant block list
                    left_block_type_info = aligned_haplotypes[c_path_idx].discordant_block_assignment[left_block_idx]
                    left_block_type = left_block_type_info.split(',')[0]
                    if left_block_type == 'balanced_translocation':
                        next_block_idx = left_block_idx
                        next_block_event_id = int(left_block_type_info.split(',')[1])
                right_block_idx = c_block_idx + 1
                if right_block_idx in aligned_haplotypes[c_path_idx].discordant_block_assignment:
                    # note: if left_block_idx >= n_blocks, won't be in discordant block list
                    right_block_type_info = aligned_haplotypes[c_path_idx].discordant_block_assignment[right_block_idx]
                    right_block_type = right_block_type_info.split(',')[0]
                    if right_block_type == 'balanced_translocation':
                        if next_block_idx != -1:
                            raise RuntimeError(
                                'balanced translocation chaining appeared on both sides, require additional implementation for better implementation')
                        next_block_idx = right_block_idx
                        next_block_event_id = int(right_block_type_info.split(',')[1])

                ## navigate to the other side of the pair
                if next_block_idx == -1:
                    # all event_id_visited are non-loop closed balanced translocations
                    break
                # find the located block in event_block pairing
                c_event_block_list = event_blocks[next_block_event_id]
                # either the first/second is the paired block
                path_idx1 = int(c_event_block_list[0].split('.')[0])
                block_idx1 = int(c_event_block_list[0].split('.')[1])
                path_idx2 = int(c_event_block_list[1].split('.')[0])
                block_idx2 = int(c_event_block_list[1].split('.')[1])
                if path_idx1 == c_path_idx and block_idx1 == next_block_idx:
                    if next_block_event_id in event_id_visited and next_block_event_id != event_id_visited[0]:
                        raise RuntimeError('duplicate event found, loop closed not with the origin block')
                    elif next_block_event_id not in event_id_visited:
                        # thus, when we are back to the first event_id, do not append
                        event_id_visited.append(next_block_event_id)
                    c_path_idx = path_idx2
                    c_block_idx = block_idx2
                elif path_idx2 == c_path_idx and block_idx2 == next_block_idx:
                    if next_block_event_id in event_id_visited and next_block_event_id != event_id_visited[0]:
                        raise RuntimeError('duplicate event found, loop closed not with the origin block')
                    elif next_block_event_id not in event_id_visited:
                        # thus, when we are back to the first event_id, do not append
                        event_id_visited.append(next_block_event_id)
                    c_path_idx = path_idx1
                    c_block_idx = block_idx1
                else:
                    raise RuntimeError('bug in event finding, block pair located was not found')

                ## base-case: returned to the original start, loop closed
                if c_path_idx == origin_path_idx and c_block_idx == origin_block_idx:
                    loop_closed = True
                    break

            if not loop_closed:
                for event_id_itr in event_id_visited:
                    if conglomerated_event_types[event_id_itr] != 'balanced_translocation':
                        raise RuntimeError('illegal labeling of balanced translocation')
                    conglomerated_event_types[event_id_itr] = 'balanced_translocation_unassociated'
            else:
                for event_id_idx, event_id_itr in enumerate(event_id_visited):
                    if event_id_idx >= len(event_id_visited) - 1:
                        next_event_id_idx = 0
                    else:
                        next_event_id_idx = event_id_idx + 1
                    if conglomerated_event_types[event_id_itr] != 'balanced_translocation':
                        raise RuntimeError('illegal labeling of balanced translocation')
                    conglomerated_event_types[event_id_itr] = 'balanced_translocation_associated' + '<{}>'.format(event_id_visited[next_event_id_idx])

        elif event_type.startswith('balanced_translocation_associated') or \
                event_type.startswith('balanced_translocation_unassociated'):
            # already labeled
            continue

    output_list = []
    for event_id in sorted_event_id:
        print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
        output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    return output_list, aligned_haplotypes


def continuous_extension(input_hap, idx_ptr):
    """
    start from the idx_ptr location, moving leftward until no longer can form a continuous section
    :param idx_ptr:
    :param input_hap:
    :return: final_ptr, where [idx_ptr, final_ptr) is a continuous section
    """
    final_ptr = idx_ptr + 1
    hap_len = len(input_hap)
    if idx_ptr >= hap_len:
        raise ValueError('idx_ptr out of bound')

    current_seg = input_hap[idx_ptr]
    section_sign = current_seg[-1]
    while True:
        if final_ptr == hap_len:
            break

        next_seg = input_hap[final_ptr]
        if next_seg == '-':
            break
        if section_sign == '+':
            if next_seg[-1] != '+':
                break
            elif int(next_seg[:-1]) != int(current_seg[:-1]) + 1:
                break
        elif section_sign == '-':
            if next_seg[-1] != '-':
                break
            elif int(next_seg[:-1]) != int(current_seg[:-1]) - 1:
                break
        else:
            raise ValueError('sign error')
        current_seg = input_hap[final_ptr]
        final_ptr += 1

    return final_ptr


def is_seeded(supergroup_section, cont_section, size_dict, indel_direction, d, eps):
    """
    whether cont_section is seeded in input_hap, with significant size (>d)
    :param supergroup_section:
    :param cont_section:
    :param size_dict:
    :param indel_direction: 'both', 'left', OR 'right': which direction to compute indels size on
    :param d: required seed overlap size
    :param eps: epsilon, allowed 1/2 indel size (i.e. 2 * eps >= indel size)
    :return: [start, end) indices in the supergroup
    """

    def sublist_idx(list1, list2):
        # list 1 is small, list 2 is large (presumed superlist)
        for i in range(len(list2) - len(list1) + 1):
            if list2[i:i + len(list1)] == list1:
                return i
        return None

    all_sublists = [cont_section[i:j + 1] for i in range(len(cont_section)) for j in range(i, len(cont_section))]
    sublist_sizes = [section_size(sublist, size_dict) for sublist in all_sublists]

    sublists_found = []
    seed_start_location_in_supergroup = []
    max_size = -1
    max_size_sublist_idx = -1
    for idx, sublist in enumerate(all_sublists):
        current_size = sublist_sizes[idx]
        if current_size < d:
            seed_start_location_in_supergroup.append(-1)
            continue
        current_start_location = sublist_idx(sublist, supergroup_section)
        if current_start_location is not None:
            sublists_found.append(sublist)
            seed_start_location_in_supergroup.append(current_start_location)
            if current_size > max_size:
                max_size = current_size
                max_size_sublist_idx = idx
        else:
            seed_start_location_in_supergroup.append(-1)

    if max_size_sublist_idx == -1:
        # not found
        return -1, -1

    max_size_sublist = all_sublists[max_size_sublist_idx]
    max_size_sublist_start_location = seed_start_location_in_supergroup[max_size_sublist_idx]
    max_size_sublist_end_location = seed_start_location_in_supergroup[max_size_sublist_idx] + len(max_size_sublist)

    del_size = section_size(cont_section, size_dict) - section_size(max_size_sublist, size_dict)
    left_ins_size = section_size(supergroup_section[:max_size_sublist_start_location], size_dict)
    right_ins_size = section_size(supergroup_section[max_size_sublist_end_location:], size_dict)

    if indel_direction == 'both':
        indel_size = del_size + left_ins_size + right_ins_size
        if indel_size <= 2 * eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    elif indel_direction == 'left':
        indel_size = del_size + left_ins_size
        if indel_size <= eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    elif indel_direction == 'right':
        indel_size = del_size + right_ins_size
        if indel_size <= eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    else:
        raise ValueError('indel_direction parameter input illegal')

    # indel size too large
    return -1, -1


def section_size(input_section, size_dict):
    tot_size = 0
    for seg in input_section:
        tot_size += size_dict[seg[:-1]]
    return tot_size


def populate_wt_indexed_lists(mt_path_chrs, wt_path_dict):
    wt_indexed_lists = []
    for path_chr in mt_path_chrs:
        wt_indexed_lists.append(wt_path_dict[path_chr])
    return wt_indexed_lists


def form_dependent_clusters(input_events, i_aligned_haplotypes):
    """
    :param input_events:
    :return: list of list, inner-list contains the path_id of the cluster, outer-list sorted by min-value path_id in the list
    """
    dependent_clusters = []
    cluster_events = []
    cluster_chrs = []
    for event_info in input_events:
        event_paths = set()
        for block in event_info[2]:
            event_paths.add(int(block.split('.')[0]))
        event_paths = list(event_paths)
        # see if there exists a cluster to add
        cluster_found = -1
        for cluster_idx, cluster in enumerate(dependent_clusters):
            # search in all clusters
            for event_path in event_paths:
                # if there exists one path overlap
                path_chr = i_aligned_haplotypes[event_path].chrom
                if path_chr in cluster_chrs[cluster_idx]:
                    cluster_found = cluster_idx
                    break

        if cluster_found >= 0:
            for event_path in event_paths:
                located_cluster = dependent_clusters[cluster_found]
                if event_path not in located_cluster:
                    located_cluster.append(event_path)
                    cluster_chrs[cluster_found].add(i_aligned_haplotypes[event_path].chrom)
            cluster_events[cluster_found].append(event_info)
        elif cluster_found == -1:
            dependent_clusters.append(event_paths)
            cluster_events.append([event_info])
            c_path_chr = set()
            for event_path in event_paths:
                c_path_chr.add(i_aligned_haplotypes[event_path].chrom)
            cluster_chrs.append(c_path_chr)
        else:
            raise RuntimeError


    ## sort clusters by min_element
    min_values = [min(inner_list) for inner_list in dependent_clusters]
    permutations_indices = sorted(range(len(min_values)), key=lambda k: min_values[k])
    sorted_dependent_clusters = [dependent_clusters[i] for i in permutations_indices]
    sorted_cluster_events = [cluster_events[i] for i in permutations_indices]

    return sorted_dependent_clusters, sorted_cluster_events


########################PERIPHERAL UTILS#######################
def gather_breakpoints(breakpoints: []):
    all_breakpoints = []
    for c_breakpoint in breakpoints:
        if c_breakpoint[1] is None:
            continue
        all_breakpoints.append((c_breakpoint[0], c_breakpoint[1]))
    return all_breakpoints


def get_genes_near_breakpoints(breakpoints: [(str, int)], proximity=50000):
    breakpoint_ranges = []
    for c_breakpoint in breakpoints:
        breakpoint_ranges.append((c_breakpoint[0],
                                  c_breakpoint[1] - proximity,
                                  c_breakpoint[1] + proximity))
    genes_in_regions = set()
    for breakpoint_range in breakpoint_ranges:
        genes_in_regions = genes_in_regions.union(get_genes_in_region(*breakpoint_range))
    return list(genes_in_regions)


def report_on_genes_based_on_breakpoints(breakpoints):
    breakpoints = gather_breakpoints(breakpoints)
    genes_near_bp = get_genes_near_breakpoints(breakpoints)
    DDG_df = get_DDG_overlapped_genes(genes_near_bp)
    DDG_gene_list, DDG_disease_list = tostring_gene_disease_omim(DDG_df)
    return genes_near_bp, list(zip(DDG_gene_list, DDG_disease_list))


def report_cnv_genes_on_region(chrom, start, end, proximity=50000):
    start = max(0, start - proximity)
    end = end + proximity
    genes = get_genes_in_region('chr' + chrom, start, end)
    DDG_df = get_DDG_overlapped_genes(genes)
    DDG_gene_list, DDG_disease_list = tostring_gene_disease_omim(DDG_df)
    return genes, list(zip(DDG_gene_list, DDG_disease_list))


def get_ucsc_url(chrom, start_pos, end_pos, db='hg38'):
    prefix = 'https://genome.ucsc.edu/cgi-bin/hgTracks?db={}' \
             '&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position='.format(db)
    suffix = '{}%3A{}%2D{}'.format(chrom.lower(), start_pos, end_pos)
    return prefix + suffix


def sort_events(input_events):
    def sort_key(event):
        event_type = event[1]
        if 'translocation' in event_type:
            return -1.0
        else:
            path_cluster_idx = '.'.join(event[2][0].split('.')[:2]).split(',')[0]
            return float(path_cluster_idx)

    return sorted(input_events, key=sort_key)


def test_interpreter():
    i_mt_hap1 = ['1+', '1+', '2+', '3-', '4+', '9+', '10+']
    i_mt_hap2 = ['7+', '8+', '5+', '6+']
    i_mt_hap3 = ['11+', '12+', '11-']
    i_mt_hap4 = ['14-', '13+', '14+']
    i_mt_hap5 = ['15+', '16+', '17-', '16-', '18+']
    i_wt_hap1 = ['1+', '2+', '3+', '4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+', '10+']
    i_wt_hap3 = ['11+', '12+']
    i_wt_hap4 = ['13+', '14+']
    i_wt_hap5 = ['15+', '16+', '17+', '18+']
    i_mt_list = [i_mt_hap1, i_mt_hap2, i_mt_hap3, i_mt_hap4, i_mt_hap5]
    i_wt_list = [i_wt_hap1, i_wt_hap2, i_wt_hap3, i_wt_hap4, i_wt_hap5]
    i_size_dict = {str(i): 1 for i in range(19)}
    # print(lcs(i_mt_hap, i_wt_hap, i_size_dict))
    out = interpret_haplotypes(i_mt_list, i_wt_list, i_size_dict, 1, 1)
    print(out)


def test_3break_qter():
    i_wt_hap0 = ['1+', '2+', '3+']
    i_wt_hap1 = ['4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+']
    i_mt_hap0 = ['1+', '2+', '9+']
    i_mt_hap1 = ['4+', '5+', '3+']
    i_mt_hap2 = ['7+', '8+', '6+']
    chrom_id = ['Chr1', 'Chr2', 'Chr3']
    i_mt_list = [i_mt_hap0, i_mt_hap1, i_mt_hap2]
    i_wt_list = [i_wt_hap0, i_wt_hap1, i_wt_hap2]
    i_size_dict = {str(i): 1 for i in range(19)}
    out = interpret_haplotypes(i_mt_list, i_wt_list, chrom_id, i_size_dict, 1, 1)
    print(out)


def test_reciprocal_trans():
    i_wt_hap0 = ['1+', '2+', '3+']
    i_wt_hap1 = ['4+', '5+', '6+']
    i_mt_hap0 = ['1+', '2+', '6+']
    i_mt_hap1 = ['3+', '4+', '5+']
    chrom_id = ['Chr1', 'Chr2']
    i_mt_list = [i_mt_hap0, i_mt_hap1]
    i_wt_list = [i_wt_hap0, i_wt_hap1]
    i_size_dict = {str(i): 1 for i in range(19)}
    out = interpret_haplotypes(i_mt_list, i_wt_list, chrom_id, i_size_dict, 1, 1)
    print(out)


def test_ucsc_url():
    x = get_ucsc_url('Chr2', 130271298, 131386307)
    print(x)


# def test_segs_union():
#     l1 = ['7+', '8+']
#     l2 = ['7+']
#     l3 = segs_union(l1, l2)
#     print(l3)

if __name__ == "__main__":
    test_ucsc_url()
