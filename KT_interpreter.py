class Aligned_Haplotype:
    id: int  # for debugging
    alignment_score: int
    mt_aligned: [str]
    wt_aligned: [str]
    block_indices: [(int, int)]  # [(a, b)] where event block are [a, b)
    mt_blocks: {(str,): [int]}  # block segments: block index (n-th block) in block_indices
    wt_blocks: {(str,): [int]}
    mt_blocks_inv: {int: (str,)}
    wt_blocks_inv: {int: (str,)}
    block_assignment: {int: str}  # block index: event name, exactly one event can be assigned to one block
    size_dict: {str: int}  # seg_name: bp size

    def __init__(self, mt_aligned, wt_aligned, alignment_score, size_dict, id):
        self.id = id
        self.mt_aligned = mt_aligned
        self.wt_aligned = wt_aligned
        self.alignment_score = alignment_score
        self.size_dict = size_dict

        ## get event_blocks
        event_block = []
        block_types = []  # [str], a block is either ins/del
        # we know len(wt_hap) == len(mt_hap) as they are alignments
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
        self.block_indices = event_block

        ## extract mt/wt blocks
        # since blocks are cont., it is either a whole section of INS or a whole DEL
        self.wt_blocks = {}
        self.mt_blocks = {}
        self.block_assignment = {}
        for block_idx, block_info in enumerate(self.block_indices):
            self.block_assignment[block_idx] = ""
            c_block_start = block_info[0]
            if self.mt_aligned[c_block_start] == '-':
                block_segs = self.wt_aligned[block_info[0]: block_info[1]]
                if tuple(block_segs) in self.wt_blocks:
                    self.wt_blocks[tuple(block_segs)].append(block_idx)
                else:
                    self.wt_blocks[tuple(block_segs)] = [block_idx]
            elif self.wt_aligned[c_block_start] == '-':
                block_segs = self.mt_aligned[block_info[0]: block_info[1]]
                if tuple(block_segs) in self.mt_blocks:
                    self.mt_blocks[tuple(block_segs)].append(block_idx)
                else:
                    self.mt_blocks[tuple(block_segs)] = [block_idx]
            else:
                raise ValueError('mismatch not allowed')

        self.wt_blocks_inv = {}
        self.mt_blocks_inv = {}
        for block, block_idx_list in self.wt_blocks.items():
            for block_idx in block_idx_list:
                self.wt_blocks_inv[block_idx] = block
        for block, block_idx_list in self.mt_blocks.items():
            for block_idx in block_idx_list:
                self.mt_blocks_inv[block_idx] = block

    def __str__(self):
        return "ID<{}>".format(self.id)

    def report_SV(self, event_blocks, event_types):
        """
        :param event_blocks: cumulative dict, where value is a list of event_block, this helps to group paired-events that were on different chromosomes
        :param event_types: cumulative dict, where value is a list of event_type, same reason as above
        :return:
        """
        for mt_block, mt_block_idx_list in self.mt_blocks.items():
            mt_block_str = ','.join(list(mt_block))
            for mt_block_idx in mt_block_idx_list:
                event_id = int(self.block_assignment[mt_block_idx].split(',')[-1])
                event_type = self.block_assignment[mt_block_idx].split(',')[0]
                if event_id in event_blocks:
                    event_blocks[event_id].append('{}.mt({})'.format(self.id, mt_block_str))
                    event_types[event_id].append(event_type)
                else:
                    event_blocks[event_id] = ['{}.mt({})'.format(self.id, mt_block_str)]
                    event_types[event_id] = [event_type]
        for wt_block, wt_block_idx_list in self.wt_blocks.items():
            wt_block_str = ','.join(list(wt_block))
            for wt_block_idx in wt_block_idx_list:
                event_id = int(self.block_assignment[wt_block_idx].split(',')[-1])
                event_type = self.block_assignment[wt_block_idx].split(',')[0]
                if event_id in event_blocks:
                    event_blocks[event_id].append('{}.wt({})'.format(self.id, wt_block_str))
                    event_types[event_id].append(event_type)
                else:
                    event_blocks[event_id] = ['{}.wt({})'.format(self.id, wt_block_str)]
                    event_types[event_id] = [event_type]
        return event_blocks, event_types

    def search_seed(self, query_section, seed_type, d=200000, eps=200000):
        """
        :param query_section:
        :param seed_type:
        :param d:
        :param eps:
        :return: block_idx if found, -1 if not
        """
        # search in the opposite type
        if seed_type == 'ins':
            query_block = self.wt_blocks
        elif seed_type == 'del':
            query_block = self.mt_blocks
        else:
            raise ValueError('invalid seed type')
        for c_block, c_block_idx_list in query_block.items():
            for c_block_idx in c_block_idx_list:
                if len(self.block_assignment[c_block_idx]) != 0:
                    # already has assignment
                    # TODO: can improve this by bipartite matching, instead of greedy matching
                    continue
                if is_seeded(c_block, query_section, self.size_dict, d, eps):
                    return c_block_idx
        return -1


def interpret_haplotypes(mt_hap_list: [[str]], wt_hap_list: [[str]], segment_size_dict: {str: int}):
    """
    Input a haplotype and a WT, report SVs; hap from mt list must correspond to hap from wt list
    :param segment_size_dict: {str: int} segment mapped to size of the segment (eg. 11: 12,345)
    :param mt_hap_list:
    :param wt_hap_list:
    :return:
    """
    event_id = 1  # static variable, keep different events separate, keep paired events together (eg. balanced trans)
    hap_id = 0  # used for naming haplotypes, preserves the order they come-in in mt_hap/wt_hap lists
    aligned_haplotypes = []
    for idx, mt_hap in enumerate(mt_hap_list):
        wt_hap = wt_hap_list[idx]
        a1, a2, a3 = lcs(mt_hap, wt_hap, segment_size_dict)
        aligned_haplotypes.append(Aligned_Haplotype(a2, a3, a1, segment_size_dict, hap_id))
        hap_id += 1

    def genomewide_seed_search(query_section, query_block_idx, query_hap_idx, query_type):
        seed_found = False
        for aligned_hap_idx1, aligned_hap1 in enumerate(aligned_haplotypes):
            if aligned_hap_idx1 == query_hap_idx:
                # first get inter-chr event
                continue
            # TODO: reset d,eps from debug mode
            seeded_block_idx = aligned_hap1.search_seed(query_section, query_type, d=2, eps=1)
            if seeded_block_idx != -1:
                aligned_hap1.block_assignment[seeded_block_idx] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].block_assignment[query_block_idx] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
                break
        if not seed_found:
            # at last, check intra-chr event
            seeded_block_idx = aligned_haplotypes[query_hap_idx].search_seed(query_section, query_type, d=2, eps=1)
            if seeded_block_idx != -1:
                aligned_haplotypes[query_hap_idx].block_assignment[seeded_block_idx] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].block_assignment[query_block_idx] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
        if seed_found:
            return True
        else:
            return False

    ### Order of resolving: balanced-translocation, inv, dup-inv, tandem-dup, del, unbalanced-translocation

    ## resolve all balanced translocations
    for aligned_hap_idx, aligned_hap in enumerate(aligned_haplotypes):
        for c_wt_block, c_wt_block_idx_list in aligned_hap.wt_blocks.items():
            for c_wt_block_idx in c_wt_block_idx_list:
                if len(aligned_hap.block_assignment[c_wt_block_idx]) > 0:
                    # already has assignment
                    continue
                aligned_hap.block_assignment[c_wt_block_idx] = 'under investigation'  # block it from being matched in intra-chr event
                balanced_translocation_found = genomewide_seed_search(c_wt_block, c_wt_block_idx, aligned_hap_idx, 'del')
                if balanced_translocation_found:
                    event_id += 1
                else:
                    # reset block assignment for query
                    aligned_hap.block_assignment[c_wt_block_idx] = ''
        for c_mt_block, c_mt_block_idx_list in aligned_hap.mt_blocks.items():
            for c_mt_block_idx in c_mt_block_idx_list:
                if len(aligned_hap.block_assignment[c_mt_block_idx]) > 0:
                    # already has assignment
                    continue
                aligned_hap.block_assignment[c_mt_block_idx] = 'under investigation'  # block it from being matched in intra-chr event
                balanced_translocation_found = genomewide_seed_search(c_mt_block, c_mt_block_idx, aligned_hap_idx, 'ins')
                if balanced_translocation_found:
                    event_id += 1
                else:
                    # reset block assignment for query
                    aligned_hap.block_assignment[c_mt_block_idx] = ''

    ## inversion: mt{- k-} wt{k+ -} OR mt{k- -} wt{- k+}
    for aligned_hap in aligned_haplotypes:
        # case: mt{k- -} wt{- k+}
        for c_mt_block, c_mt_block_idx_list in aligned_hap.mt_blocks.items():
            for c_mt_block_idx in c_mt_block_idx_list:
                if c_mt_block[0][-1] == "+":
                    # not inverted
                    continue
                if len(aligned_hap.block_assignment[c_mt_block_idx]) > 0:
                    continue
                if c_mt_block_idx + 1 >= len(aligned_hap.block_assignment) or len(aligned_hap.block_assignment[c_mt_block_idx + 1]) > 0:
                    # next block unavailable
                    continue
                uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
                if c_mt_block_idx + 1 in aligned_hap.wt_blocks_inv and list(aligned_hap.wt_blocks_inv[c_mt_block_idx + 1]) == uninverted_segs:
                    aligned_hap.block_assignment[c_mt_block_idx] = 'inversion,{}'.format(event_id)
                    aligned_hap.block_assignment[c_mt_block_idx + 1] = 'inversion,{}'.format(event_id)
                    event_id += 1
        # case: mt{- k-} wt{k+ -}
        for c_wt_block, c_wt_block_idx_list in aligned_hap.wt_blocks.items():
            for c_wt_block_idx in c_wt_block_idx_list:
                if len(aligned_hap.block_assignment[c_wt_block_idx]) > 0:
                    continue
                if c_wt_block_idx + 1 >= len(aligned_hap.block_assignment) or len(aligned_hap.block_assignment[c_wt_block_idx + 1]) > 0:
                    # next block unavailable
                    continue
                inverted_segs = [seg[:-1] + '-' for seg in c_wt_block[::-1]]
                if c_wt_block_idx + 1 in aligned_hap.mt_blocks_inv and list(aligned_hap.mt_blocks_inv[c_wt_block_idx + 1]) == inverted_segs:
                    aligned_hap.block_assignment[c_wt_block_idx] = 'inversion,{}'.format(event_id)
                    aligned_hap.block_assignment[c_wt_block_idx + 1] = 'inversion,{}'.format(event_id)
                    event_id += 1

    ## duplication inversion
    for aligned_hap in aligned_haplotypes:
        for c_mt_block, c_mt_block_idx_list in aligned_hap.mt_blocks.items():
            for c_mt_block_idx in c_mt_block_idx_list:
                if c_mt_block[0][-1] == "+":
                    # not inverted
                    continue
                if len(aligned_hap.block_assignment[c_mt_block_idx]) > 0:
                    continue
                block_len = len(c_mt_block)
                uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
                # left-dup-inv
                nonblock_start = aligned_hap.block_indices[c_mt_block_idx][1]
                nonblock_end = nonblock_start + block_len
                if nonblock_end > len(aligned_hap.wt_aligned):
                    continue
                if aligned_hap.wt_aligned[nonblock_start: nonblock_end] == uninverted_segs and aligned_hap.mt_aligned[nonblock_start: nonblock_end]== uninverted_segs:
                    aligned_hap.block_assignment[c_mt_block_idx] = 'left_duplication_inversion,{}'.format(event_id)
                    event_id += 1
                    continue
                # right-dup-inv
                nonblock_end = aligned_hap.block_indices[c_mt_block_idx][0]
                nonblock_start = nonblock_end - block_len
                if nonblock_start < 0:
                    continue
                if aligned_hap.wt_aligned[nonblock_start: nonblock_end] == uninverted_segs and aligned_hap.mt_aligned[nonblock_start: nonblock_end]== uninverted_segs:
                    aligned_hap.block_assignment[c_mt_block_idx] = 'right_duplication_inversion,{}'.format(event_id)
                    event_id += 1

    ## tandem-dup: mt{k+ k+}, wt{- k+} OR wt{k+ -}
    for aligned_hap in aligned_haplotypes:
        for c_mt_block, c_mt_block_idx_list in aligned_hap.mt_blocks.items():
            for c_mt_block_idx in c_mt_block_idx_list:
                if len(aligned_hap.block_assignment[c_mt_block_idx]) > 0:
                    continue
                block_len = len(c_mt_block)
                block_segs = list(c_mt_block)
                # case: mt{k+ k+}, wt{- k+}
                nonblock_start = aligned_hap.block_indices[c_mt_block_idx][1]
                nonblock_end = nonblock_start + block_len
                if nonblock_end > len(aligned_hap.wt_aligned):
                    continue
                if aligned_hap.wt_aligned[nonblock_start: nonblock_end] == block_segs and aligned_hap.mt_aligned[nonblock_start: nonblock_end] == block_segs:
                    aligned_hap.block_assignment[c_mt_block_idx] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1
                    continue
                # case: mt{k+ k+}, wt{k+ -}
                nonblock_end = aligned_hap.block_indices[c_mt_block_idx][0]
                nonblock_start = nonblock_end - block_len
                if nonblock_start < 0:
                    continue
                if aligned_hap.wt_aligned[nonblock_start: nonblock_end] == block_segs and aligned_hap.mt_aligned[nonblock_start: nonblock_end]== block_segs:
                    aligned_hap.block_assignment[c_mt_block_idx] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1

    ## deletion: all remaining wt_blocks are deletions
    for aligned_hap in aligned_haplotypes:
        for c_wt_block, c_wt_block_idx_list in aligned_hap.wt_blocks.items():
            for c_wt_block_idx in c_wt_block_idx_list:
                if len(aligned_hap.block_assignment[c_wt_block_idx]) > 0:
                    continue
                else:
                    aligned_hap.block_assignment[c_wt_block_idx] = 'deletion,{}'.format(event_id)
                    event_id += 1

    ## unbalanced translocation: all remaining mt_blocks are insertions (i.e. unbalanced translocations)
    for aligned_hap in aligned_haplotypes:
        for c_mt_block, c_mt_block_idx_list in aligned_hap.mt_blocks.items():
            for c_mt_block_idx in c_mt_block_idx_list:
                if len(aligned_hap.block_assignment[c_mt_block_idx]) > 0:
                    continue
                else:
                    aligned_hap.block_assignment[c_mt_block_idx] = 'unbalanced_translocation,{}'.format(event_id)
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

    output_list = []
    for event_id in sorted_event_id:
        print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
        output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    return output_list


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
            if int(next_seg[:-1]) != int(current_seg[:-1]) + 1:
                break
        elif section_sign == '-':
            if int(next_seg[:-1]) != int(current_seg[:-1]) - 1:
                break
        else:
            raise ValueError('sign error')
        final_ptr += 1

    return final_ptr


def is_seeded(supergroup_section, cont_section, size_dict, d=200000, eps=200000.0):
    """
    whether cont_section is seeded in input_hap, with significant size (>d)
    :param supergroup_section:
    :param cont_section:
    :param size_dict:
    :param d: required seed overlap size
    :param eps: epsilon, allowed 1/2 indel size
    :return:
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
    max_size = -1
    max_size_sublist_idx = -1
    for idx, sublist in enumerate(all_sublists):
        current_size = sublist_sizes[idx]
        if current_size < d:
            continue
        if sublist_idx(sublist, supergroup_section) is not None:
            sublists_found.append(sublist)
            if current_size > max_size:
                max_size = current_size
                max_size_sublist_idx = idx

    if max_size_sublist_idx == -1:
        # not found
        return False

    max_size_sublist = all_sublists[max_size_sublist_idx]
    del_size = section_size(cont_section, size_dict) - section_size(max_size_sublist, size_dict)
    ins_size = section_size(supergroup_section, size_dict) - section_size(max_size_sublist, size_dict)

    if del_size + ins_size <= 2 * eps:
        return True
    else:
        # introduced too much indel
        return False


def section_size(input_section, size_dict):
    tot_size = 0
    for seg in input_section:
        tot_size += size_dict[seg[:-1]]
    return tot_size


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


if __name__ == "__main__":
    i_mt_hap1 = ['1+', '1+', '2+', '3-', '4+', '9+', '10+']
    i_mt_hap2 = ['7+', '8+', '5+', '6+']
    i_mt_hap3 = ['11+', '11-', '12+']
    i_mt_hap4 = ['13+', '14-', '14+']
    i_wt_hap1 = ['1+', '2+', '3+', '4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+', '10+']
    i_wt_hap3 = ['11+', '12+']
    i_wt_hap4 = ['13+', '14+']
    i_mt_list = [i_mt_hap1, i_mt_hap2, i_mt_hap3, i_mt_hap4]
    i_wt_list = [i_wt_hap1, i_wt_hap2, i_wt_hap3, i_wt_hap4]
    i_size_dict = {str(i): 1 for i in range(15)}
    # print(lcs(i_mt_hap, i_wt_hap, i_size_dict))
    out = interpret_haplotypes(i_mt_list, i_wt_list, i_size_dict)
    print(out)
