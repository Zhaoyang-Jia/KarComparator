def interpret_haplotypes(mt_hap_list: [[str]], wt_hap_list: [[str]], segment_size_dict: {str: int}):
    """
    Input a haplotype and a WT, report SVs; hap from mt list must correspond to hap from wt list
    :param segment_size_dict: {str: int} segment mapped to size of the segment (eg. 11: 12,345)
    :param mt_hap_list:
    :param wt_hap_list:
    :return:
    """
    mt_aligned = []
    wt_aligned = []
    for idx, mt_hap in enumerate(mt_hap_list):
        wt_hap = wt_hap_list[idx]
        _, a1, a2 = lcs(mt_hap, wt_hap, segment_size_dict)
        mt_aligned.append(a1)
        wt_aligned.append(a2)

    event_blocks = []  # [[(int, int)]]; [[(a, b)]] where event block are [a, b) of each hap
    block_types = []  # [[str]], a block is either ins/del
    for hap_idx1, mt_hap in enumerate(mt_aligned):
        # we know len(wt_hap) == len(mt_hap) as they are alignments
        wt_hap = wt_aligned[hap_idx1]
        hap_len = len(mt_hap)
        seg_idx = 0
        hap_blocks = []
        while seg_idx < hap_len:
            mt_seg = mt_hap[seg_idx]
            wt_seg = wt_hap[seg_idx]
            if mt_seg == wt_seg:
                seg_idx += 1
                continue
            else:
                # this will be the first seg of (potentially) a section, that is indel
                if wt_seg == '-':
                    # insertion
                    continuous_end = continuous_extension(mt_hap, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if wt_hap[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('ins')
                elif mt_seg == '-':
                    # deletion
                    continuous_end = continuous_extension(wt_hap, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if mt_hap[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('del')
                else:
                    raise ValueError('mismatch not allowed')
                hap_blocks.append((seg_idx, block_end))
                seg_idx = block_end
        event_blocks.append(hap_blocks)

    def genomewide_seed_search(query_section, query_chr_idx, query_type):
        for hap_idx2, hap_block_list2 in enumerate(event_blocks):
            if hap_idx2 == hap_idx1:
                # first get inter-chr event
                continue
            for block_idx2, current_block2 in enumerate(hap_block_list2):
                if block_types[block_idx2] == query_type:
                    # we only search seed in the opposite type
                    continue
    #xxx
    for hap_idx1, hap_block_list1 in enumerate(event_blocks):
        for block_idx1, current_block1 in enumerate(hap_block_list1):
            current_type = block_types[block_idx1]
            block_start = current_block1[0]
            block_end = current_block1[1]
            if current_type == 'ins':
                cont_section = mt_aligned[hap_idx1][block_start: block_end]
                # ins search for del seeds in all other haps, then its own hap



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


def is_seeded(input_hap, cont_section, size_dict, d=200000):
    """
    whether cont_section is seeded in input_hap, with significant size (>d)
    :param input_hap:
    :param cont_section:
    :param size_dict:
    :param d:
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
        if sublist_idx(sublist, input_hap) is not None:
            sublists_found.append(sublist)
            if current_size > max_size:
                max_size = current_size
                max_size_sublist_idx = idx

    if max_size_sublist_idx != -1:
        return sublist_idx(all_sublists[max_size_sublist_idx], input_hap)
    else:
        return 0


def section_size(input_section, size_dict):
    tot_size = 0
    for seg in input_section:
        tot_size += size_dict[seg]
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
    i_mt_hap1 = ['1+', '2+', '3-', '4+', '9+', '10+']
    i_mt_hap2 = ['7+', '8+', '5+', '6+']
    i_wt_hap1 = ['1+', '2+', '3+', '4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+', '10+']
    i_mt_list = [i_mt_hap1, i_mt_hap2]
    i_wt_list = [i_wt_hap1, i_wt_hap2]
    i_size_dict = {'1': 1,
                   '2': 1,
                   '3': 1,
                   '4': 1,
                   '5': 1,
                   '6': 1,
                   '7': 1,
                   '8': 1,
                   '9': 1,
                   '10': 1}
    # print(lcs(i_mt_hap, i_wt_hap, i_size_dict))
    print(interpret_haplotypes(i_mt_list, i_wt_list, i_size_dict))
