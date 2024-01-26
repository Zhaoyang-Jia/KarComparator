from Structures import *


def align_paths(segment_list1, segment_list2):
    forbidden_comparison_region_types = ['acrocentric', 'telomere1', 'telomere2', 'acrocentric-telomere1', 'acrocentric-centromere']
    indel_penalty_per_nt = -1
    alignment_1 = []
    alignment_2 = []
    scoring_matrix = [[0 for i in range(0, len(segment_list2) + 1)] for j in
                      range(0, len(segment_list1) + 1)]
    backtrack_matrix = [["" for i in range(0, len(segment_list2) + 1)] for j in
                        range(0, len(segment_list1) + 1)]

    # initialize starting grid
    scoring_matrix[0][0] = 0
    for row_index in range(1, len(segment_list1) + 1):
        current_segment = segment_list1[row_index - 1]
        if current_segment.segment_type not in forbidden_comparison_region_types:
            scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0] + len(current_segment) * indel_penalty_per_nt
        else:
            # no penalty if gaped on forbidden regions
            scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0]
        backtrack_matrix[row_index][0] = "down"

    for col_index in range(1, len(segment_list2) + 1):
        current_segment = segment_list2[col_index - 1]
        if current_segment.segment_type not in forbidden_comparison_region_types:
            scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1] + len(current_segment) * indel_penalty_per_nt
        else:
            scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1]
        backtrack_matrix[0][col_index] = "rigt"

    # DP recursion
    for row_index in range(1, len(segment_list1) + 1):
        for col_index in range(1, len(segment_list2) + 1):
            if segment_list1[row_index - 1].segment_type in forbidden_comparison_region_types:
                # forbidden region indel
                down_value = scoring_matrix[row_index - 1][col_index]
            else:
                # standard indel
                down_value = scoring_matrix[row_index - 1][col_index] \
                             + len(segment_list1[row_index - 1]) * indel_penalty_per_nt

            if segment_list2[col_index - 1].segment_type in forbidden_comparison_region_types:
                # forbidden region indel
                right_value = scoring_matrix[row_index][col_index - 1]
            else:
                # standard indel
                right_value = scoring_matrix[row_index][col_index - 1] \
                              + len(segment_list2[col_index - 1]) * indel_penalty_per_nt

            if segment_list1[row_index - 1] != segment_list2[col_index - 1]:
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
    final_score = scoring_matrix[len(segment_list1)][len(segment_list2)]
    current_row = len(segment_list1)
    current_col = len(segment_list2)

    while True:
        if current_row == 0 and current_col == 0:
            break
        if backtrack_matrix[current_row][current_col] == "diag":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, segment_list2[current_col - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, segment_list1[current_row - 1])
            alignment_2.insert(0, "-")
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "rigt":
            alignment_1.insert(0, "-")
            alignment_2.insert(0, segment_list2[current_col - 1])
            current_col -= 1
        else:
            raise RuntimeError("error in backtrack matrix")

    return final_score, alignment_1, alignment_2


def tostring_alignment(index_to_segment, alignment_1, alignment_2):
    """
    :param index_to_segment: segment dict, key=int, value=Segment
    :param alignment_1: KarSim
    :param alignment_2: OMKar
    :return:
    """
    forbidden_comparison_region_types = ['acrocentric',
                                         'telomere1',
                                         'telomere2',
                                         'acrocentric-telomere1',
                                         'acrocentric-centromere']

    segment_to_index = flip_dict(index_to_segment)
    str1 = []
    str2 = []

    for segment in alignment_1:
        if segment == "-":
            str1.append("-")
        elif segment == 's':
            str1.append('skipped')
        else:
            appending_str = ""
            # to denote no penalty for certain regions
            if any(substring in segment.segment_type for substring in forbidden_comparison_region_types):
                appending_str += "x"
            # get segment index
            if segment not in segment_to_index:
                new_segment = segment.duplicate()
                new_segment.invert()
                appending_str += str(segment_to_index[new_segment]) + "-"
            else:
                appending_str += str(segment_to_index[segment]) + "+"
            str1.append(appending_str)

    for segment in alignment_2:
        if segment == "-":
            str2.append("-")
        elif segment == 's':
            str2.append('skipped')
        else:
            appending_str = ""
            # to denote no penalty for certain regions
            if any(substring in segment.segment_type for substring in forbidden_comparison_region_types):
                appending_str += "x"
            # get segment index
            if segment not in segment_to_index:
                new_segment = segment.duplicate()
                new_segment.invert()
                appending_str += str(segment_to_index[new_segment]) + "-"
            else:
                appending_str += str(segment_to_index[segment]) + "+"
            str2.append(appending_str)

    return "\t".join(str1), "\t".join(str2)
