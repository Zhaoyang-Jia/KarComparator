from NW_aligner import *
from Structures import *

from scipy.optimize import linear_sum_assignment
import numpy as np


def hungarian_algorithm_for_cluster(path_list1, path_list2, index_to_segment, verbose=False):
    """
    Use Hungarian Algorithm to find min_cost perfect bipartite matching
    If input not the same size,
    the smaller side will have dummy nodes with edge weight = len(path excluding forbidden regions)
    :return:
    """
    # add dummy, 0 if none, 1 if on path_list1 side, 2 if on 2's side
    dummy_path_index = []
    if len(path_list1) == len(path_list2):
        dummy_side = 0
    elif len(path_list1) < len(path_list2):
        dummy_side = 1
    else:
        dummy_side = 2

    n_dummy = abs(len(path_list2) - len(path_list1))
    if dummy_side == 1:
        dummy_path_index = [len(path_list1) + i for i in range(n_dummy)]
        for i in range(n_dummy):
            path_list1.append(Path(Arm([], 'dummy'), path_name='dummy'))
    elif dummy_side == 2:
        dummy_path_index = [len(path_list2) + i for i in range(n_dummy)]
        for i in range(n_dummy):
            path_list2.append(Path(Arm([], 'dummy'), path_name='dummy'))

    if len(path_list2) - len(path_list1) != 0:
        raise RuntimeError('bipartite sizes not equal')

    # form cost_matrix
    cost_matrix = [[-1 for col_ind in range(len(path_list2))] for row_ind in range(len(path_list1))]
    alignment_matrix = [[['', ''] for col_ind in range(len(path_list2))] for row_ind in range(len(path_list1))]
    for row_ind in range(len(path_list1)):
        for col_ind in range(len(path_list2)):
            # it is not possible that both sides have dummies
            if path_list1[row_ind].path_name == 'dummy':
                score = path_list2[col_ind].nonforbidden_len()
                alignment_1 = "skipped"
                alignment_2 = "-"
            elif path_list2[col_ind].path_name == 'dummy':
                score = path_list1[row_ind].nonforbidden_len()
                alignment_1 = "-"
                alignment_2 = "skipped"
            else:
                score, alignment_1, alignment_2 = align_paths(path_list1[row_ind].linear_path.segments,
                                                              path_list2[col_ind].linear_path.segments)
                score = score * -1

            cost_matrix[row_ind][col_ind] = score
            alignment_matrix[row_ind][col_ind][0] = alignment_1
            alignment_matrix[row_ind][col_ind][1] = alignment_2

    # optimize for bipartite matching
    np_cost_matrix = np.array(cost_matrix)
    opt_row_ind, opt_col_ind = linear_sum_assignment(np_cost_matrix)

    # tostring
    if verbose:
        print('total cost: {}'.format(np_cost_matrix[opt_row_ind, opt_col_ind].sum()))
        for i in range(len(opt_col_ind)):
            print("alignment {}: {}, {} \t cost: {}".format(i,
                                                            path_list1[opt_row_ind[i]].path_name,
                                                            path_list2[opt_col_ind[i]].path_name,
                                                            cost_matrix[opt_row_ind[i]][opt_col_ind[i]]))
            str1, str2 = tostring_alignment(index_to_segment,
                                            alignment_matrix[opt_row_ind[i]][opt_col_ind[i]][0],
                                            alignment_matrix[opt_row_ind[i]][opt_col_ind[i]][1])
            print(str1)
            print(str2)

    return opt_row_ind, opt_col_ind, cost_matrix, alignment_matrix

