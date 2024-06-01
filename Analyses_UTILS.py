from COMPARISON_with_graphs import *
from Karsimulator_Start_Genome import get_event_chr, get_history_events
from debug_omkar import *
from read_KarSimulator_output import *
from Structures import *

import os
import pandas as pd
import numpy as np
import re

### Can be overwritten during IMPORT
data_folder = '/media/zhaoyang-new/workspace/KarSim/KarComparator/omkar_analyses_pipeline/builds/b14/cluster_files/'
omkar_log_dir = '/media/zhaoyang-new/workspace/KarSim/KarComparator/omkar_analyses_pipeline/builds/b14/omkar_output/'
karsim_file_prefix = 'new_data_files/KarSimulator/'
karsim_history_edges_folder = 'packaged_data/Karsimulator_history_intermediate_terminal_labeled/'
forbidden_region_file = 'Metadata/acrocentric_telo_cen.bed'


def form_graph(input_df_row):
    file_prefix = data_folder
    file_suffix = '.txt'
    cluster_file_name = file_prefix + input_df_row['file_name'] + 'cluster_' + str(input_df_row['cluster']) + file_suffix
    graph = form_graph_from_cluster(cluster_file_name)
    return graph


def iterative_get_segment_distance(df_row):
    graph = form_graph(df_row)
    graph.remove_forbidden_nodes(forbidden_region_file)

    return graph.get_segment_distance()


def iterative_remove_approximated_edges(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()
    graph.remove_forbidden_nodes(forbidden_region_file)

    n_approximated_edges = graph.karsim_n_transition_approximated + graph.omkar_n_transition_approximated
    approximated_cnv_distance = graph.approximated_cnv

    return approximated_cnv_distance, n_approximated_edges


def iterative_get_missed_transition_edges(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()
    graph.remove_forbidden_nodes(forbidden_region_file)

    missed_transition_edges = graph.get_missed_transition_edges()
    return missed_transition_edges[0], len(missed_transition_edges[0]), len(missed_transition_edges[1])


def iterative_get_cnv(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()
    graph.remove_forbidden_nodes(forbidden_region_file)

    return graph.get_segment_distance()


def iterative_check_labeled_edges_in_residual_graph(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()
    graph.remove_forbidden_nodes(forbidden_region_file)
    node_name_to_pos_dict = reverse_dict(graph.node_name)

    omkar_residual_segments, omkar_residual_transitions = graph.gather_edges('omkar')
    karsim_residual_segments, karsim_residual_transitions = graph.gather_edges('karsim')
    labeled_edges = graph.edges_of_interest

    total_dummies_introduced = 0
    dummies_in_residual = 0
    significant_dummies = 0
    significant_dummies_in_residual = 0
    significant_dummies_edges = []
    total_seg_introduced = 0
    significant_seg_introduced = 0
    seg_in_residual = 0
    significant_seg_in_residual = 0
    significant_seg_edges = []  # TODO: also report the CNV change by ILP/Graph formation
    for edge in labeled_edges:
        start_node = edge[0]
        end_node = edge[1]
        multiplicity = edge[2]
        edge_type = edge[3]
        edge_distance = edge[4]
        if edge_type == 'D':
            total_dummies_introduced += multiplicity
            if edge_distance > 200000 or edge_distance == -1:
                significant_dummies += 1
            if (start_node, end_node) in omkar_residual_transitions:
                dummies_in_residual += min(multiplicity, omkar_residual_transitions[(start_node, end_node)])
                if edge_distance > 200000 or edge_distance == -1:
                    significant_dummies_in_residual += 1
            else:
                if edge_distance > 200000 or edge_distance == -1:
                    start_node_chr, start_node_pos = node_name_to_pos_dict[start_node]
                    end_node_chr, end_node_pos = node_name_to_pos_dict[end_node]
                    significant_dummies_edges.append((start_node_chr, end_node_chr, start_node_pos, end_node_pos))
        elif edge_type == 'S':
            total_seg_introduced += abs(multiplicity)
            if edge_distance > 200000:
                significant_seg_introduced += abs(multiplicity)
            if multiplicity > 0:
                # we added SEG during ILP
                if (start_node, end_node) in omkar_residual_segments:
                    seg_in_residual += min(multiplicity, omkar_residual_segments[(start_node, end_node)])
                    if edge_distance > 200000:
                        significant_seg_in_residual += min(multiplicity, omkar_residual_segments[(start_node, end_node)])
            elif multiplicity < 0:
                # we removed SEG during ILP
                if (start_node, end_node) in karsim_residual_segments:
                    multiplicity = -1 * multiplicity
                    seg_in_residual += min(multiplicity, karsim_residual_segments[(start_node, end_node)])
                    if edge_distance > 200000:
                        significant_seg_in_residual += min(multiplicity, karsim_residual_segments[(start_node, end_node)])
            else:
                raise RuntimeError('multiplicity == 0 makes no sense')

    return total_dummies_introduced, \
           dummies_in_residual, \
           total_seg_introduced, \
           seg_in_residual, \
           significant_dummies, \
           significant_dummies_in_residual, \
           significant_seg_introduced, \
           significant_seg_in_residual, \
           significant_dummies_edges


def iterative_get_dummy_lengths(df_row):
    """
    return all lengths (-1 denote of between different chrs)
    :param df_row:
    :return:
    """
    graph = form_graph(df_row)
    return_arr = []
    for edge in graph.edges_of_interest:
        if edge[3] == 'D':
            edge_distance = edge[4]
            return_arr.append(edge_distance)
    return return_arr


def iterative_get_initial_n_SV(df_row):
    graph = form_graph(df_row)
    graph.remove_approximate_transition_edges()
    graph.remove_forbidden_nodes(forbidden_region_file)
    karsim_residual_segments, karsim_residual_transitions = graph.gather_edges('karsim')

    n_transitions = 0
    for edge, multiplicity in karsim_residual_transitions.items():
        n_transitions += multiplicity

    return n_transitions


def iterative_count_events(df_row):
    event_dict = df_row['histories']
    tot = 0
    for key, value in event_dict.items():
        tot += value
    return tot


def iterative_check_missed_SV_in_preILP(input_file_name, missed_SVs, d=200000):
    node_file = omkar_log_dir + input_file_name + '.1/' + input_file_name + '.1.preILP_nodes.txt'
    edge_file = omkar_log_dir + input_file_name + '.1/' + input_file_name + '.1.preILP_edges.txt'
    V = get_vertices_pre_ILP(node_file)
    E = get_edges_pre_ILP(edge_file)

    translated_sv_edges = []
    for edge in E:
        start_node = V[edge[0]]
        end_node = V[edge[1]]
        edge_type = edge[3]
        if edge_type != 'SV':
            continue

        start_chr = start_node.origin_chr
        end_chr = end_node.origin_chr
        start_pos = start_node.pos
        end_pos = end_node.pos
        if start_chr == '23':
            start_chr = 'ChrX'
        elif start_chr == '24':
            start_chr = 'ChrY'
        else:
            start_chr = 'Chr' + start_chr
        if end_chr == '23':
            end_chr = 'ChrX'
        elif end_chr == '24':
            end_chr = 'ChrY'
        else:
            end_chr = 'Chr' + end_chr

        translated_sv_edges.append((start_chr, start_pos, end_chr, end_pos))

    similar_sv = []
    for missed_SV in missed_SVs:
        e1_start_chr = missed_SV[0]
        e1_start_pos = missed_SV[1]
        e1_end_chr = missed_SV[2]
        e1_end_pos = missed_SV[3]
        for preILP_SV in translated_sv_edges:
            e2_start_chr = preILP_SV[0]
            e2_start_pos = preILP_SV[1]
            e2_end_chr = preILP_SV[2]
            e2_end_pos = preILP_SV[3]
            if (e1_start_chr != e2_start_chr) or (e1_end_chr != e2_end_chr):
                continue
            if (abs(e1_start_pos - e2_start_pos) < d * 2) and (abs(e1_end_pos - e2_end_pos) < d * 2):
                similar_sv.append(preILP_SV)
                break

    return similar_sv, len(similar_sv)


def prep_df():
    files = []
    for file in os.listdir(data_folder):
        files.append(file)
    files.sort()

    ## Extract basic cluster file's info
    data = []
    for file in files:
        file_name = file.split('cluster')[0]
        cluster_number = file.split('cluster_')[1].replace('.txt', '')
        with open(data_folder + file) as fp_read:
            line1 = fp_read.readline()
            matches = re.findall(r'<(.*?)>', line1)
            origins = eval(matches[1])
            new_data = {'file_name': file_name,
                        'cluster': matches[0],
                        'n_origin_chr': len(origins),
                        'origin_chr': origins,
                        'n_path_karsim': int(matches[2]),
                        'n_path_omkar': int(matches[3])}
        # alignment_file = file.split('.')[0] + '.alignment.txt'
        # with open(alignment_file_prefix + alignment_file) as fp_read:
        #     line1 = fp_read.readline()
        #     line1 = line1.replace('\n', '').split(': ')[1]
        #     new_data['total_alignment_cost'] = int(line1)
        #     alignment_costs = []
        #     for line in fp_read:
        #         if line.startswith('alignment'):
        #             line = line.replace('\n', '').split('cost: ')[1]
        #             alignment_costs.append(int(line))
        #     new_data['alignment_costs'] = alignment_costs
        data.append(new_data)

    df = pd.DataFrame(data)

    ## get event_chr
    df['event_chr'] = df['file_name'].apply(lambda x: list(get_event_chr(karsim_file_prefix + x + '.kt.txt')))
    df['event_chr'] = df['event_chr'].apply(lambda x: [entry[:-1] for entry in x])

    df['histories'] = df.apply(lambda x: get_history_events(karsim_file_prefix + x['file_name'] + '.kt.txt',
                                                            x['origin_chr']), axis=1)
    df['n_events'] = df.apply(lambda row: iterative_count_events(row), axis=1)
    df['is_event_cluster'] = df.apply(lambda row: any(origin_chr in row['event_chr'] for origin_chr in row['origin_chr']), axis=1)

    ## n_files with different number of Chr
    df['n_path_diff'] = df['n_path_omkar'] - df['n_path_karsim']

    return df


def process_comparison(df):
    df[['approximated_distance', 'n_approximated_edges']] \
        = df.apply(lambda row: pd.Series(iterative_remove_approximated_edges(row)), axis=1)
    df[['karsim_missed_transition', 'n_karsim_missed_transition', 'n_omkar_missed_transition']] \
        = df.apply(lambda row: pd.Series(iterative_get_missed_transition_edges(row)), axis=1)
    df['SV_missed'] = df['n_karsim_missed_transition'] + df['n_omkar_missed_transition']
    df['CNV_missed'] = df.apply(lambda row: pd.Series(iterative_get_cnv(row)), axis=1)
    df['log10_CNV_missed'] = np.log10(df['CNV_missed'])
    df[['n_dummies', 'n_leftover_dummies', 'n_seg_changed', 'n_leftover_segchange', 'n_significant_dummies', 'n_leftover_significant_dummies',
        'n_significant_seg_change',
        'n_leftover_significant_seg_change',
        'cancelled_significant_dummies']] \
        = df.apply(lambda row: pd.Series(iterative_check_labeled_edges_in_residual_graph(row)), axis=1)

    df['dummy_distance'] = df.apply(lambda row: iterative_get_dummy_lengths(row), axis=1)
    df['n_initial_SVs'] = df.apply(lambda row: pd.Series(iterative_get_initial_n_SV(row)), axis=1)

    return df


def iterative_label_missed_SV_edges(df_row):
    karsim_filename = df_row['file_name']
    karsim_history_edges_filepath = karsim_history_edges_folder + karsim_filename + '.history_sv.txt'
    event_sv_edges = read_history_edges_intermediate_file(karsim_history_edges_filepath)
    missed_sv_edges = df_row['karsim_missed_transition']
    labeled_event_type = []
    for missed_sv_edge in missed_sv_edges:
        event_found = False
        for entry in event_sv_edges:
            event_type = entry[0]
            event_edges = entry[1]
            if missed_sv_edge in event_edges:
                labeled_event_type.append(event_type)
                event_found = True
                break
        if not event_found:
            labeled_event_type.append('ENF')
    return labeled_event_type


def label_missed_SV_edges(df):
    df['karsim_missed_transition_event_type'] = df.apply(lambda row: iterative_label_missed_SV_edges(row), axis=1)
    return df


def iterative_get_cn_with_bins(df_row, cn_file_name='Metadata/cn_bins_200kbp.txt'):
    cn_bins = read_cn_bin_file(cn_file_name)
    graph = form_graph(df_row)
    karsim_cn, omkar_cn = graph_assign_cn_bin(graph, cn_bins)
    return karsim_cn, omkar_cn


def sum_history_dicts(history_dicts):
    summed_dict = {}
    for history_dict in history_dicts:
        for key, value in history_dict.items():
            if key in summed_dict:
                summed_dict[key] += value
            else:
                summed_dict[key] = value
    return summed_dict


#################CN###############
def read_cn_bin_file(cn_file_name):
    cn_bins = []
    with open(cn_file_name) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            cn_bins.append({'chrom': line[0],
                            'start': int(line[1]),
                            'end': int(line[2])})
    return cn_bins


def graph_assign_cn_bin(input_graph: Graph, input_cn_bins):
    karsim_cn = np.array([0.0 for _ in input_cn_bins])
    omkar_cn = np.array([0.0 for _ in input_cn_bins])
    for start_node, end_nodes in input_graph.karsim_dict.items():
        for end_node in end_nodes:
            edge_type = end_node[2]
            if edge_type == 'segment':
                c_seg = Segment(end_node[0], start_node[1], end_node[1])
                karsim_cn += c_seg.assign_cn_bin(input_cn_bins)
    for start_node, end_nodes in input_graph.omkar_dict.items():
        for end_node in end_nodes:
            edge_type = end_node[2]
            if edge_type == 'segment':
                c_seg = Segment(end_node[0], start_node[1], end_node[1])
                omkar_cn += c_seg.assign_cn_bin(input_cn_bins)
    return karsim_cn, omkar_cn


def cn_bin_value_boolean_conversion(input_cn, expected_count=2, rounding_allowance=0.01):
    """
    a bin will be called '1'/changed if deviated from WT
    :param rounding_allowance:
    :param expected_count: WT CN
    :param input_cn:
    :return:
    """
    output_boolean_cn = []
    for bin_itr in input_cn:
        if abs(bin_itr - expected_count) <= rounding_allowance:
            output_boolean_cn.append(0)
        else:
            output_boolean_cn.append(1)
    return np.array(output_boolean_cn)


def cn_jaccard_similarity(cn1, cn2):
    bool_array1 = cn1.astype(bool)
    bool_array2 = cn2.astype(bool)

    intersection = np.sum(np.logical_and(bool_array1, bool_array2))
    union = np.sum(np.logical_or(bool_array1, bool_array2))
    if union == 0:
        return 1.0  # If both arrays are all zeros, define Jaccard similarity as 1
    else:
        return intersection / union


def cn_jaccard_sim_by_case(dict1: {str: [float]}, dict2):
    if dict1.keys() != dict2.keys():
        raise RuntimeError
    case_jaccard_sim = {}
    for c_file_name, file1_cn in dict1.items():
        file2_cn = dict2[c_file_name]
        file1_bool_cn = cn_bin_value_boolean_conversion(file1_cn)
        file2_bool_cn = cn_bin_value_boolean_conversion(file2_cn)
        case_jaccard_sim[c_file_name] = cn_jaccard_similarity(file1_bool_cn, file2_bool_cn)
    return case_jaccard_sim


def plot_cn_jaccared_sim_by_case(case_jaccard_sim, output_path):
    values = list(case_jaccard_sim.values())
    print('mean jaccard similarity: ', np.array(values).mean())
    plt.figure(figsize=(8, 5))
    plt.hist(values, bins=10, edgecolor=(0, 0, 0, 0.5), histtype='bar')
    plt.xlabel('Jaccard Similarity', fontsize=12)
    plt.ylabel('Frequency', fontsize=12)
    plt.title('Histogram of Copy Number Jaccard Similarity, from Simulations (n=28)', fontsize=15, pad=14)
    plt.savefig(output_path, dpi=300)


def plot_cn_jaccared_sim_cos_sim_by_case(case_jaccard_sim, case_cos_sim, output_path):
    values1 = list(case_jaccard_sim.values())
    values2 = list(case_cos_sim.values())
    print('mean jaccard similarity: ', np.array(values1).mean())
    print('mean cos similarity: ', np.array(values2).mean())
    plt.figure(figsize=(6, 5))
    plt.hist(values1, bins=10, alpha=0.7, label='Jaccard Similarity', edgecolor=(0, 0, 0, 0.5))
    plt.hist(values2, bins=10, alpha=0.7, label='Cosine Similarity', color='orange', edgecolor=(0, 0, 0, 0.5))
    plt.xlabel('Similarity Score', fontsize=14)
    plt.ylabel('Frequency', fontsize=14)
    plt.title('Histogram of Affected Copy Number Similarities\nfrom Simulations (n=28)', fontsize=15, pad=9)
    plt.legend(fontsize=14)
    plt.savefig(output_path, dpi=300)


def output_cn_dict_by_case(case_jaccard_sim, output_file_path):
    with open(output_file_path, 'w') as fp_write:
        for case_name, case_jaccard_score in case_jaccard_sim.items():
            fp_write.write('{}\t{}\n'.format(case_name, case_jaccard_score))


def read_cn_jaccard_sim_by_case_and_plot(jaccard_sim_file_path, output_pic_path):
    case_jaccard_sim = {}
    with open(jaccard_sim_file_path) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            case_jaccard_sim[line[0]] = float(line[1])
    plot_cn_jaccared_sim_by_case(case_jaccard_sim, output_pic_path)


def read_cn_jaccard_sim_cos_sim_and_plot(jaccard_sim_file_path, cos_sim_file_path, output_pic_path):
    case_jaccard_sim = {}
    with open(jaccard_sim_file_path) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            case_jaccard_sim[line[0]] = float(line[1])
    case_cos_sim = {}
    with open(cos_sim_file_path) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            case_cos_sim[line[0]] = float(line[1])
    plot_cn_jaccared_sim_cos_sim_by_case(case_jaccard_sim, case_cos_sim, output_pic_path)


def vector_addition(lists):
    result = [sum(x) for x in zip(*lists)]
    return result


def cos_similarity(arr1, arr2):
    dot_product = np.dot(arr1, arr2)
    norm1 = np.linalg.norm(arr1)
    norm2 = np.linalg.norm(arr2)

    cosine_sim = dot_product / (norm1 * norm2)
    return cosine_sim


def cos_sim_by_case(dict1: {str: [float]}, dict2):
    if dict1.keys() != dict2.keys():
        raise RuntimeError
    case_cos_sim = {}
    cos_sim_sum = 0.0
    for c_file_name, file1_cn in dict1.items():
        file2_cn = dict2[c_file_name]
        c_cos_sim = cos_similarity(np.array(file1_cn), np.array(file2_cn))
        case_cos_sim[c_file_name] = c_cos_sim
        cos_sim_sum += c_cos_sim
    avg_cos_sim = cos_sim_sum / len(dict1)
    return case_cos_sim, avg_cos_sim


def union_altered_cn_bin(cn1, cn2):
    cn1_bool = cn_bin_value_boolean_conversion(cn1)
    cn2_bool = cn_bin_value_boolean_conversion(cn2)
    print('cn1: ', cn1_bool.sum())
    print('cn2: ', cn2_bool.sum())
    return np.logical_or(cn1_bool, cn2_bool)


def cos_sim_with_jaccard_union_by_case(dict1: {str: [float]}, dict2):
    """
    for each case, find the bins as the union of the mutated bins, and compute the cosine similarity
    :param dict1:
    :param dict2:
    :return:
    """
    if dict1.keys() != dict2.keys():
        raise RuntimeError
    case_cos_sim = {}
    for c_file_name, file1_cn in dict1.items():
        file2_cn = dict2[c_file_name]
        bin_filter = union_altered_cn_bin(file1_cn, file2_cn)
        print('union: ', bin_filter.sum())
        file1_cn_sublist = [val for idx, val in enumerate(file1_cn) if bin_filter[idx] == 1]
        file2_cn_sublist = [val for idx, val in enumerate(file2_cn) if bin_filter[idx] == 1]
        case_cos_sim[c_file_name] = cos_similarity(np.array(file1_cn_sublist),
                                                   np.array(file2_cn_sublist))
    return case_cos_sim


def cos_sim_stacked(dict1: {str: [float]}, dict2):
    if dict1.keys() != dict2.keys():
        raise RuntimeError
    file1_concat_cn = []
    file2_concat_cn = []
    for c_file_name, file1_cn in dict1.items():
        file2_cn = dict2[c_file_name]
        file1_concat_cn += file1_cn
        file2_concat_cn += file2_cn
    return cos_similarity(np.array(file1_concat_cn), np.array(file2_concat_cn))


##############################################
################ACCURACY BY SV TYPES##########
def iterative_label_graph_with_sv_types(df_row):
    karsim_filename = df_row['file_name']
    karsim_history_edges_filepath = karsim_history_edges_folder + karsim_filename + '.history_sv.txt'
    graph = form_graph(df_row)
    graph.add_sv_label_to_karsim_edges(karsim_history_edges_filepath)
    sv_edge_event_types_tally1 = graph.tally_karsim_edge_event_types()

    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()
    sv_edge_event_types_tally2 = graph.tally_karsim_edge_event_types()

    graph.remove_forbidden_nodes(forbidden_region_file)
    sv_edge_event_types_tally3 = graph.tally_karsim_edge_event_types()

    return sv_edge_event_types_tally1, sv_edge_event_types_tally2, sv_edge_event_types_tally3


def accuracy_by_event_types():
    df = prep_df()
    df[['karsim_edge_type_tally_prematching',
        'karsim_edge_type_tally_postmatching',
        'karsim_edge_type_tally_postforbiddenremoval']] = \
        df.apply(lambda row: iterative_label_graph_with_sv_types(row), axis=1, result_type='expand')

    tally_dicts1 = sum_history_dicts(df['karsim_edge_type_tally_prematching'].tolist())
    tally_dicts2 = sum_history_dicts(df['karsim_edge_type_tally_postmatching'].tolist())
    tally_dicts3 = sum_history_dicts(df['karsim_edge_type_tally_postforbiddenremoval'].tolist())

    print(tally_dicts1)
    print(tally_dicts2)
    print(tally_dicts3)


################################################

def test_label_missed_sv_edges():
    i_df = prep_df()
    i_df = process_comparison(i_df)
    SV_missed = i_df['n_karsim_missed_transition'].sum()
    SV_added = i_df['n_omkar_missed_transition'].sum()
    i_df[['preILP_similar_SV', 'n_preILP_similar_SV']] \
        = i_df.apply(lambda row: pd.Series(iterative_check_missed_SV_in_preILP(row['file_name'], row['karsim_missed_transition'])), axis=1)
    i_df = i_df.sort_values('n_karsim_missed_transition', ascending=False)
    i_df = label_missed_SV_edges(i_df)
    i_df.to_csv('missed_sv_edges_labeled.csv')


def test_cnv_plot():
    read_cn_jaccard_sim_cos_sim_and_plot(
        '/media/zhaoyang-new/workspace/KarSim/KarComparator/omkar_analyses_pipeline/builds/b14/analyses_summary/cn_jaccard_scores.txt',
        '/media/zhaoyang-new/workspace/KarSim/KarComparator/omkar_analyses_pipeline/builds/b14/analyses_summary/cn_cos_scores.txt',
        '/media/zhaoyang-new/workspace/KarSim/KarComparator/omkar_analyses_pipeline/builds/b14/analyses_summary/cn_jaccard_cos_plot.png')


def test_graph_label_sv_types():
    i_df = prep_df()


if __name__ == "__main__":
    accuracy_by_event_types()
