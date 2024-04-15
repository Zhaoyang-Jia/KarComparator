from COMPARISON_with_graphs import *
from Karsimulator_Start_Genome import get_event_chr, get_history_events
from debug_omkar import *

import os
import pandas as pd
import re

data_folder = 'new_data_files/cluster_files_testbuild8/'


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

    omkar_residual_segments, omkar_residual_transitions = graph.gather_edges('omkar')
    karsim_residual_segments, karsim_residual_transitions = graph.gather_edges('karsim')
    labeled_edges = graph.edges_of_interest

    total_dummies_introduced = 0
    dummies_in_residual = 0
    significant_dummies = 0
    significant_dummies_in_residual = 0
    total_seg_introduced = 0
    seg_in_residual = 0
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
        elif edge_type == 'S':
            total_seg_introduced += abs(multiplicity)
            if multiplicity > 0:
                # we added SEG during ILP
                if (start_node, end_node) in omkar_residual_segments:
                    seg_in_residual += min(multiplicity, omkar_residual_segments[(start_node, end_node)])
            elif multiplicity < 0:
                # we removed SEG during ILP
                if (start_node, end_node) in karsim_residual_segments:
                    multiplicity = -1 * multiplicity
                    seg_in_residual += min(multiplicity, karsim_residual_segments[(start_node, end_node)])
            else:
                raise RuntimeError('multiplicity == 0 makes no sense')

    return total_dummies_introduced, dummies_in_residual, total_seg_introduced, seg_in_residual, significant_dummies, significant_dummies_in_residual


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
    omkar_log_dir = 'batch_processing/OMKar_testbuild3/'
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
        new_data = {}
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
        alignment_file = file.split('.')[0] + '.alignment.txt'
        with open('new_data_files/alignment_files/' + alignment_file) as fp_read:
            line1 = fp_read.readline()
            line1 = line1.replace('\n', '').split(': ')[1]
            new_data['total_alignment_cost'] = int(line1)
            alignment_costs = []
            for line in fp_read:
                if line.startswith('alignment'):
                    line = line.replace('\n', '').split('cost: ')[1]
                    alignment_costs.append(int(line))
            new_data['alignment_costs'] = alignment_costs
        data.append(new_data)

    df = pd.DataFrame(data)

    ## get event_chr
    karsim_file_prefix = 'new_data_files/KarSimulator/'
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
    df[['n_dummies', 'n_leftover_dummies', 'n_seg_changed', 'n_leftover_segchange', 'n_significant_dummies', 'n_leftover_significant_dummies']] \
        = df.apply(lambda row: pd.Series(iterative_check_labeled_edges_in_residual_graph(row)), axis=1)

    df['dummy_distance'] = df.apply(lambda row: iterative_get_dummy_lengths(row), axis=1)
    df['n_initial_SVs'] = df.apply(lambda row: pd.Series(iterative_get_initial_n_SV(row)), axis=1)

    return df


if __name__ == "__main__":
    f = '23X_15q26_overgrowth_r1'
    iterative_check_missed_SV_in_preILP(f, [])
