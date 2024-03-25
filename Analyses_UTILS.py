from COMPARISON_with_graphs import *

data_folder = 'cluster_files_testbuild5/'


def form_graph(input_df_row):
    file_prefix = data_folder
    file_suffix = '.txt'
    cluster_file_name = file_prefix + input_df_row['file_name'] + 'cluster_' + str(input_df_row['cluster']) + file_suffix
    graph = form_graph_from_cluster(cluster_file_name)
    return graph


def iterative_get_segment_distance(df_row):
    graph = form_graph(df_row)

    return graph.get_segment_distance()


def iterative_remove_approximated_edges(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()

    n_approximated_edges = graph.karsim_n_transition_approximated \
                           + graph.omkar_n_transition_approximated
    approximated_cnv_distance = graph.approximated_cnv

    return approximated_cnv_distance, n_approximated_edges


def iterative_get_missed_transition_edges(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()

    missed_transition_edges = graph.get_missed_transition_edges()
    return len(missed_transition_edges[0]), len(missed_transition_edges[1])


def iterative_get_cnv(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()

    return graph.get_segment_distance()


def iterative_check_labeled_edges_in_residual_graph(df_row):
    graph = form_graph(df_row)
    graph.prune_same_edges()
    graph.remove_approximate_transition_edges()
    graph.match_transition_edges()

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
