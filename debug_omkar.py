import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
import argparse
import ast


from forbidden_region_processing import read_forbidden_regions
from Structures import *


dir_name = '/media/zhaoyang-new/workspace/KarSim/KarComparator/batch_processing/omkar_output_temp/'
forbidden_region_file = '/media/zhaoyang-new/workspace/KarSim/KarComparator/Metadata/acrocentric_telo_cen.bed'

uniform_dist = 1 / 6
graph_size_width_multiplier = 20
graph_width = uniform_dist * graph_size_width_multiplier
width_ratio_multiplier = 2.5
height_ratio_multiplier = 1


class Vertex:
    name: int
    origin_chr: str
    pos: int
    copy_num: int
    neighbors: [int]

    def __init__(self, name, origin_chr, pos):
        self.name = name
        self.origin_chr = origin_chr
        self.pos = pos


def get_vertices_pre_ILP(log_file):
    V = {}
    with open(log_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            pos = int(line[2].split('.')[0])
            V[line[0]] = Vertex(int(line[0]), line[1], pos)
    return V


def get_edges_pre_ILP(log_file):
    E = []
    with open(log_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').replace('(', '').replace(')', '').replace("'", '').split(', ')
            E.append(tuple(line))
    return E


def get_edges_post_ILP(log_file):
    E = []
    with open(log_file) as fp_read:
        fp_read.readline()  # skip first line
        for line in fp_read:
            line = line.replace('\n', '').replace('(', '').replace(')', '').replace("'", '').split(', ')
            E.append(tuple(line))
    return E


def get_dummy_edges(log_file):
    E = []
    with open(log_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').replace('(', '').replace(')', '').replace("'", '').split(', ')
            if line[3] == 'D':
                # only include the dummy edges
                E.append(tuple(line))
    return E


def consolidate_edge_with_dummy_edge(E, dummy_E):
    for dummy_E_itr in dummy_E:
        val0 = dummy_E_itr[0]
        val1 = dummy_E_itr[1]

        same_edge_found = False
        for E_itr in E:
            if E_itr[0] == val0 and E_itr[1] == val1:
                same_edge_found = True
                break

        if not same_edge_found:
            E.append(dummy_E_itr)
        else:
            # test if the reverse of the dummy edge is also in E
            for E_itr in E:
                if E_itr[0] == val1 and E_itr[1] == val0:
                    raise RuntimeError('reverse of dummy edge also in file, preventing overwrite')
            reversed_dummy = tuple([dummy_E_itr[1], dummy_E_itr[0], dummy_E_itr[2], dummy_E_itr[3]])
            E.append(reversed_dummy)
    return E


def filter_nodes(chrs_of_interest, V):
    key_to_pop = []
    for key in V:
        if V[key].origin_chr not in chrs_of_interest:
            key_to_pop.append(key)
    for key in key_to_pop:
        V.pop(key)
    return V


def filter_edges(filtered_V, E):
    vertex_of_interest = list(filtered_V.keys())
    ind_to_keep = []
    for edge_ind in range(len(E)):
        node1 = E[edge_ind][0]
        node2 = E[edge_ind][1]
        if node1 in vertex_of_interest or node2 in vertex_of_interest:
            ind_to_keep.append(edge_ind)

    filtered_E = [E[i] for i in range(len(E)) if i in ind_to_keep]
    return filtered_E


def iterative_add_edge(edges_list: [(str, str, str, str)], graph):
    for edge in edges_list:
        node1 = edge[0]
        node2 = edge[1]
        multiplicity = int(edge[2])
        edge_type = edge[3]

        edge_color = None
        if edge_type == 'S':
            edge_color = 'black'
        elif edge_type == 'R':
            edge_color = 'blue'
        elif edge_type == 'SV':
            edge_color = 'red'
        elif edge_type == 'D':
            edge_color = 'orange'

        graph.add_edge(node1, node2, color=edge_color, weight=multiplicity)


def generate_uniform_linear_coordinates(n, fixed_y=0.5, fixed_distance=uniform_dist):
    coordinates = []
    for i in range(n):
        x = round(0 + i * fixed_distance, 4)
        y = fixed_y
        coordinates.append((x, y))

    return coordinates


def label_centromere_nodes(V, forbidden_file):
    """
    archived, should label ref edge as CEN edge,instead of nodes
    :param V:
    :param forbidden_file:
    :return:
    """
    forbidden_region_segments = read_forbidden_regions(forbidden_file).segments
    centromere_segments = []
    for segment in forbidden_region_segments:
        if 'centromere' in segment.segment_type:
            centromere_segments.append(segment)

    centromere_V = []
    for node_id, node in V.items():
        current_chr = 'Chr' + node.origin_chr
        current_chr_centromere_segment = None
        for segment in centromere_segments:
            if segment.chr_name == current_chr:
                current_chr_centromere_segment = segment
                break
        if current_chr_centromere_segment.start <= node.pos <= current_chr_centromere_segment.end:
            centromere_V.append(node_id)
    print('x')


def graph_pre_ILP(chrs_of_interest, header_name, output_dir):
    node_file = dir_name + header_name + '.1/' + header_name + '.1.preILP_nodes.txt'
    edge_file = dir_name + header_name + '.1/' + header_name + '.1.preILP_edges.txt'

    G_preILP = nx.DiGraph()
    V = get_vertices_pre_ILP(node_file)
    filtered_V = filter_nodes(chrs_of_interest, V)
    G_preILP.add_nodes_from(filtered_V)

    E = get_edges_pre_ILP(edge_file)
    filtered_E = filter_edges(filtered_V, E)
    iterative_add_edge(filtered_E, G_preILP)

    V_pos = {}
    V_names = list(V.keys())
    cor = generate_uniform_linear_coordinates(len(V))
    for node_ind in range(len(V)):
        V_pos[V_names[node_ind]] = cor[node_ind]

    E_colors = nx.get_edge_attributes(G_preILP, 'color')
    E_weights = nx.get_edge_attributes(G_preILP, 'weight')

    plt.figure(figsize=(graph_width * width_ratio_multiplier, graph_width * height_ratio_multiplier))
    plot_karsim = ng.InteractiveGraph(G_preILP,
                                      node_layout=V_pos,
                                      node_labels=True,
                                      edge_color=E_colors,
                                      edge_layout='arc',
                                      edge_labels=E_weights,
                                      arrows=True,
                                      node_size=6,
                                      node_label_offset=0.001,
                                      node_label_font_dict=dict(size=10),
                                      edge_label_fontdict=dict(size=9),
                                      scale=(graph_width, 1))
    plt.savefig(output_dir + header_name + str(chrs_of_interest) + '.preILP.png')


def graph_post_ILP(chrs_of_interest, header_name, output_dir):
    node_file = dir_name + header_name + '.1/' + header_name + '.1.preILP_nodes.txt'
    metadata_file = dir_name + header_name + '.1/postILP_components/' + header_name + '.1.postILP.metadata.txt'

    V = get_vertices_pre_ILP(node_file)
    filtered_V = filter_nodes(chrs_of_interest, V)
    # label_centromere_nodes(V, forbidden_region_file)

    ## find the right cluster file
    cluster_metadata = {}
    with open(metadata_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            cluster_number = line[0]
            cluster_nodes = [str(i) for i in eval(line[1])]
            cluster_metadata[cluster_number] = cluster_nodes
    # assumes that only one cluster is in play
    cluster_to_graph = None
    for key, value in cluster_metadata.items():
        if list(filtered_V.keys())[0] in value:
            cluster_to_graph = key
            break

    edge_file = dir_name + header_name + '.1/postILP_components/' + header_name + '.1.postILP_component_' + cluster_to_graph + '.txt'
    E = get_edges_post_ILP(edge_file)

    G_postILP = nx.DiGraph()
    G_postILP.add_nodes_from(filtered_V)
    iterative_add_edge(E, G_postILP)

    V_pos = {}
    V_names = list(V.keys())
    cor = generate_uniform_linear_coordinates(len(V))
    for node_ind in range(len(V)):
        V_pos[V_names[node_ind]] = cor[node_ind]

    E_colors = nx.get_edge_attributes(G_postILP, 'color')
    E_weights = nx.get_edge_attributes(G_postILP, 'weight')

    plt.figure(figsize=(graph_width * width_ratio_multiplier, graph_width * height_ratio_multiplier))
    plot_karsim = ng.InteractiveGraph(G_postILP,
                                      node_layout=V_pos,
                                      node_labels=True,
                                      edge_color=E_colors,
                                      edge_layout='arc',
                                      edge_labels=E_weights,
                                      arrows=True,
                                      node_size=6,
                                      node_label_offset=0.001,
                                      node_label_font_dict=dict(size=10),
                                      edge_label_fontdict=dict(size=9),
                                      scale=(graph_width, 1))
    plt.savefig(output_dir + header_name + str(chrs_of_interest) + '.postILP.png')


def graph_post_ILP_with_dummies(chrs_of_interest, header_name, output_dir):
    node_file = dir_name + header_name + '.1/' + header_name + '.1.preILP_nodes.txt'
    metadata_file = dir_name + header_name + '.1/postILP_components/' + header_name + '.1.postILP.metadata.txt'

    V = get_vertices_pre_ILP(node_file)
    filtered_V = filter_nodes(chrs_of_interest, V)
    # label_centromere_nodes(V, forbidden_region_file)

    ## find the right cluster file
    cluster_metadata = {}
    with open(metadata_file) as fp_read:
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            cluster_number = line[0]
            cluster_nodes = [str(i) for i in eval(line[1])]
            cluster_metadata[cluster_number] = cluster_nodes
    # assumes that only one cluster is in play
    cluster_to_graph = None
    for key, value in cluster_metadata.items():
        if list(filtered_V.keys())[0] in value:
            cluster_to_graph = key
            break

    edge_file = dir_name + header_name + '.1/postILP_components/' + header_name + '.1.postILP_component_' + cluster_to_graph + '.txt'
    dummy_file = dir_name + header_name + '.1/all_edges_with_dummies/' + header_name + '.1.with_dummies_component_' + cluster_to_graph + '.txt'
    E = get_edges_post_ILP(edge_file)
    dummy_E = get_dummy_edges(dummy_file)

    # prevent dummy_edge overwriting regular edge by having duplicate (multi-graph not allowed)
    E = consolidate_edge_with_dummy_edge(E, dummy_E)

    G_postILP = nx.DiGraph()
    G_postILP.add_nodes_from(filtered_V)
    iterative_add_edge(E, G_postILP)

    V_pos = {}
    V_names = list(V.keys())
    cor = generate_uniform_linear_coordinates(len(V))
    for node_ind in range(len(V)):
        V_pos[V_names[node_ind]] = cor[node_ind]

    E_colors = nx.get_edge_attributes(G_postILP, 'color')
    E_weights = nx.get_edge_attributes(G_postILP, 'weight')

    plt.figure(figsize=(graph_width * width_ratio_multiplier, graph_width * height_ratio_multiplier))
    plot_karsim = ng.InteractiveGraph(G_postILP,
                                      node_layout=V_pos,
                                      node_labels=True,
                                      edge_color=E_colors,
                                      edge_layout='arc',
                                      edge_labels=E_weights,
                                      arrows=True,
                                      node_size=6,
                                      node_label_offset=0.001,
                                      node_label_font_dict=dict(size=10),
                                      edge_label_fontdict=dict(size=9),
                                      scale=(graph_width, 1))
    plt.savefig(output_dir + header_name + str(chrs_of_interest) + '.with_dummy.png')


if __name__ == "__main__":
    def parse_list(s):
        x = ast.literal_eval(s)
        new_list = []
        for item in x:
            new_list.append(str(item))
        return new_list

    parser = argparse.ArgumentParser(description='plotting omkar preILP and postILP graphs')
    parser.add_argument('--output_dir', type=str)
    parser.add_argument('--case_file', type=str)
    parser.add_argument('--chr_of_int', type=parse_list)

    args = parser.parse_args()
    chr_of_int = args.chr_of_int
    current_output_dir = args.output_dir
    header = args.case_file

    # chr_of_int = ['13']
    # current_output_dir = '/media/zhaoyang-new/workspace/KarSim/KarComparator/debug_omkar/'
    # header = '23X_15q26_overgrowth_r1'
    graph_pre_ILP(chr_of_int, header, current_output_dir)
    graph_post_ILP(chr_of_int, header, current_output_dir)
    graph_post_ILP_with_dummies(chr_of_int, header, current_output_dir)
