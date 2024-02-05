import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
import math
import numpy as np
from scipy.optimize import linear_sum_assignment
from collections import Counter

from Structures import *
from read_cluster_file import *


class Graph:
    node_name: {(str, int): str}  # chr, position: node name
    edge_type: {(str, int, str, int): str}  # chr1, position1, chr2, position2: edge type
    segment_edge_type: {(str, int, str, int): str}  # documents for comparison forbidden regions
    karsim_dict: {(str, int): [(str, int)]}
    omkar_dict: {(str, int): [(str, int)]}

    def __init__(self):
        self.node_name = {}
        self.edge_type = {}
        self.segment_edge_type = {}
        self.karsim_dict = {}
        self.omkar_dict = {}

    def add_edge_to_dict(self, chr1, pos1, chr2, pos2, target_dict: str):
        if target_dict == 'omkar':
            target_dict = self.omkar_dict
        elif target_dict == 'karsim':
            target_dict = self.karsim_dict
        else:
            raise ValueError()

        if (chr1, pos1) in target_dict:
            target_dict[(chr1, pos1)].append((chr2, pos2))
        else:
            target_dict[(chr1, pos1)] = [(chr2, pos2)]

    def add_segment_edge(self, input_segment: Segment, target_graph: str):
        """
        :param input_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        # add node
        self.node_name[(input_segment.chr_name, input_segment.start)] = input_segment.kt_index[:-1] + "$"
        self.node_name[(input_segment.chr_name, input_segment.end)] = input_segment.kt_index[:-1] + "#"

        # document edge type
        self.edge_type[(input_segment.chr_name, input_segment.start,
                        input_segment.chr_name, input_segment.end)] = 'segment'
        self.segment_edge_type[(input_segment.chr_name, input_segment.start,
                                input_segment.chr_name, input_segment.end)] = input_segment.segment_type

        # add edge to dict
        self.add_edge_to_dict(input_segment.chr_name, input_segment.start,
                              input_segment.chr_name, input_segment.end,
                              target_graph)

    def add_transition_edge(self, start_segment: Segment, end_segment: Segment, target_graph: str):
        """
        transition edge from the end of the start_segment to the start of the end_segment
        :param start_segment:
        :param end_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        # document edge type
        self.edge_type[(start_segment.chr_name, start_segment.end,
                        end_segment.chr_name, end_segment.start)] = 'transition'

        # add edge to dict
        self.add_edge_to_dict(start_segment.chr_name, start_segment.end,
                              end_segment.chr_name, end_segment.start,
                              target_graph)

    def prune_same_edges(self):
        dict1 = self.karsim_dict
        dict2 = self.omkar_dict

        # Find the common keys between the dictionaries
        common_keys = dict1.keys() & dict2.keys()

        # Iterate over common keys and remove instances
        for key in common_keys:
            occurrences_in_dict1 = Counter(dict1[key])
            occurrences_in_dict2 = Counter(dict2[key])

            for value in occurrences_in_dict1:
                # Determine the minimum number of occurrences between the two dicts
                if value in occurrences_in_dict2:
                    min_occurrences = min(occurrences_in_dict1[value], occurrences_in_dict2[value])
                    for _ in range(min_occurrences):
                        dict1[key].remove(value)
                        dict2[key].remove(value)
                else:
                    continue

        # Remove keys that have become empty
        self.karsim_dict = {k: v for k, v in dict1.items() if v}
        self.omkar_dict = {k: v for k, v in dict2.items() if v}

    def gather_edges(self, target_graph: str) -> ({(str, str): int}, {(str, str): int}):
        if target_graph == 'karsim':
            target_dict = self.karsim_dict
        elif target_graph == 'omkar':
            target_dict = self.omkar_dict
        else:
            raise ValueError

        E_segment = {}
        E_transition = {}
        for node1, value in target_dict.items():
            for node2 in value:
                node1_name = self.node_name[(node1[0], node1[1])]
                node2_name = self.node_name[(node2[0], node2[1])]
                edge_type = self.edge_type[(node1[0], node1[1], node2[0], node2[1])]

                if edge_type == 'segment':
                    if (node1_name, node2_name) in E_segment:
                        E_segment[(node1_name, node2_name)] += 1
                    else:
                        E_segment[(node1_name, node2_name)] = 1
                elif edge_type == 'transition':
                    if (node1_name, node2_name) in E_transition:
                        E_transition[(node1_name, node2_name)] += 1
                    else:
                        E_transition[(node1_name, node2_name)] = 1
                else:
                    raise ValueError

        return E_segment, E_transition

    def custom_sort_node(self, node_tuple):
        chr_info = node_tuple[0]
        pos_info = int(node_tuple[1])
        chr_index = chr_info[3:]
        if chr_index == "X":
            chr_index = 23
        elif chr_index == "Y":
            chr_index = 24
        else:
            chr_index = int(chr_index)

        return chr_index, pos_info

    def get_segment_distance(self):
        # reverse the dict
        node_name_to_node_position = {}
        for position, name in self.node_name.items():
            node_name_to_node_position[name] = position

        total_distance = 0
        for node1, value in self.karsim_dict.items():
            for node2 in value:
                if self.edge_type[(node1[0], node1[1], node2[0], node2[1])] == 'transition':
                    continue
                # all segment edge are between the same chromosome nodes
                edge_type = self.segment_edge_type[(node1[0], node1[1], node2[0], node2[1])]
                if edge_type.startswith('telomere') or edge_type.startswith('acrocentric'):
                    continue
                else:
                    total_distance += abs(node2[1] - node1[1] + 1)

        for node1, value in self.omkar_dict.items():
            for node2 in value:
                if self.edge_type[(node1[0], node1[1], node2[0], node2[1])] == 'transition':
                    continue
                # all segment edge are between the same chromosome nodes
                edge_type = self.segment_edge_type[(node1[0], node1[1], node2[0], node2[1])]
                if edge_type.startswith('telomere') or edge_type.startswith('acrocentric'):
                    continue
                else:
                    total_distance += abs(node2[1] - node1[1] + 1)

        return total_distance

    def get_missed_transition_edges(self):
        karsim_transition_edges = []
        for node1, value in self.karsim_dict.items():
            for node2 in value:
                if self.edge_type[(node1[0], node1[1], node2[0], node2[1])] == 'transition':
                    # skip reference edge
                    if node1[0] == node2[0] and node1[1] + 1 == node2[1]:
                        continue
                    else:
                        karsim_transition_edges.append((node1[0], node1[1], node2[0], node2[1]))

        omkar_transition_edges = []
        for node1, value in self.omkar_dict.items():
            for node2 in value:
                if self.edge_type[(node1[0], node1[1], node2[0], node2[1])] == 'transition':
                    # skip reference edge
                    if node1[0] == node2[0] and node1[1] + 1 == node2[1]:
                        continue
                    else:
                        omkar_transition_edges.append((node1[0], node1[1], node2[0], node2[1]))

        return karsim_transition_edges, omkar_transition_edges

    def get_chr_start_end_nodes(self):
        """
        find the two nodes that create a chr-breakpoint (i.e. belong to two diff. chr)
        :return: list of all chr-breakpoint nodes' name
        """
        nodes = self.node_name.keys()
        sorted_nodes = sorted(nodes, key=self.custom_sort_node)  # sorted by chr, then by positions -> chr always start to end sorted in the list
        terminal_nodes = [self.node_name[sorted_nodes[0]]]  # first node always the start of a chr
        node_ind = 1
        while node_ind <= len(sorted_nodes) - 3:
            current_node = sorted_nodes[node_ind]
            next_node = sorted_nodes[node_ind + 1]
            if current_node[0] != next_node[0]:
                terminal_nodes.append(self.node_name[current_node])
                terminal_nodes.append(self.node_name[next_node])
            node_ind += 2  # because of the $# paired structure, only check the # against the next $
        terminal_nodes.append(self.node_name[sorted_nodes[-1]])  # last node always the end of a chr

        return terminal_nodes

    def visualize_graph(self, output_prefix):
        # create sorted nodes (Endpoints)

        nodes = self.node_name.keys()
        sorted_nodes = sorted(nodes, key=self.custom_sort_node)
        V = []
        for node in sorted_nodes:
            V.append(self.node_name[node])

        E_karsim_segment, E_karsim_transition = self.gather_edges('karsim')
        E_omkar_segment, E_omkar_transition = self.gather_edges('omkar')

        def translate_segment_edge_type_to_edge_name():
            new_dict = {}
            for edge, edge_type in self.segment_edge_type.items():
                node1_name = self.node_name[(edge[0], edge[1])]
                node2_name = self.node_name[(edge[2], edge[3])]
                new_dict[(node1_name, node2_name)] = edge_type
            return new_dict

        def iterative_add_edge(edges_dict: {(str, str): int}, group_color, graph, forbidden_segment_edge_labels=False):
            """
            if forbidden_segment_edge_labels=True, prepare "translated_segment_edge_type"
            :param edges_dict:
            :param group_color:
            :param graph:
            :param forbidden_segment_edge_labels:
            :return:
            """
            for edge, edge_weight in edges_dict.items():
                edge_color = group_color
                if forbidden_segment_edge_labels:
                    this_edge_type = translated_segment_edge_type[(edge[0]), (edge[1])]
                    if this_edge_type.startswith('telomere') or this_edge_type.startswith('acrocentric'):
                        edge_color = 'blue'
                graph.add_edge(edge[0], edge[1], color=edge_color, weight=edge_weight)

        G_karsim = nx.DiGraph()
        G_karsim.add_nodes_from(V)

        G_omkar = nx.DiGraph()
        G_omkar.add_nodes_from(V)

        # add all edges
        translated_segment_edge_type = translate_segment_edge_type_to_edge_name()
        iterative_add_edge(E_karsim_segment, 'black', G_karsim, forbidden_segment_edge_labels=True)
        iterative_add_edge(E_karsim_transition, 'red', G_karsim)
        iterative_add_edge(E_omkar_segment, 'black', G_omkar, forbidden_segment_edge_labels=True)
        iterative_add_edge(E_omkar_transition, 'red', G_omkar)

        def generate_circular_coordinates(n, radius=0.5, center=(0.5, 0.5)):
            cx, cy = center
            coordinates = []

            for k in range(n):
                angle = (k / n) * 2 * math.pi
                x = cx + radius * math.cos(angle)
                y = cy + radius * math.sin(angle)
                coordinates.append((x, y))

            return coordinates

        def generate_linear_coordinates(n, left_boundary=0, right_boundary=5, fixed_y=0.5):
            coordinates = []
            for i in range(n):
                x = 0 + (i / (n - 1)) * (right_boundary - left_boundary)
                y = fixed_y
                coordinates.append((x, y))

            return coordinates

        uniform_dist = 5 / 24
        graph_width = uniform_dist * len(V)

        def generate_uniform_linear_coordinates(n, fixed_y=0.5, fixed_distance=uniform_dist):
            coordinates = []
            for i in range(n):
                x = 0 + i * fixed_distance
                y = fixed_y
                coordinates.append((x, y))

            return coordinates

        ## node parameters: same for the two graphs (nodes are always the same)
        V_pos = {}
        # cor = generate_circular_coordinates(len(V))
        cor = generate_uniform_linear_coordinates(len(V))
        for node_itr_ind in range(len(V)):
            V_pos[V[node_itr_ind]] = cor[node_itr_ind]

        V_colors = {}
        all_nodes = list(self.node_name.values())
        terminal_nodes = self.get_chr_start_end_nodes()
        for node_itr in all_nodes:
            if node_itr in terminal_nodes:
                V_colors[node_itr] = 'orange'
            else:
                V_colors[node_itr] = 'white'

        ## edges parameters: diff. for the two graphs
        karsim_E_weights = {}
        karsim_E_colors = nx.get_edge_attributes(G_karsim, 'color')
        for edge_itr in E_karsim_segment:
            karsim_E_weights[edge_itr] = E_karsim_segment[edge_itr]
        for edge_itr in E_karsim_transition:
            karsim_E_weights[edge_itr] = E_karsim_transition[edge_itr]

        omkar_E_weights = {}
        omkar_E_colors = nx.get_edge_attributes(G_omkar, 'color')
        for edge_itr in E_omkar_segment:
            omkar_E_weights[edge_itr] = E_omkar_segment[edge_itr]
        for edge_itr in E_omkar_transition:
            omkar_E_weights[edge_itr] = E_omkar_transition[edge_itr]

        ## plotting
        plt.figure(figsize=(graph_width * 4, 4))
        plot_karsim = ng.InteractiveGraph(G_karsim,
                                          node_color=V_colors,
                                          node_layout=V_pos,
                                          node_labels=True,
                                          edge_color=karsim_E_colors,
                                          edge_layout='arc',
                                          edge_labels=karsim_E_weights,
                                          arrows=True,
                                          node_size=6,
                                          node_label_offset=0.001,
                                          node_label_font_dict=dict(size=10),
                                          edge_label_fontdict=dict(size=9),
                                          scale=(graph_width, 1))
        plt.savefig(output_prefix + '.karsim.graph.png')

        plt.figure(figsize=(graph_width * 4, 4))
        plot_omkar = ng.InteractiveGraph(G_omkar,
                                         node_color=V_colors,
                                         node_layout=V_pos,
                                         node_labels=True,
                                         edge_color=omkar_E_colors,
                                         edge_layout='arc',
                                         edge_labels=omkar_E_weights,
                                         arrows=True,
                                         node_size=6,
                                         node_label_offset=0.001,
                                         node_label_font_dict=dict(size=10),
                                         edge_label_fontdict=dict(size=9),
                                         scale=(graph_width, 1))
        plt.savefig(output_prefix + '.omkar.graph.png')


def form_graph_from_cluster(cluster_file):
    index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file(cluster_file)
    graph = Graph()

    def iterative_add_segment_edge(path_list, target_graph):
        for path in path_list:
            for segment in path.linear_path.segments:
                graph.add_segment_edge(segment, target_graph)

    def iterative_add_transition_edge(path_list, target_graph):
        for path in path_list:
            for segment_ind in range(len(path.linear_path.segments) - 1):
                current_segment = path.linear_path.segments[segment_ind]
                next_segment = path.linear_path.segments[segment_ind + 1]
                graph.add_transition_edge(current_segment, next_segment, target_graph)

    # segment edge
    iterative_add_segment_edge(karsim_path_list, 'karsim')
    iterative_add_segment_edge(omkar_path_list, 'omkar')

    # transition edge
    iterative_add_transition_edge(karsim_path_list, 'karsim')
    iterative_add_transition_edge(omkar_path_list, 'omkar')

    return graph


def draw_graph(cluster_file):
    file_basename = cluster_file.split('/')[-1].split('.')[0]
    graph = form_graph_from_cluster(cluster_file)

    # graph.visualize_graph('new_data_files/complete_graphs/' + file_basename)
    print(graph.get_segment_distance())
    graph.prune_same_edges()
    # graph.visualize_graph('new_data_files/complete_graphs/' + file_basename + '.pruned')

    print(graph.get_segment_distance())

    graph.get_missed_transition_edge_count()
