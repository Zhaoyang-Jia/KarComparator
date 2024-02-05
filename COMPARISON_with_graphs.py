from typing import Set
import networkx as nx
import matplotlib.pyplot as plt

from Structures import *
from read_cluster_file import *


class End_Point:
    origin_chr: str
    position: int
    name: str

    def __init__(self, origin_chr: str, position: int, segment_name: str):
        self.origin_chr = origin_chr
        self.position = position
        self.name = segment_name

    def __hash__(self):
        return hash((self.origin_chr, self.position))

    def __str__(self):
        return "[{},{}:{}]".format(self.name, self.origin_chr, self.position)

    def __eq__(self, other):
        if (self.origin_chr, self.position) == (other.origin_chr, other.position):
            return True
        else:
            return False

    def __lt__(self, other):
        self_chr = self.origin_chr[3:]
        if self_chr == "X":
            self_chr = 23
        elif self_chr == "Y":
            self_chr = 24
        else:
            self_chr = int(self_chr)

        other_chr = other.origin_chr[3:]
        if other_chr == "X":
            other_chr = 23
        elif other_chr == "Y":
            other_chr = 24
        else:
            other_chr = int(other_chr)

        if self_chr != other_chr:
            return self_chr < other_chr
        else:
            return self.position < other.position


class Edge:
    endpoint1 = End_Point
    endpoint2 = End_Point
    edge_type = str

    def __init__(self, endpoint1: End_Point, endpoint2: End_Point, edge_type: str):
        self.endpoint1 = endpoint1
        self.endpoint2 = endpoint2
        self.edge_type = edge_type

    def __str__(self):
        return "({},{},{})".format(self.endpoint1, self.endpoint2, self.edge_type)

    def distance(self, other):
        if (self.endpoint1.origin_chr != other.endpoint1.origin_chr) or (self.endpoint2.origin_chr != other.endpoint2.origin_chr):
            return float('inf')
        else:
            endpoint1_dist = abs(self.endpoint1.position - other.endpoint1.position)
            endpoint2_dist = abs(self.endpoint2.position - other.endpoint2.position)
            return float(endpoint1_dist + endpoint2_dist)


class Graph:
    node_name: {(str, int): str}  # chr, position: node name
    edge_type: {(str, int, str, int): str}  # chr1, position1, chr2, position2: edge type
    karsim_dict: {(str, int): (str, int)}
    omkar_dict: {(str, int): (str, int)}

    def __init__(self):
        self.node_name = {}
        self.edge_type = {}
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

    def visualize_graph(self):
        # create sorted nodes (Endpoints)
        def custom_sort_node(node_tuple):
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

        nodes = self.node_name.keys()
        sorted_nodes = sorted(nodes, key=custom_sort_node)
        V = []
        for node in sorted_nodes:
            V.append(self.node_name[node])

        E_karsim_segment, E_karsim_transition = self.gather_edges('karsim')
        E_omkar_segment, E_omkar_transition = self.gather_edges('omkar')

        def iterative_add_edge(edges_dict: {(str, str): int}, group_color, graph):
            for edge, edge_weight in edges_dict.items():
                graph.add_edge(edge[0], edge[1], color=group_color, weight=edge_weight)

        G_karsim = nx.DiGraph()
        G_karsim.add_nodes_from(V)

        G_omkar = nx.DiGraph()
        G_omkar.add_nodes_from(V)

        # add all edges
        iterative_add_edge(E_karsim_segment, 'black', G_karsim)
        iterative_add_edge(E_karsim_transition, 'blue', G_karsim)
        iterative_add_edge(E_omkar_segment, 'black', G_omkar)
        iterative_add_edge(E_omkar_transition, 'blue', G_omkar)

        karsim_color_mapping = nx.get_edge_attributes(G_karsim, 'color')
        karsim_pos = nx.circular_layout(G_karsim, scale=2)
        karsim_weight_labels = nx.get_edge_attributes(G_karsim, 'weight')

        karsim_colors = []
        for edge_itr, color_itr in karsim_color_mapping.items():
            karsim_colors.append(color_itr)

        # node_sizes = [100 for _ in V]

        # pos = {}
        # y_value = 0
        # x_value = 0
        # for node in V:
        #     pos[node] = [x_value, y_value]
        #     x_value += 300

        # plt.figure(figsize=(100, 100))
        nx.draw(G_karsim, karsim_pos, edge_color=karsim_colors, with_labels=True)
        nx.draw_networkx_edge_labels(G_karsim, karsim_pos, edge_labels=karsim_weight_labels)
        plt.savefig("test_graph.png")


def draw_graph(cluster_file):
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

    graph.visualize_graph()

