import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
import math
from collections import Counter
from scipy.optimize import linear_sum_assignment
import numpy as np

from Structures import *
from read_cluster_file import *
from utils import *
from forbidden_region_processing import read_forbidden_regions

d = 500000


def intra_transition_edge_distance(chr1, pos1, chr2, pos2):
    if chr1 != chr2:
        return -1
    else:
        return abs(pos2 - pos1) - 1


def inter_transition_edge_distance(edge1_chr1, edge1_pos1, edge1_chr2, edge1_pos2,
                                   edge2_chr1, edge2_pos1, edge2_chr2, edge2_pos2):
    if edge1_chr1 != edge2_chr1 or edge1_chr2 != edge2_chr2:
        return -1
    return abs(edge1_pos1 - edge2_pos1) + abs(edge1_pos2 - edge2_pos2)


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


def pop_edge(chr1, pos1, chr2, pos2, edge_type, target_dict):
    node2_list = target_dict[(chr1, pos1)]
    if (chr2, pos2, edge_type) not in node2_list:
        raise ValueError('edge does not exist')
    else:
        node2_list.remove((chr2, pos2, edge_type))


class Graph:
    node_name: {(str, int): str}  # chr, position: node name
    segment_edge_type: {(str, int, str, int): str}  # documents for comparison forbidden regions
    karsim_dict: {(str, int): [(str, int, str)]}  # chr1, pos1: chr2, pos2, edge_type
    omkar_dict: {(str, int): [(str, int, str)]}  # chr1, pos1: chr2, pos2, edge_type
    approximated_cnv: int
    karsim_n_transition_approximated: int
    omkar_n_transition_approximated: int
    edges_of_interest: [(str, str, int, str, int)]  # start_node, end_node, multiplicity, edge_type, edge_distance
    events: {str: int}  # event_name: multiplicity

    def __init__(self):
        self.node_name = {}
        self.segment_edge_type = {}
        self.karsim_dict = {}
        self.omkar_dict = {}
        self.approximated_cnv = 0
        self.karsim_n_transition_approximated = 0
        self.omkar_n_transition_approximated = 0
        self.edges_of_interest = []
        self.events = {}

    def add_edge_to_dict(self, chr1, pos1, chr2, pos2, edge_type, target_dict: str):
        if target_dict == 'omkar':
            target_dict = self.omkar_dict
        elif target_dict == 'karsim':
            target_dict = self.karsim_dict
        else:
            raise ValueError()

        if (chr1, pos1) in target_dict:
            target_dict[(chr1, pos1)].append((chr2, pos2, edge_type))
        else:
            target_dict[(chr1, pos1)] = [(chr2, pos2, edge_type)]

    def add_segment_edge(self, input_segment: Segment, target_graph: str):
        """
        :param input_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        # add node
        if input_segment.direction():
            self.node_name[(input_segment.chr_name, input_segment.start)] = input_segment.kt_index[:-1] + "s"
            self.node_name[(input_segment.chr_name, input_segment.end)] = input_segment.kt_index[:-1] + "t"
        else:
            self.node_name[(input_segment.chr_name, input_segment.end)] = input_segment.kt_index[:-1] + "s"
            self.node_name[(input_segment.chr_name, input_segment.start)] = input_segment.kt_index[:-1] + "t"

        # document edge type
        self.segment_edge_type[(input_segment.chr_name, input_segment.start,
                                input_segment.chr_name, input_segment.end)] = input_segment.segment_type

        # add edge to dict
        self.add_edge_to_dict(input_segment.chr_name, input_segment.start,
                              input_segment.chr_name, input_segment.end,
                              'segment', target_graph)

    def add_transition_edge(self, start_segment: Segment, end_segment: Segment, target_graph: str):
        """
        transition edge from the end of the start_segment to the start of the end_segment
        :param start_segment:
        :param end_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        # add edge to dict
        self.add_edge_to_dict(start_segment.chr_name, start_segment.end,
                              end_segment.chr_name, end_segment.start,
                              'transition', target_graph)

    # def add_start_end_transition_edges(self):
    #     terminal_nodes = self.get_chr_start_end_nodes()
    #     node_name_to_coordinates = reverse_dict(self.node_name)
    #     terminal_node_coordinates = {}
    #     for node_name in terminal_nodes:
    #         node_coordinate = node_name_to_coordinates[node_name]
    #         node_chr = node_coordinate[0]
    #         node_index = node_coordinate[1]
    #         if node_chr in terminal_node_coordinates:
    #             terminal_node_coordinates[node_chr].append(node_index)
    #         else:
    #             terminal_node_coordinates[node_chr] = [node_index]
    #
    #     # print('x')
    #
    #     for start_node in self.omkar_start_node:
    #         start_node_chr = start_node[0]
    #         if start_node_chr not in terminal_node_coordinates:
    #             continue
    #         graph_start_node_coordinate = min(terminal_node_coordinates[start_node_chr])  # contains 2 coor, min is the start
    #         if start_node == (start_node_chr, graph_start_node_coordinate):
    #             continue
    #
    #         if (start_node_chr, graph_start_node_coordinate) in self.omkar_dict:
    #             self.omkar_dict[(start_node_chr, graph_start_node_coordinate)].append((*start_node, 'transition'))
    #         else:
    #             self.omkar_dict[(start_node_chr, graph_start_node_coordinate)] = [(*start_node, 'transition')]
    #
    #     for end_node in self.omkar_end_node:
    #         end_node_chr = end_node[0]
    #         if end_node_chr not in terminal_node_coordinates:
    #             # FIXME: why would this even occur?
    #             continue
    #         graph_end_node_coordinate = max(terminal_node_coordinates[end_node_chr])  # contains 2 coor, max is the end
    #         if end_node == (end_node_chr, graph_end_node_coordinate):
    #             continue
    #
    #         if end_node in self.omkar_dict:
    #             self.omkar_dict[end_node].append((end_node_chr, graph_end_node_coordinate, 'transition'))
    #         else:
    #             self.omkar_dict[end_node] = [(end_node_chr, graph_end_node_coordinate, 'transition')]

    def prune_same_edges(self):
        dict1 = self.karsim_dict
        dict2 = self.omkar_dict

        # Find the common keys between the dictionaries
        common_keys = dict1.keys() & dict2.keys()

        # Iterate over common keys and remove instances
        for key in common_keys:
            # TODO: test if this will distinguish between different edge types of the same chr, pos
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
        """
        :param target_graph: 'karsim' or 'omkar'
        :return:
        """
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
                edge_type = node2[2]

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

    def get_segment_distance(self):
        """
        return the total distance of the prunned graphs' remaining segment edges
        :return:
        """
        # reverse the dict
        node_name_to_node_position = reverse_dict(self.node_name)

        total_distance = 0
        for node1, value in self.karsim_dict.items():
            for node2 in value:
                if node2[2] == 'transition':
                    continue
                # all segment edge are between the same chromosome nodes
                edge_type = self.segment_edge_type[(node1[0], node1[1], node2[0], node2[1])]
                # total_distance += abs(node2[1] - node1[1] + 1)
                if edge_type.startswith('telomere') or edge_type.startswith('acrocentric'):
                    continue
                else:
                    total_distance += abs(node2[1] - node1[1] + 1)

        for node1, value in self.omkar_dict.items():
            for node2 in value:
                if node2[2] == 'transition':
                    continue
                # all segment edge are between the same chromosome nodes
                edge_type = self.segment_edge_type[(node1[0], node1[1], node2[0], node2[1])]

                # total_distance += abs(node2[1] - node1[1] + 1)
                # TODO: recover this
                if edge_type.startswith('telomere') or edge_type.startswith('acrocentric'):
                    continue
                else:
                    total_distance += abs(node2[1] - node1[1] + 1)

        return total_distance

    def get_missed_transition_edges(self):
        """
        archived
        :return:
        """
        karsim_transition_edges = []
        for node1, value in self.karsim_dict.items():
            for node2 in value:
                if node2[2] == 'transition':
                    # skip reference edge
                    if node1[0] == node2[0] and node1[1] + 1 == node2[1]:
                        continue
                    else:
                        karsim_transition_edges.append((node1[0], node1[1], node2[0], node2[1]))

        omkar_transition_edges = []
        for node1, value in self.omkar_dict.items():
            for node2 in value:
                if node2[2] == 'transition':
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
        sorted_nodes = sorted(nodes, key=custom_sort_node)  # sorted by chr, then by positions -> chr always start to end sorted in the list
        terminal_nodes = [self.node_name[sorted_nodes[0]]]  # first node always the start of a chr
        node_ind = 1
        while node_ind <= len(sorted_nodes) - 3:
            current_node = sorted_nodes[node_ind]
            next_node = sorted_nodes[node_ind + 1]
            if current_node[0] != next_node[0]:
                terminal_nodes.append(self.node_name[current_node])
                terminal_nodes.append(self.node_name[next_node])
            node_ind += 2  # because of the st paired structure, only check the t against the next s
        terminal_nodes.append(self.node_name[sorted_nodes[-1]])  # last node always the end of a chr

        return terminal_nodes

    def remove_approximate_transition_edges(self):
        def mark_approximated_transition_edges(target_graph):
            if target_graph == 'karsim':
                target_dict = self.karsim_dict
            elif target_graph == 'omkar':
                target_dict = self.omkar_dict
            else:
                raise ValueError

            # mark all SMALL edges for removal
            edges_to_remove = []
            for node1, value in target_dict.items():
                for node2 in value:
                    chr1 = node1[0]
                    pos1 = node1[1]
                    chr2 = node2[0]
                    pos2 = node2[1]
                    edge_type = node2[2]

                    if edge_type == 'segment':
                        continue

                    intra_distance = intra_transition_edge_distance(chr1, pos1, chr2, pos2)
                    if intra_distance == -1:
                        # -1 encodes +inf
                        continue
                    elif intra_distance < d:
                        edges_to_remove.append((chr1, pos1, chr2, pos2, edge_type))
            return edges_to_remove

        # removal
        karsim_edges_to_remove = mark_approximated_transition_edges('karsim')
        omkar_edges_to_remove = mark_approximated_transition_edges('omkar')
        for current_edge in karsim_edges_to_remove:
            distance_param = current_edge[:4]
            current_distance = intra_transition_edge_distance(*distance_param)
            self.approximated_cnv += current_distance
            if current_distance > 0:
                self.karsim_n_transition_approximated += 1
            pop_edge(*current_edge, self.karsim_dict)
        for current_edge in omkar_edges_to_remove:
            distance_param = current_edge[:4]
            current_distance = intra_transition_edge_distance(*distance_param)
            self.approximated_cnv += current_distance
            if current_distance > 0:
                self.omkar_n_transition_approximated += 1
            pop_edge(*current_edge, self.omkar_dict)

    def match_transition_edges(self):
        large_value = 10000000
        node_name_to_coordinate = reverse_dict(self.node_name)

        ## gather all transition edges
        _, karsim_transition_edge_dict = self.gather_edges('karsim')
        _, omkar_transition_edge_dict = self.gather_edges('omkar')
        karsim_transition_edges = []
        omkar_transition_edges = []

        for edge, multiplicity in karsim_transition_edge_dict.items():
            node1 = node_name_to_coordinate[edge[0]]
            node2 = node_name_to_coordinate[edge[1]]
            chr1 = node1[0]
            pos1 = node1[1]
            chr2 = node2[0]
            pos2 = node2[1]
            for itr in range(multiplicity):
                karsim_transition_edges.append((chr1, pos1, chr2, pos2))

        for edge, multiplicity in omkar_transition_edge_dict.items():
            node1 = node_name_to_coordinate[edge[0]]
            node2 = node_name_to_coordinate[edge[1]]
            chr1 = node1[0]
            pos1 = node1[1]
            chr2 = node2[0]
            pos2 = node2[1]
            for itr in range(multiplicity):
                omkar_transition_edges.append((chr1, pos1, chr2, pos2))

        if len(karsim_transition_edges) == 0 or len(omkar_transition_edges) == 0:
            return

        ## add dummies
        # len_diff = len(karsim_transition_edges) - len(omkar_transition_edges)
        # if len_diff > 0:
        #     for itr in range(len_diff):
        #         omkar_transition_edges.append('dummy')
        # elif len_diff < 0:
        #     for itr in range(abs(len_diff)):
        #         karsim_transition_edges.append('dummy')

        ## populate cost_matrix
        cost_matrix = [[large_value for col in range(len(omkar_transition_edges))] for row in range(len(karsim_transition_edges))]
        for row_index in range(len(karsim_transition_edges)):
            for col_index in range(len(omkar_transition_edges)):
                karsim_edge = karsim_transition_edges[row_index]
                omkar_edge = omkar_transition_edges[col_index]

                # if the start/end direction do not match, transition edges are different
                karsim_node1_sign = self.node_name[(karsim_edge[0], karsim_edge[1])][-1]
                karsim_node2_sign = self.node_name[(karsim_edge[2], karsim_edge[3])][-1]
                omkar_node1_sign = self.node_name[(omkar_edge[0], omkar_edge[1])][-1]
                omkar_node2_sign = self.node_name[(omkar_edge[2], omkar_edge[3])][-1]
                if karsim_node1_sign != omkar_node1_sign or karsim_node2_sign != omkar_node2_sign:
                    continue

                distance = inter_transition_edge_distance(*karsim_edge, *omkar_edge)
                if distance == -1:
                    # -1 signals no match
                    continue

                if distance < 2 * d:
                    # only update cost if within threshold
                    cost_matrix[row_index][col_index] = distance

        ## bipartite matching
        np_cost_matrix = np.array(cost_matrix)
        karsim_assignment, omkar_assignment = linear_sum_assignment(np_cost_matrix)

        ## remove all approximate matching
        for ind in range(len(karsim_assignment)):
            karsim_ind = karsim_assignment[ind]
            omkar_ind = omkar_assignment[ind]
            current_distance = cost_matrix[karsim_ind][omkar_ind]

            if current_distance < large_value:
                # an approximate matching is found, pop both edges in the matching
                print("matched: " + str(current_distance))
                self.approximated_cnv += current_distance
                karsim_edge = karsim_transition_edges[karsim_ind]
                omkar_edge = omkar_transition_edges[omkar_ind]
                self.karsim_dict[(karsim_edge[0], karsim_edge[1])].remove((karsim_edge[2], karsim_edge[3], 'transition'))
                self.omkar_dict[(omkar_edge[0], omkar_edge[1])].remove((omkar_edge[2], omkar_edge[3], 'transition'))
                self.karsim_n_transition_approximated += 1
                self.omkar_n_transition_approximated += 1

    def remove_forbidden_nodes(self, forbidden_region_file):
        """
        remove all nodes associated with forbidden regions, and then remove all edges associated with these nodes
        :param forbidden_region_file:
        :return:
        """
        def remove_node_and_associated_edge(input_node):
            self.node_name.pop(input_node)
            # remove out going edges
            if input_node in self.karsim_dict:
                self.karsim_dict.pop(input_node)
            if input_node in self.omkar_dict:
                self.omkar_dict.pop(input_node)
            # remove incoming edges
            for node1, node2_list in self.karsim_dict.items():
                node2_to_remove = []
                for idx, node2 in enumerate(node2_list):
                    if node2[0] == input_node[0] and node2[1] == input_node[1]:
                        node2_to_remove.append(idx)
                new_node2_list = [value for i, value in enumerate(node2_list) if i not in node2_to_remove]
                self.karsim_dict[node1] = new_node2_list
            for node1, node2_list in self.omkar_dict.items():
                node2_to_remove = []
                for idx, node2 in enumerate(node2_list):
                    if node2[0] == input_node[0] and node2[1] == input_node[1]:
                        node2_to_remove.append(idx)
                new_node2_list = [value for i, value in enumerate(node2_list) if i not in node2_to_remove]
                self.omkar_dict[node1] = new_node2_list

        forbidden_segments = read_forbidden_regions(forbidden_region_file).segments
        for segment in forbidden_segments:
            start_node = (segment.chr_name, segment.start)
            end_node = (segment.chr_name, segment.end)
            if start_node in self.node_name:
                remove_node_and_associated_edge(start_node)
            if end_node in self.node_name:
                remove_node_and_associated_edge(end_node)

    def visualize_graph(self, output_prefix, merged=False):
        # create sorted nodes (Endpoints)

        nodes = self.node_name.keys()
        sorted_nodes = sorted(nodes, key=custom_sort_node)
        V = []
        for node in sorted_nodes:
            V.append(self.node_name[node])

        E_karsim_segment, E_karsim_transition = self.gather_edges('karsim')
        E_omkar_segment, E_omkar_transition = self.gather_edges('omkar')

        def translate_segment_edge_type_to_edge_name():
            new_dict = {}
            for edge, edge_type in self.segment_edge_type.items():
                if (edge[0], edge[1]) in self.node_name and (edge[2], edge[3]) in self.node_name:
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

                if graph.has_edge(edge[0], edge[1]):
                    # FIXME: right now, if different types of the same edge exist, the multiplicity of the latters will be added to the first added type
                    graph[edge[0]][edge[1]]['weight'] += edge_weight
                else:
                    graph.add_edge(edge[0], edge[1], color=edge_color, weight=edge_weight)

        def iterative_add_edge_trace(edges_dict: {(str, str): int},
                                     trace_dict: {(str, str): (float, float)},
                                     peak_height_nonadjacent,
                                     peak_height_adjacent,
                                     adjacent_x_dist):
            """
            :param edges_dict:
            :param trace_dict:
            :param peak_height_nonadjacent: for nodes not right next to each other
            :param peak_height_adjacent: for nodes right next to each other
            :param adjacent_x_dist: nodes x-dist if they are adjacent
            :return:
            """
            for edge in edges_dict:
                node1_x = V_pos[edge[0]][0]
                node2_x = V_pos[edge[1]][0]

                if edge in trace_dict:
                    # FIXME: look into why this happens
                    continue
                    raise RuntimeError('edge already exists, multi-graph not supported')

                if edge[0] == edge[1]:
                    # self edge
                    trace_dict[edge] = generate_circle(node1_x, peak_height_nonadjacent, 7)
                else:
                    # regular edge
                    if abs(node2_x - node1_x) <= adjacent_x_dist + 0.001:
                        peak_x = round((node1_x + node2_x) / 2, 6)
                        trace_dict[edge] = generate_parabola(node1_x, node2_x, peak_x, peak_height_adjacent)
                    else:
                        peak_x = round((node1_x + node2_x) / 2, 6)
                        trace_dict[edge] = generate_parabola(node1_x, node2_x, peak_x, peak_height_nonadjacent)

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

        uniform_dist = 1 / 6
        graph_width = uniform_dist * 200
        width_multiplier = 2.5
        height_multiplier = 1

        def generate_uniform_linear_coordinates(n, fixed_y=0.5, fixed_distance=uniform_dist):
            coordinates = []
            for i in range(n):
                x = round(0 + i * fixed_distance, 4)
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
        if not merged:
            karsim_E_pos = {}
            iterative_add_edge_trace(E_karsim_segment, karsim_E_pos, 1, 0.58, uniform_dist)
            iterative_add_edge_trace(E_karsim_transition, karsim_E_pos, 1, 0.58, uniform_dist)

            omkar_E_pos = {}
            iterative_add_edge_trace(E_omkar_segment, omkar_E_pos, 0, 0.42, uniform_dist)
            iterative_add_edge_trace(E_omkar_transition, omkar_E_pos, 0, 0.42, uniform_dist)

            plt.figure(figsize=(graph_width * width_multiplier, graph_width * height_multiplier))
            plot_karsim = ng.InteractiveGraph(G_karsim,
                                              node_color=V_colors,
                                              node_layout=V_pos,
                                              node_labels=True,
                                              edge_color=karsim_E_colors,
                                              edge_layout=karsim_E_pos,
                                              edge_labels=karsim_E_weights,
                                              arrows=True,
                                              node_size=6,
                                              node_label_offset=0.001,
                                              node_label_font_dict=dict(size=20),
                                              edge_label_fontdict=dict(size=15),
                                              scale=(graph_width, 1))
            plt.savefig(output_prefix + '.karsim.graph.png')

            plt.figure(figsize=(graph_width * width_multiplier, graph_width * height_multiplier))
            plot_omkar = ng.InteractiveGraph(G_omkar,
                                             node_color=V_colors,
                                             node_layout=V_pos,
                                             node_labels=True,
                                             edge_color=omkar_E_colors,
                                             edge_layout=omkar_E_pos,
                                             edge_labels=omkar_E_weights,
                                             arrows=True,
                                             node_size=6,
                                             node_label_offset=0.001,
                                             node_label_font_dict=dict(size=20),
                                             edge_label_fontdict=dict(size=15),
                                             scale=(graph_width, 1))
            plt.savefig(output_prefix + '.omkar.graph.png')
        else:
            G_merged = nx.DiGraph()
            G_merged.add_nodes_from(V)

            iterative_add_edge(E_karsim_segment, 'black', G_merged, forbidden_segment_edge_labels=True)
            iterative_add_edge(E_karsim_transition, 'red', G_merged)
            iterative_add_edge(E_omkar_segment, 'black', G_merged, forbidden_segment_edge_labels=True)
            iterative_add_edge(E_omkar_transition, 'red', G_merged)

            merged_E_pos = {}
            iterative_add_edge_trace(E_karsim_segment, merged_E_pos, 1, 0.58, uniform_dist)
            iterative_add_edge_trace(E_karsim_transition, merged_E_pos, 1, 0.58, uniform_dist)
            iterative_add_edge_trace(E_omkar_segment, merged_E_pos, 0, 0.42, uniform_dist)
            iterative_add_edge_trace(E_omkar_transition, merged_E_pos, 0, 0.42, uniform_dist)

            merged_E_colors = {**karsim_E_colors, **omkar_E_colors}
            merged_E_weights = {**karsim_E_weights, **omkar_E_weights}
            plt.figure(figsize=(graph_width * width_multiplier, graph_width * height_multiplier * 2))  # merged graph is going to be twice as high
            plot_merged = ng.InteractiveGraph(G_merged,
                                              node_color=V_colors,
                                              node_layout=V_pos,
                                              node_labels=True,
                                              edge_color=merged_E_colors,
                                              edge_layout=merged_E_pos,
                                              edge_labels=merged_E_weights,
                                              arrows=True,
                                              node_size=6,
                                              node_label_offset=0.001,
                                              node_label_font_dict=dict(size=10),
                                              edge_label_fontdict=dict(size=9),
                                              scale=(graph_width, 1))
            plt.savefig(output_prefix + '.merged.graph.png')


def form_graph_from_cluster(cluster_file, forbidden_region_file='/media/zhaoyang-new/workspace/KarSim/KarComparator/Metadata/acrocentric_telo_cen.bed'):
    index_to_segment, karsim_path_list, omkar_path_list, labeled_edges = read_cluster_file(cluster_file)
    graph = Graph()
    graph.edges_of_interest = labeled_edges

    def iterative_add_segment_edge(path_list, target_graph):
        for path in path_list:
            for segment in path.linear_path.segments:
                graph.add_segment_edge(segment, target_graph)

    def iterative_add_transition_edge(path_list, target_graph):
        for path in path_list:
            for segment_ind in range(len(path.linear_path.segments) - 1):
                current_segment = path.linear_path.segments[segment_ind]
                next_segment = path.linear_path.segments[segment_ind + 1]
                if current_segment.end == next_segment.start:
                    # FIXME: minor bug in the segment indexing exists in .kt file
                    continue
                graph.add_transition_edge(current_segment, next_segment, target_graph)

                # if segment_ind == 0:
                #     if target_graph == 'omkar':
                #         graph.omkar_start_node.append((current_segment.chr_name, current_segment.start))
                # elif segment_ind == len(path.linear_path.segments) - 2:
                #     if target_graph == 'omkar':
                #         graph.omkar_end_node.append((next_segment.chr_name, next_segment.end))

    # segment edge
    iterative_add_segment_edge(karsim_path_list, 'karsim')
    iterative_add_segment_edge(omkar_path_list, 'omkar')

    # transition edge
    iterative_add_transition_edge(karsim_path_list, 'karsim')
    iterative_add_transition_edge(omkar_path_list, 'omkar')

    # graph.add_start_end_transition_edges()

    return graph


def draw_graph(cluster_file, output_dir):
    import os
    import shutil

    file_basename = cluster_file.split('/')[-1].split('.')[0]
    file_basename_no_cluster = file_basename.split('cluster')[0]
    folder = output_dir + file_basename + '/'

    os.makedirs(folder, exist_ok=True)

    shutil.copyfile('new_data_files/OMKar/' + file_basename_no_cluster + '.1.txt',
                    folder + file_basename_no_cluster + '.omkar_paths.txt')
    shutil.copyfile('new_data_files/KarSimulator/' + file_basename_no_cluster + '.kt.txt',
                    folder + file_basename_no_cluster + '.karsim_paths.txt')
    # shutil.copyfile('new_data_files/alignment_files/' + file_basename + '.alignment.txt',
    #                 folder + file_basename + '.alignment.txt')
    shutil.copyfile(cluster_dir + file_basename + '.txt',
                    folder + file_basename + '.cluster.txt')

    base_file_name = cluster_file.split('/')[-1].split('cluster')[0]
    node_file = omkar_output_dir + base_file_name + '.1/' + base_file_name + '.1.preILP_nodes.txt'
    edge_file = omkar_output_dir + base_file_name + '.1/' + base_file_name + '.1.preILP_edges.txt'
    shutil.copyfile(node_file, folder + file_basename + '.preILP_nodes.txt')
    shutil.copyfile(edge_file, folder + file_basename + '.preILP_edges.txt')
    print(node_file)
    print(folder)
    print(file_basename)

    graph = form_graph_from_cluster(cluster_file)

    graph.visualize_graph(folder + 'raw')
    print('initial segment distance: ' + str(graph.get_segment_distance()))
    graph.prune_same_edges()
    graph.visualize_graph(folder + 'pruned', merged=True)
    print('pruned segment distance: ' + str(graph.get_segment_distance()))
    graph.remove_approximate_transition_edges()
    graph.visualize_graph(folder + 'approximated', merged=True)
    print('approximated segment distance: ' + str(graph.get_segment_distance()))
    print('approximated_cnv: ' + str(graph.approximated_cnv))
    graph.match_transition_edges()
    graph.visualize_graph(folder + 'matched', merged=True)
    print('post-matching segment distance: ' + str(graph.get_segment_distance()))
    print('approximated_cnv: ' + str(graph.approximated_cnv))

    # graph.get_missed_transition_edge_count()


def export_graph_vertices(cluster_file):
    graph = form_graph_from_cluster(cluster_file)
    return graph.node_name


omkar_output_dir = 'batch_processing/omkar_output_temp/'
cluster_dir = 'new_data_files/cluster_files_testbuild6/'

if __name__ == "__main__":
    import subprocess
    file_name = '23Y_1q21_recurrent_microdeletion_r2'
    cluster_number = '5'

    input_cluster_file = 'new_data_files/cluster_files_testbuild6/' + file_name + 'cluster_' + cluster_number + '.txt'
    draw_graph(input_cluster_file, 'new_data_files/complete_graphs/')

    print('debug-omkar')
    chr_of_int = ['10']  # figure this out by looking at the cluster file first, required for running the debug_omkar
    debug_omkar_script = './debug_omkar.py'
    omkar_V_rename_dict = export_graph_vertices(input_cluster_file)
    input_file_basename = input_cluster_file.split('/')[-1].split('.')[0]
    subprocess.run(['python', debug_omkar_script,
                    '--output_dir', 'new_data_files/complete_graphs/' + input_file_basename + '/',
                    '--case_file', file_name,
                    '--chr_of_int', str(chr_of_int),
                    '--rename_dict', str(omkar_V_rename_dict)])


