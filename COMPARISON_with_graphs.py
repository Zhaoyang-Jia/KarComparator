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
    nodes: {(str, int): str}
    edges: Set[Edge]
    karsim_dict: {End_Point: [End_Point]}
    omkar_dict: {End_Point: [End_Point]}

    def __init__(self):
        self.nodes = {}
        self.edges = set()
        self.karsim_dict = {}
        self.omkar_dict = {}

    def add_edge_to_dict(self, edge: Edge, target_dict: str):
        if target_dict == 'omkar':
            target_dict = self.omkar_dict
        elif target_dict == 'karsim':
            target_dict = self.karsim_dict
        else:
            raise ValueError()
        if edge.endpoint1 in target_dict:
            target_dict[edge.endpoint1].append(edge.endpoint2)
        else:
            target_dict[edge.endpoint1] = [edge.endpoint2]

    def add_segment_edge(self, input_segment: Segment, target_graph: str):
        """
        :param input_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        new_endpoint1 = End_Point(input_segment.chr_name, input_segment.start, input_segment.kt_index[:-1] + "$")
        new_endpoint2 = End_Point(input_segment.chr_name, input_segment.end, input_segment.kt_index[:-1] + "#")
        new_edge = Edge(new_endpoint1, new_endpoint2, 'segment')

        # add node
        self.nodes[(new_endpoint1.origin_chr, new_endpoint1.position)] = new_endpoint1.name
        self.nodes[(new_endpoint2.origin_chr, new_endpoint2.position)] = new_endpoint2.name
        # add edge to set
        self.edges.add(new_edge)
        # add edge to dict
        self.add_edge_to_dict(new_edge, target_graph)

    def add_transition_edge(self, start_segment: Segment, end_segment: Segment, target_graph: str):
        """
        transition edge from the end of the start_segment to the start of the end_segment
        :param start_segment:
        :param end_segment:
        :param target_graph: "omkar" or "karsim"
        :return:
        """
        new_endpoint1 = End_Point(start_segment.chr_name, start_segment.end, "transition_start")
        new_endpoint2 = End_Point(end_segment.chr_name, end_segment.start, "transition_end")
        new_edge = Edge(new_endpoint1, new_endpoint2, 'transition')

        # add edge to set
        self.edges.add(new_edge)
        # add edge to dict
        self.add_edge_to_dict(new_edge, target_graph)

    def get_all_segment_edges(self):
        segment_edges = []
        for edge in self.edges:
            if edge.edge_type == 'segment':
                segment_edges.append(edge)
        return segment_edges

    def locate_edge_object(self, endpoint1, endpoint2):
        for edge in self.edges:
            if edge.endpoint1 == endpoint1 and endpoint2 == endpoint2:
                return edge

    def visualize_graph(self):
        # create sorted nodes (Endpoints)
        all_segment_edges = self.get_all_segment_edges()
        all_segment_endpoints = []
        for segment_edge in all_segment_edges:
            all_segment_endpoints.append(segment_edge.endpoint1)
            all_segment_endpoints.append(segment_edge.endpoint2)
        all_segment_endpoints = sorted(all_segment_endpoints)

        visualized_nodes = []
        for segment_endpoint in all_segment_endpoints:
            visualized_nodes.append(segment_endpoint.name)

        G = nx.DiGraph()
        G.add_nodes_from(visualized_nodes)

        # add all edges
        # TODO: edge weights
        for endpoint1 in self.karsim_dict:
            for endpoint2 in self.karsim_dict[endpoint1]:
                current_edge = self.locate_edge_object(endpoint1, endpoint2)
                if current_edge.edge_type == 'segment':
                    G.add_edge(endpoint1.name, endpoint2.name, color='red')
                elif current_edge.edge_type == 'transition':
                    endpoint1_name = self.nodes[(endpoint1.origin_chr, endpoint1.position)]
                    endpoint2_name = self.nodes[(endpoint2.origin_chr, endpoint2.position)]
                    G.add_edge(endpoint1_name, endpoint2_name, color='orange')
                else:
                    raise ValueError

        for endpoint1 in self.omkar_dict:
            for endpoint2 in self.omkar_dict[endpoint1]:
                current_edge = self.locate_edge_object(endpoint1, endpoint2)
                if current_edge.edge_type == 'segment':
                    G.add_edge(endpoint1.name, endpoint2.name, color='blue')
                elif current_edge.edge_type == 'transition':
                    endpoint1_name = self.nodes[(endpoint1.origin_chr, endpoint1.position)]
                    endpoint2_name = self.nodes[(endpoint2.origin_chr, endpoint2.position)]
                    G.add_edge(endpoint1_name, endpoint2_name, color='cyan')
                else:
                    raise ValueError

        pos = nx.circular_layout(G)  # TODO: replace with linear layout
        nx.draw(G, pos, with_labels=True)
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


