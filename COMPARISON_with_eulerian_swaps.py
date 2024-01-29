from Structures import *
from read_cluster_file import *

from typing import Set
import re


class Node:
    segment: Segment
    out_node: {Segment: int}

    def __init__(self, segment):
        self.segment = segment
        self.out_node = {}
        self.reachable = set()

    def __str__(self):
        return_str = "str(self.segment) + \n"
        for segment in self.out_node:
            return_str += '\t'
            return_str += str(segment) + ": " + str(self.out_node[segment])
            return_str += '\n'
        return return_str

    def __hash__(self):
        """
        Node with the same Segment ident will be treated the same in a set
        :return:
        """
        return hash(self.segment)

    def add_edge(self, node2):
        if node2 in self.out_node:
            self.out_node[node2] += 1
        else:
            self.out_node[node2] = 1

    def is_junction_node(self):
        if len(self.out_node) >= 2:
            return True
        else:
            return False


class Graph:
    nodes: {Segment: Node}
    source_out_nodes: {Segment: int}
    sink_in_nodes: {Segment: int}
    junction_nodes: Set[Node]

    def __init__(self):
        self.nodes = {}
        self.source_out_nodes = {}
        self.sink_in_nodes = {}
        self.junction_nodes = set()

    def __str__(self):
        return_str = ""
        for key, value in self.nodes.items():
            return_str += str(value)
        return return_str

    def add_node(self, node):
        self.nodes[node.segment] = node

    def get_node(self, segment):
        return

    def update_junction_nodes(self):
        """
        update to include all the junction nodes in the Graph
        :return: number of junction nodes
        """
        count = 0
        for key, value in self.nodes.items():
            if value.is_junction_node():
                count += 1
                self.junction_nodes.add(value)
        return count

    def update_reachable_nodes(self):
        """
        for all junction nodes, update their reachable nodes
        :return:
        """
        # TODO: optimize
        pass

    def append_source(self, node2):
        if node2 in self.source_out_nodes:
            self.source_out_nodes[node2] += 1
        else:
            self.source_out_nodes[node2] = 1

    def append_sink(self, node1):
        if node1 in self.sink_in_nodes:
            self.sink_in_nodes[node1] += 1
        else:
            self.sink_in_nodes[node1] = 1


def form_graph(path_list: [Path]):
    graph = Graph()
    for path in path_list:
        # add the first node to source
        if path.linear_path.segments[0] not in graph.nodes:
            first_node = Node(path.linear_path.segments[0])
            graph.add_node(first_node)
        else:
            first_node = graph.nodes[path.linear_path.segments[0]]
        graph.append_source(first_node)

        # for the rest of the graph
        left_node = first_node  # start with the first node as the left node
        for segment_ind in range(len(path.linear_path.segments) - 1):
            right_node_segment = path.linear_path.segments[segment_ind + 1]
            # left_node is always already present in the Graph
            if right_node_segment not in graph.nodes:
                right_node = Node(right_node_segment)
                graph.add_node(right_node)
            else:
                right_node = graph.nodes[right_node_segment]
            left_node.add_edge(right_node)

            left_node = right_node

        # add the last node to sink
        last_node = left_node
        graph.append_sink(last_node)

    return Graph

def read_alignment_file(alignment_file, cluster_file):
    """
    combine alignment file's bipartite matching + scoring with the paths from the cluster_file
    :param alignment_file:
    :param cluster_file:
    :return:
    """
    total_cost = -1
    path_cost = []  # associated costs for each path in the path_list

    # index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file

    with open(alignment_file) as fp_read:
        line1 = fp_read.readline()
        line1 = line1.replace('\n', '').split(': ')[1]
        total_cost = int(line1)

        karsim_path_name = ''
        omkar_path_name = ''
        for line in fp_read:
            if line.startswith('alignment'):
                info = re.match(r"alignment (\d+): (.+), (.+) \t cost: (\d+)", line)
                print(info)