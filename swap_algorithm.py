class Edge:
    in_node: str
    out_node: str

    def __init__(self, in_node, out_node):
        self.in_node = in_node
        self.out_node = out_node


class Graph:
    edges: [Edge]

    def __init__(self):
        edges = []
        neighbors = {}
        reachable = {}


def form_graph(path_list):
    pass
