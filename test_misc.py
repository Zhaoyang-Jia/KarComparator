from scipy.optimize import linear_sum_assignment
import numpy as np
import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
from utils import *

# uniform_dist = 1 / 6
# graph_width = uniform_dist * 100
# width_multiplier = 2.5
# height_multiplier = 1
#
# G = nx.DiGraph()
# G.add_nodes_from(['0$', '0#'])
# G.add_edge('0$', '0#', weight=2)
#
# V_pos = {'0$': (0.0, 0.5), '0#': (1.0, 0.5)}
# E_weights = {('0$', '0#'): 2}
# E_pos = {('0$', '0#'): generate_parabola(0.0, 1.0, 0.5, 2)}
#
# plt.figure(figsize=(graph_width * width_multiplier, graph_width * height_multiplier))
# plot_karsim = ng.InteractiveGraph(G,
#                                   node_layout=V_pos,
#                                   node_labels=True,
#                                   edge_layout=E_pos,
#                                   edge_labels=E_weights,
#                                   arrows=True,
#                                   node_size=6,
#                                   node_label_offset=0.001,
#                                   node_label_font_dict=dict(size=10),
#                                   edge_label_fontdict=dict(size=9),
#                                   scale=(graph_width, 1))
# plt.savefig('test_merged_plot.png')

x = ""

if __name__ == "__main__":
    print(len(x))

