from scipy.optimize import linear_sum_assignment
import numpy as np
import networkx as nx
import netgraph as ng
import matplotlib.pyplot as plt
from utils import *
import os

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

# x = ''
# print(len(x))

import re


input_string0 = "Chr1: 237,668,707-248,943,333"
input_string1 = "Balanced translocation between Chr1 and Chr17, between segments Chr1: 237,668,707-248,943,333 (q43 - q44) and Chr17: 80,544,491-83,246,392 (q25.3 - q25.3)"
input_string2 = "Deletion on Chr2: 391-97,531,656 (q11.2 - q11.2)"

pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'

# Extract integers using regex
matches = re.finditer(pattern, input_string1)

for match_itr in matches:
    print(match_itr.group(1))
    print(match_itr.start())
    print(match_itr.end())
    print(input_string1[match_itr.start(): match_itr.end()])


# if __name__ == "__main__":
#     path = 'latex_reports/paul_dremsek_plots/002.pdf'
#     # path = '/media/zhaoyang-new/workspace/KarSim/KarComparator/latex_reports/sunnyside_plots/2.pdf'
#     if os.path.exists(path):
#         print('x')
#     else:
#         print('y')

