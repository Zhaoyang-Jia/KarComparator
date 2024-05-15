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

def extract_segments(input_string):
    pattern = r"Chr(\d+): (\d+),(\d+),(\d+)-(\d+),(\d+),(\d+) \((q\d+\.\d+) - (q\d+\.\d+)\)"
    matches = re.findall(pattern, input_string)
    segments = {}
    for match in matches:
        chr_num, start1, start2, start3, end1, end2, end3, q_start, q_end = match
        segment_key = f"Chr{chr_num}: {start1},{start2},{start3}-{end1},{end2},{end3} (q{q_start} - q{q_end})"
        segments[segment_key] = (f"chr{chr_num}", int(start1+start2+start3), int(end1+end2+end3))
    return segments

input_string = "Balanced translocation between Chr1 and Chr17, between segments Chr1: 237,668,707-248,943,333 (q43 - q44) and Chr17: 80,544,491-83,246,392 (q25.3 - q25.3)"
segments = extract_segments(input_string)
print(segments)

# if __name__ == "__main__":
#     path = 'latex_reports/paul_dremsek_plots/002.pdf'
#     # path = '/media/zhaoyang-new/workspace/KarSim/KarComparator/latex_reports/sunnyside_plots/2.pdf'
#     if os.path.exists(path):
#         print('x')
#     else:
#         print('y')

