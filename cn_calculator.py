from COMPARISON_with_graphs import *
from Structures import *
import numpy as np

#
# def read_cn_bin_file(cn_file_name):
#     cn_bins = []
#     with open(cn_file_name) as fp_read:
#         for line in fp_read:
#             line = line.replace('\n', '').split('\t')
#             cn_bins.append({'chrom': line[0],
#                             'start': int(line[1]),
#                             'end': int(line[2])})
#     return cn_bins
#
#
# def graph_assign_cn_bin(input_graph: Graph, input_cn_bins):
#     karsim_cn = np.array([0.0 for _ in input_cn_bins])
#     omkar_cn = np.array([0.0 for _ in input_cn_bins])
#     for start_node, end_nodes in input_graph.karsim_dict.items():
#         for end_node in end_nodes:
#             edge_type = end_node[2]
#             if edge_type == 'segment':
#                 c_seg = Segment(end_node[0], start_node[1], end_node[1])
#                 karsim_cn += c_seg.assign_cn_bin(input_cn_bins)
#     for start_node, end_nodes in input_graph.omkar_dict.items():
#         for end_node in end_nodes:
#             edge_type = end_node[2]
#             if edge_type == 'segment':
#                 c_seg = Segment(end_node[0], start_node[1], end_node[1])
#                 omkar_cn += c_seg.assign_cn_bin(input_cn_bins)
#     return karsim_cn, omkar_cn


# test_cluster_file =
# '/media/zhaoyang-new/workspace/KarSim/KarComparator/omkar_analyses_pipeline/builds/b14/cluster_files/23X_1q21_recurrent_microduplication_r1cluster_1.txt'
# graph = form_graph_from_cluster(test_cluster_file)
# print('x')
# cn_bins = read_cn_bin_file('Metadata/cn_bins_200kbp.txt')
#
# seg = Segment('Chr1', 8000, 220000)
# print(seg.assign_cn_bin(cn_bins))
