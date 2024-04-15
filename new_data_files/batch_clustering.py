from Structures import *
from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from read_cluster_file import *
from bipartite_matching import *
from utils import *

import os
import sys

karsim_folder = 'KarSimulator/'
omkar_folder = 'OMKar_testbuild4/'
omkar_output_folder = '../batch_processing/OMKar_testbuild4/'  # to access the dummy edge and change in CNV
forbidden_region_file = '../Metadata/acrocentric_telo_cen.bed'
output_dir = 'cluster_files_testbuild7/'

os.makedirs(output_dir, exist_ok=True)

with open('clustering.log', 'w') as fp_write:
    sys.stdout = fp_write
    for file in os.listdir(omkar_folder):
        file_name = file.split('.')[0]
        omkar_file_path = omkar_folder + file
        karsim_file_path = karsim_folder + file_name + '.kt.txt'
        omkar_case_folder = omkar_output_folder + '/' + file.split('/')[-1].replace('.txt', '')
        print(file)

        omkar_log_index_to_breakpoint_dict = omkar_log_get_node(omkar_case_folder)
        omkar_log_diff_edges = omkar_log_get_diff_edges(omkar_case_folder)
        omkar_log_diff_edges_coordinates = translate_indexed_edge_to_coordinates(omkar_log_diff_edges, omkar_log_index_to_breakpoint_dict)
        print(omkar_log_diff_edges_coordinates)

        karsim_segment_to_index_dict, karsim_path_list = read_KarSimulator_output_to_path(karsim_file_path, forbidden_region_file)
        omkar_index_to_segment_dict, omkar_path_list = read_OMKar_output_to_path(omkar_file_path, forbidden_region_file)
        omkar_segment_to_index_dict = reverse_dict(omkar_index_to_segment_dict)
        genome_wide_mutual_breaking(karsim_path_list, omkar_path_list)
        form_dependent_clusters(
            karsim_path_list,
            omkar_path_list,
            karsim_segment_to_index_dict,
            omkar_segment_to_index_dict,
            forbidden_region_file,
            output_dir,
            prefix=file_name,
            omkar_modified_edges=omkar_log_diff_edges_coordinates
        )
        print()

sys.stdout = sys.__stdout__
