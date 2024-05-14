from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from utils import *

import os
import sys
import argparse

karsim_folder = '../new_data_files/KarSimulator/'
forbidden_region_file = '../Metadata/acrocentric_telo_cen.bed'

parser = argparse.ArgumentParser()
parser.add_argument('--omkar_output', type=str, help='path to DIR of OMKar outputs')
parser.add_argument('--omkar_paths_output', type=str, help="path to DIR containing the OMKar path files")
parser.add_argument('--output_dir', type=str, help= "path to this program's output")
args = parser.parse_args()

omkar_output_folder = args.omkar_output
omkar_folder = args.omkar_paths_output
output_dir = args.output_dir
# omkar_folder = 'OMKar_testbuild11/'
# omkar_output_folder = '../batch_processing/OMKar_testbuild11/'  # to access the dummy edge and change in CNV
# output_dir = 'cluster_files_testbuild14/'

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
