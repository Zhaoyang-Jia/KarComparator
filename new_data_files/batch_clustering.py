from Structures import *
from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from read_cluster_file import *
from bipartite_matching import *

import os
import sys

karsim_folder = '/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/KarSimulator/'
omkar_folder = '/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/OMKar/'
forbidden_region_file = '/media/zhaoyang-new/workspace/KarSim/KarComparator/Metadata/acrocentric_telo_cen.bed'
output_dir = '/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/cluster_files/'

with open('clustering.log', 'w') as fp_write:
    sys.stdout = fp_write
    for file in os.listdir(omkar_folder):
        file_name = file.split('.')[0]
        omkar_file_path = omkar_folder + file
        karsim_file_path = karsim_folder + file_name + '.kt.txt'

        print(file)
        karsim_path_list = read_KarSimulator_output_to_path(karsim_file_path, forbidden_region_file)
        omkar_path_list = read_OMKar_output_to_path(omkar_file_path, forbidden_region_file)
        genome_wide_mutual_breaking(karsim_path_list, omkar_path_list)
        form_dependent_clusters(karsim_path_list, omkar_path_list, output_dir, prefix=file_name)
        print()

sys.stdout = sys.__stdout__