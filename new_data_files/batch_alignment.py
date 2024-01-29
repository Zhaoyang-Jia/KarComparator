from Structures import *
from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from read_cluster_file import *
from NW_aligner import *
from bipartite_matching import *

import os
import sys

for file in os.listdir('cluster_files/'):
    new_file_name = file.split('.')[0]
    new_file_name += '.alignment.txt'
    new_file_path = 'alignment_files/' + new_file_name
    with open(new_file_path, 'w') as fp_write:
        sys.stdout = fp_write
        index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file('cluster_files/' + file)
        hungarian_algorithm_for_cluster(karsim_path_list, omkar_path_list, index_to_segment, verbose=True)
        print()

sys.stdout = sys.__stdout__
