from Structures import *
from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from read_cluster_file import *
from NW_aligner import *
from bipartite_matching import *

import os
import sys


with open('alignment.txt', 'w') as fp_write:
    sys.stdout = fp_write
    for file in os.listdir('cluster_files/'):
        print(file.split('.')[0])
        index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file('cluster_files/' + file)
        hungarian_algorithm_for_cluster(karsim_path_list, omkar_path_list, index_to_segment, verbose=True)
        print()
sys.stdout = sys.__stdout__

