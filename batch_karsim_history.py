import os
from read_KarSimulator_output import *


karsim_folder = 'new_data_files/KarSimulator/'
karsim_history_output_folder = 'packaged_data/Karsimulator_history_intermediate/'
forbidden_region_file = 'Metadata/acrocentric_telo_cen.bed'

for file in os.listdir(karsim_folder):
    file_name = file.split('.')[0]
    karsim_file_path = karsim_folder + file
    index_dict, path_list, event_histories = read_KarSimulator_output(karsim_file_path, forbidden_region_file)
    sv_edge_list = label_event_sv_edge(index_dict, path_list, event_histories)
    history_intermediate_file = karsim_history_output_folder + file_name + '.history_sv.txt'
    generate_history_SV_edge_labels(event_histories, sv_edge_list, history_intermediate_file)
