from read_OMKar_output import *
from KT_interpreter import *

import os

# paul_dremsek_data_dir = '/media/zhaoyang-new/workspace/paul_dremsek/omkar_output/'
paul_dremsek_data_dir = '/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/Bafna_Lab/Paul_Dremsek_OMKar_output/'
forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"


def batch_interpret(omkar_output_dir):
    files = [file for file in os.listdir(omkar_output_dir)]
    files = sorted(files, key=dremsek_file_keys)

    for file in files:
        file_path = omkar_output_dir + file
        print(file)
        mt_indexed_lists, mt_path_chrs, segment_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
        mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
        diploid_chr_report = interpret_chromosomal_anomalies(mt_path_chrs)
        print(diploid_chr_report)
        wt_path_dict = generate_wt_from_OMKar_output(segment_dict)
        wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
        interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, segment_size_dict)
        print()


def populate_wt_indexed_lists(mt_path_chrs, wt_path_dict):
    wt_indexed_lists = []
    for path_chr in mt_path_chrs:
        wt_indexed_lists.append(wt_path_dict[path_chr])
    return wt_indexed_lists


def dremsek_file_keys(f):
    return int(f.split('.')[0])


if __name__ == "__main__":
    batch_interpret(paul_dremsek_data_dir)
