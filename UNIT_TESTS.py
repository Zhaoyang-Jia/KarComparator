from Structures import *
from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from read_cluster_file import *
from NW_aligner import *


def test_read_KarSimulator_output_to_path():
    karsim_file = 'sample_input/12q14_microdeletion_v2_r2.kt.txt'
    forbidden_region_file = 'Metadata/acrocentric_telo_cen.bed'
    karsim_path_list = read_KarSimulator_output_to_path(karsim_file, forbidden_region_file)
    for path in karsim_path_list:
        print(path)


def test_read_OMKar_output_to_path():
    omkar_file = 'sample_input/12q14_microdeletion_v2_r2.1.txt'
    forbidden_region_file = 'Metadata/acrocentric_telo_cen.bed'
    omkar_path_list = read_OMKar_output_to_path(omkar_file, forbidden_region_file)
    for path in omkar_path_list:
        print(path)
    report_centromere_anomaly(omkar_path_list)


def test_form_dependent_clusters():
    karsim_file = 'sample_input/12q14_microdeletion_v2_r2.kt.txt'
    omkar_file = 'sample_input/12q14_microdeletion_v2_r2.1.txt'
    forbidden_region_file = 'Metadata/acrocentric_telo_cen.bed'
    output_dir = 'tmp/'

    karsim_path_list = read_KarSimulator_output_to_path(karsim_file, forbidden_region_file)
    omkar_path_list = read_OMKar_output_to_path(omkar_file, forbidden_region_file)
    genome_wide_mutual_breaking(karsim_path_list, omkar_path_list)
    form_dependent_clusters(karsim_path_list, omkar_path_list, output_dir)


def test_read_cluster_file():
    file = 'tmp/cluster_6.txt'
    karsim_path_list, omkar_path_list = read_cluster_file(file)
    for path in karsim_path_list:
        print(path)
    for path in omkar_path_list:
        print(path)


def test_NW_aligner():
    file = 'tmp/cluster_6.txt'
    index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file(file)
    path1 = karsim_path_list[3].linear_path.segments
    path2 = omkar_path_list[3].linear_path.segments
    score, karsim_alignment, omkar_alignment = align_paths(path1, path2)
    print(score)
    str1, str2 = tostring_alignment(index_to_segment, karsim_alignment, omkar_alignment)
    print(str1)
    print(str2)


def test_manual_correction():
    file = 'tmp/cluster_6_runexample.txt'
    index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file(file)
    path1 = karsim_path_list[1].linear_path.segments
    path2 = omkar_path_list[1].linear_path.segments
    score, karsim_alignment, omkar_alignment = align_paths(path1, path2)
    print(score)
    str1, str2 = tostring_alignment(index_to_segment, karsim_alignment, omkar_alignment)
    print(str1)
    print(str2)


if __name__ == "__main__":
    test_NW_aligner()
    test_manual_correction()
