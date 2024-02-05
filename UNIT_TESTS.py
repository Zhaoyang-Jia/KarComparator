from Structures import *
from read_KarSimulator_output import *
from read_OMKar_output import *
from dependent_clusters_processing import *
from read_cluster_file import *
from NW_aligner import *
from bipartite_matching import *
from COMPARISON_with_eulerian_swaps import *
from COMPARISON_with_graphs import *


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
    karsim_file = '/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/KarSimulator/23X_Angelman_r1.kt.txt'
    omkar_file = '/media/zhaoyang-new/workspace/KarSim/KarComparator/new_data_files/OMKar/23X_Angelman_r1.1.txt'
    forbidden_region_file = 'Metadata/acrocentric_telo_cen.bed'
    output_dir = 'tmp/'

    karsim_path_list = read_KarSimulator_output_to_path(karsim_file, forbidden_region_file)
    omkar_path_list = read_OMKar_output_to_path(omkar_file, forbidden_region_file)
    genome_wide_mutual_breaking(karsim_path_list, omkar_path_list)
    form_dependent_clusters(karsim_path_list, omkar_path_list, output_dir)


def test_read_cluster_file():
    file = 'debug_files/cluster_6_graph_test.txt'
    _, karsim_path_list, omkar_path_list = read_cluster_file(file)
    for path in karsim_path_list:
        print(path)
    for path in omkar_path_list:
        print(path)


def test_NW_aligner():
    file = '/media/zhaoyang-new/workspace/KarSim/KarComparator/batch_processing/cluster_files/12q14_microdeletion_r1cluster_4.txt'
    index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file(file)
    path1 = karsim_path_list[1].linear_path.segments
    path2 = omkar_path_list[1].linear_path.segments
    score, karsim_alignment, omkar_alignment = align_paths(path1, path2)
    print(score)
    str1, str2 = tostring_alignment(index_to_segment, karsim_alignment, omkar_alignment)
    print(str1)
    print(str2)


def test_manual_correction():
    file = 'debug_files/cluster_6_runexample.txt'
    index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file(file)
    path1 = karsim_path_list[1].linear_path.segments
    path2 = omkar_path_list[1].linear_path.segments
    score, karsim_alignment, omkar_alignment = align_paths(path1, path2)
    print(score)
    str1, str2 = tostring_alignment(index_to_segment, karsim_alignment, omkar_alignment)
    print(str1)
    print(str2)


def test_bipartite_matching():
    file = '/media/zhaoyang-new/workspace/KarSim/KarComparator/batch_processing/cluster_files/12q14_microdeletion_r1cluster_4.txt'
    index_to_segment, karsim_path_list, omkar_path_list = read_cluster_file(file)
    hungarian_algorithm_for_cluster(karsim_path_list, omkar_path_list, index_to_segment, verbose=True)


def test_form_graph():
    # TODO: complete writing this test
    file = 'new_data_files/cluster_files/23X_Early_onset_Alzheimer_r1cluster_14.txt'
    read_alignment_file(file, '')


def test_form_graph_2():
    draw_graph('new_data_files/cluster_files/23X_1q21_recurrent_microduplication_r1cluster_0.txt')


if __name__ == "__main__":
    # test_form_dependent_clusters()
    # test_NW_aligner()
    # test_manual_correction()

    # test_read_cluster_file()

    # test_NW_aligner()

    # test_bipartite_matching()

    # test_form_graph()

    test_form_graph_2()
