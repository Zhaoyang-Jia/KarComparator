from Structures import *
from forbidden_region_processing import read_forbidden_regions
from Start_Genome import *


def read_KarSimulator_output(KT_file, masking_file):
    genome = generate_genome_from_KT(KT_file)
    index_dict = genome.segment_indexing()
    t2_segments = get_t2_segments(masking_file)  # required as t2_segments start is unknown from KT file alone
    path_list = []

    for chromosome_itr in genome:
        if chromosome_itr.deleted:
            continue

        ## T1
        t1_segment = Segment(chromosome_itr.name[:-1], 0, chromosome_itr.t1_len - 1, "telomere1")
        segment_list = [t1_segment]

        ## P-arm
        for segment_itr in chromosome_itr.p_arm.segments:
            segment_itr.segment_type = "arm_region"
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        ## CEN
        cen_segment = chromosome_itr.centromere.segments[0].duplicate()
        cen_segment.segment_type = "centromere"
        cen_segment.kt_index = get_kt_index(index_dict, cen_segment)
        segment_list.append(cen_segment)

        ## Q-arm
        for segment_itr in chromosome_itr.q_arm.segments:
            segment_itr.segment_type = "arm_region"
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        ## T2
        current_t2_segment = None
        for find_t2_segment_itr in t2_segments.segments:
            if find_t2_segment_itr.chr_name == chromosome_itr.name[:-1]:
                current_t2_segment = find_t2_segment_itr
                break
        segment_list.append(current_t2_segment)
        path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))

    return path_list


def get_kt_index(input_index_dict, input_segment):
    for key in input_index_dict:
        if input_segment.same_segment_ignore_dir(key):
            return input_index_dict[key]

def get_t2_segments(masking_file):
    masking_arm = read_forbidden_regions(masking_file)
    t2_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == "telomere2":
            t2_segments.append(segment)
    return Arm(t2_segments, "t2_segments")


def segments_are_continuous(segment1: Segment, segment2: Segment):
    if segment1.direction():
        if segment1.end + 1 == segment2.start:
            return True
    else:
        if segment1.end - 1 == segment2.start:
            return True
    return False

def test():
    return_list = read_KarSimulator_output("sample_input/23X_Cri_du_Chat_r1.kt.txt",
                                           "Metadata/merged_forbidden_regions_unique.bed")
    for path in return_list:
        print(path)


if __name__ == "__main__":
    test()
