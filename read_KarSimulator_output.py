from Structures import *
from forbidden_region_processing import *
from Karsimulator_Start_Genome import *


def read_KarSimulator_output_to_path(KarSimulator_output_file, forbidden_region_file):
    index_dict, path_list, event_histories = read_KarSimulator_output(KarSimulator_output_file, forbidden_region_file)
    label_path_with_forbidden_regions(path_list, forbidden_region_file)
    return index_dict, path_list


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

    return index_dict, path_list, genome.histories


def label_event_sv_edge(path_list, history_list):
    def find_path(path_name):
        for path_idx, path in enumerate(path_list):
            if path.path_name == path_name:
                return path_idx
        raise ValueError('Path not found')

    def find_seg_sublist(superlist_segments, sublist_segments):
        superlist_hash = [hash(seg) for seg in superlist_segments]
        sublist_hash = [hash(seg) for seg in sublist_segments]
        for i in range(len(sublist_hash) - len(sublist_hash) + 1):
            if superlist_hash[i: i+len(sublist_hash)] == sublist_hash:
                return i
        return None

    def invert_segment_list(input_segs):
        inverted_segs = []
        for seg in reversed(input_segs):
            new_seg = seg.duplicate()
            new_seg.invert()
            inverted_segs.append(new_seg)
        return inverted_segs

    sv_edges_list = []
    for hist_entry in history_list:
        event_type = hist_entry[0]
        c_sv_edges = []
        if event_type == 'inversion':
            c_path = path_list[find_path(hist_entry[1])]  # Chr_from = Chr_to
            event_segments = hist_entry[3]
            inverted_event_segments = invert_segment_list(event_segments)
            sublist_idx = find_seg_sublist(c_path.linear_path.segments, inverted_event_segments)
            if sublist_idx is not None:
                seg1 = c_path.linear_path.segments[sublist_idx - 1]
                seg2 = c_path.linear_path.segments[sublist_idx]
                seg3 = c_path.linear_path.segments[sublist_idx + len(inverted_event_segments) - 1]
                seg4 = c_path.linear_path.segments[sublist_idx + len(inverted_event_segments)]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'tandem duplication':
            c_path = path_list[find_path(hist_entry[1])]  # Chr_from = Chr_to
            event_segments = hist_entry[3] + hist_entry[3]
            sublist_idx = find_seg_sublist(c_path.linear_path.segments, event_segments)
            if sublist_idx is not None:
                seg1 = c_path.linear_path.segments[sublist_idx + len(event_segments) / 2 - 1]
                seg2 = c_path.linear_path.segments[sublist_idx + len(event_segments) / 2]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
        elif event_type == 'left duplication inversion':
            c_path = path_list[find_path(hist_entry[1])]  # Chr_from = Chr_to
            event_segments = invert_segment_list(hist_entry[3]) + hist_entry[3]
            sublist_idx = find_seg_sublist(c_path.linear_path.segments, event_segments)
            if sublist_idx is not None:
                seg1 = c_path.linear_path.segments[sublist_idx - 1]
                seg2 = c_path.linear_path.segments[sublist_idx]
                seg3 = c_path.linear_path.segments[sublist_idx + len(event_segments) / 2 - 1]
                seg4 = c_path.linear_path.segments[sublist_idx + len(event_segments) / 2]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'right duplication inversion':
            c_path = path_list[find_path(hist_entry[1])]  # Chr_from = Chr_to
            event_segments = hist_entry[3] + invert_segment_list(hist_entry[3])
            sublist_idx = find_seg_sublist(c_path.linear_path.segments, event_segments)
            if sublist_idx is not None:
                seg1 = c_path.linear_path.segments[sublist_idx + len(event_segments) / 2 - 1]
                seg2 = c_path.linear_path.segments[sublist_idx + len(event_segments) / 2]
                seg3 = c_path.linear_path.segments[sublist_idx + len(event_segments) - 1]
                seg4 = c_path.linear_path.segments[sublist_idx + len(event_segments)]
                c_sv_edges.append((seg1.chr_name, seg1.end, seg2.chr_name, seg2.start))
                c_sv_edges.append((seg3.chr_name, seg3.end, seg4.chr_name, seg4.start))
        elif event_type == 'balanced reciprocal translocation':
            # balanced reciprocal translocation always appear in pairs
            # TODO: redo sorting to make sure they happen in pairs
            pass
        sv_edges_list.append(c_sv_edges)

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
    index_dict, return_list = read_KarSimulator_output_to_path("sample_input/23X_Cri_du_Chat_r1.kt.txt",
                                           "Metadata/acrocentric_telo_cen.bed")
    for path in return_list:
        print(path)


if __name__ == "__main__":
    test()
