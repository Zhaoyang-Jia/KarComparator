from Structures import *
from read_masking_regions import read_masking_regions
from Start_Genome import *


def read_KT_to_path(KT_file, masking_file):
    genome = generate_genome_from_KT(KT_file, ordinal_info_included=True)
    history_dict = extract_history_by_chr_destination(genome)
    index_dict = genome.segment_indexing()
    t2_segments = get_t2_segments(masking_file)
    path_list = []

    history_counter = 1
    for chromosome_itr in genome:
        if chromosome_itr.deleted:
            continue
        t1_segment = Segment(chromosome_itr.name[:-1], 0, chromosome_itr.t1_len - 1, "telomere1")
        segment_list = [t1_segment]

        for segment_itr in chromosome_itr.p_arm.segments:
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        cen_segment = chromosome_itr.centromere.segments[0].duplicate()
        cen_segment.segment_type = "centromere"
        cen_segment.kt_index = get_kt_index(index_dict, cen_segment)
        segment_list.append(cen_segment)

        for segment_itr in chromosome_itr.q_arm.segments:
            segment_itr.kt_index = get_kt_index(index_dict, segment_itr)
            if segment_itr.direction():
                segment_itr.kt_index += "+"
            else:
                segment_itr.kt_index += "-"
            segment_list.append(segment_itr)

        t2_segment = None
        for find_t2_segment_itr in t2_segments.segments:
            if find_t2_segment_itr.chr_name == chromosome_itr.name[:-1]:
                t2_segment = find_t2_segment_itr
                break
        segment_list.append(t2_segment)

        # # if no event happened on this chr
        # if chromosome_itr.name not in history_dict:
        #     path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))
        #     continue

        path_list.append(Path(Arm(segment_list, "KT_path"), chromosome_itr.name, chromosome_itr.name[:-1]))

    # label SV segments from history
    # get inv/ins first
    for ordinal_history_entry_index in range(len(genome.ordinal_history)):
        current_path = None
        for path_itr in path_list:
            # access chr_to
            if path_itr.path_name == genome.history[ordinal_history_entry_index][3].name:
                current_path = path_itr
                break
        if current_path is None:
            continue  # chr deleted
        current_segment_list = current_path.linear_path.segments
        ordinal_history_entry = genome.ordinal_history[ordinal_history_entry_index]
        for ordinal_sub_entry in ordinal_history_entry:
            if ordinal_sub_entry[0] in ['inv', 'ins']:
                event_segment = ordinal_sub_entry[1]
                event_segment_in_path = current_segment_list[get_segment_index_from_ordinal(event_segment, event_segment.ordinal, current_segment_list)]
                event_segment_in_path.segment_type = ordinal_sub_entry[0] + ": SV" + str(ordinal_history_entry_index)

    # create indel ghost segments
    for ordinal_history_entry_index in range(len(genome.ordinal_history)):
        current_path = None
        for path_itr in path_list:
            # access chr_from
            if path_itr.path_name == genome.history[ordinal_history_entry_index][2].name:
                current_path = path_itr
                break
        if current_path is None:
            continue  # chr deleted
        current_segment_list = current_path.linear_path.segments
        ordinal_history_entry = genome.ordinal_history[ordinal_history_entry_index]
        for ordinal_sub_entry in ordinal_history_entry:
            if ordinal_sub_entry[0] == 'del':
                event_segments = ordinal_sub_entry[1]  # directionality is preserved as mentioned in history
                copied_event_segments = []
                for event_segment_itr in event_segments:
                    new_segment = event_segment_itr.duplicate()
                    new_segment.segment_type = 'del: SV' + str(ordinal_history_entry_index)
                    if new_segment in index_dict:
                        new_segment.kt_index = index_dict[new_segment] + '+'
                    elif new_segment.duplicate.invert(inplace=False) in index_dict:
                        new_segment.kt_index = index_dict[new_segment.duplicate.invert(inplace=False)] + '-'

                    copied_event_segments.append(new_segment)
                left_boundary_segment = ordinal_sub_entry[2]
                right_boundary_segment = ordinal_sub_entry[3]
                if left_boundary_segment == 'CEN':
                    left_boundary_index = get_centromere_index(current_segment_list)
                elif left_boundary_segment != 'T1':
                    left_boundary_index = get_segment_index_from_ordinal(left_boundary_segment, left_boundary_segment.ordinal, current_segment_list)
                else:
                    left_boundary_index = 0

                if right_boundary_segment == 'CEN':
                    right_boundary_index = get_centromere_index(current_segment_list)
                elif right_boundary_segment != 'T2':
                    right_boundary_index = get_segment_index_from_ordinal(right_boundary_segment, right_boundary_segment.ordinal, current_segment_list)
                else:
                    right_boundary_index = len(current_segment_list) - 1

                ghost_multiplicity = right_boundary_index - left_boundary_index
                if ghost_multiplicity == 0:
                    raise RuntimeError('same left and right boundary found')
                elif ghost_multiplicity == 1:
                    ghost_for_insertions = []
                    for event_segment_itr in copied_event_segments:
                        ghost_for_insertions.append(event_segment_itr.duplicate())
                    current_segment_list[left_boundary_index+1:left_boundary_index+1] = ghost_for_insertions
                else:
                    # intersperse del in-between the left and right boundary
                    ghost_for_insertions = []
                    for event_segment_itr in copied_event_segments:
                        new_segment_itr = event_segment_itr.duplicate()
                        new_segment_itr.segment_type += ": " + str(ghost_multiplicity)
                        ghost_for_insertions.append(new_segment_itr)
                    for insertion_index in range(right_boundary_index, left_boundary_index, -1):
                        current_segment_list[insertion_index:insertion_index] = ghost_for_insertions

    return path_list


def get_kt_index(input_index_dict, input_segment):
    for key in input_index_dict:
        if input_segment.same_segment_ignore_dir(key):
            return input_index_dict[key]


def extract_history_by_chr_destination(genome):
    history_dict = {}
    for history_entry in genome.history:
        history_type = history_entry[0]
        history_arm = history_entry[1]
        history_destination_chr = history_entry[3].name
        if history_destination_chr not in history_dict:
            history_dict[history_destination_chr] = []
        history_dict[history_destination_chr].append(tuple([history_type, history_arm]))
    return history_dict


def get_segment_in_SV(history_entry, segment_list):
    event_name = history_entry[0]
    if event_name == 'deletion':
        return
    # for
    # for segment in segment_list:


def get_t2_segments(masking_file):
    masking_arm = read_masking_regions(masking_file)
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


def get_segment_index_from_ordinal(input_current_segment: Segment, input_segment_ordinal: int, input_segment_list: [Segment]):
    occurrence = 0
    for finder_ptr in range(len(input_segment_list)):
        if input_segment_list[finder_ptr].same_segment_ignore_dir(input_current_segment):
            if input_segment_list[finder_ptr].segment_type is None:
                occurrence += 1
            elif not input_segment_list[finder_ptr].segment_type.startswith('del'):
                # del's ghost nodes do not count into the multiplicity
                occurrence += 1
        if occurrence == input_segment_ordinal:
            return finder_ptr

    raise RuntimeError('segment not found')


def get_centromere_index(input_segment_list: [Segment]):
    for segment_index in range(len(input_segment_list)):
        if input_segment_list[segment_index].segment_type == 'centromere':
            return segment_index


def test():
    return_list = read_KT_to_path("/Users/zhaoyangjia/Library/CloudStorage/OneDrive-UCSanDiego/Bafna_Lab/KarSimulator/test_folder/23Xe10_r1.kt.txt",
                                  "../Metadata/merged_masking_unique.bed")
    for path in return_list:
        print(path.concise_str())


if __name__ == "__main__":
    test()
