from Structures import *
from forbidden_region_processing import *


def read_OMKar_output_to_path(OMKar_output_file, forbidden_region_file):
    path_list, index_dict = read_OMKar_output(OMKar_output_file, return_segment_dict=True)
    label_path_with_forbidden_regions(path_list, forbidden_region_file)
    rotate_and_bin_path(path_list, forbidden_region_file)
    report_centromere_anomaly(path_list)
    return index_dict, path_list


def read_OMKar_output(file, return_segment_dict=False):
    segment_dict = {}
    path_list = []
    with open(file) as fp_read:
        fp_read.readline()

        for line in fp_read:
            line = line.replace("\n", "").split('\t')
            # documenting segments
            if line[0] == "Segment":
                chr_name = str(line[2])
                if chr_name == "23":
                    chr_name = "X"
                elif chr_name == "24":
                    chr_name = "Y"
                chr_name = "Chr" + chr_name
                start = int(line[3].split(".")[0])
                end = int(line[4].split(".")[0])
                segment_dict[int(line[1])] = Segment(chr_name, start, end, "OMKar_unlabeled")
            elif line[0].startswith("Path"):
                # print(line)
                line = line[0].split(" = ")
                path_name = line[0]
                line = line[1]
                line = line.split(" ")
                path_segments = []
                for segment_index_itr in line:
                    if len(segment_index_itr) == 0:
                        break
                    direction = segment_index_itr[-1]
                    segment_index_itr = int(segment_index_itr[:-1])
                    new_segment = segment_dict[segment_index_itr].duplicate()
                    new_segment.kt_index = str(segment_index_itr)
                    new_segment.kt_index += '+'
                    if direction == "+":
                        path_segments.append(new_segment)
                    elif direction == "-":
                        new_segment.invert()
                        path_segments.append(new_segment)
                    else:
                        print(direction)
                        raise ValueError("direction must be + or -")
                path_list.append(Path(Arm(path_segments, "solved_path"), path_name))

    if return_segment_dict:
        return path_list, segment_dict
    else:
        return path_list


def rotate_and_bin_path(path_list, forbidden_region_file):
    """
    only works if each path contains exactly one centromere, OW will bin according to t1+t2+centromere percentage,
    if still no, will bin according to overall chr-content percentage
    will mark path accordingly if centromere anomaly exists
    :param forbidden_region_file:
    :param path_list:
    :return: path_list
    """
    # isolate centromere
    forbidden_region_path = Path(read_forbidden_regions(forbidden_region_file), 'forbidden_regions', 'forbidden_regions')
    for path in path_list:
        path.generate_mutual_breakpoints(other_path=forbidden_region_path, mutual=False)

    # get centromere, rotate if backward, and bin path
    for path in path_list:
        path_centromere = []
        path_telomeres = []
        for segment_itr in path.linear_path.segments:
            if 'centromere' in segment_itr.segment_type:
                path_centromere.append(segment_itr.duplicate())
            elif 'telomere' in segment_itr.segment_type:
                path_telomeres.append(segment_itr.duplicate())

        path_centromere_arm = Arm(path_centromere, 'centromeres')
        path_centromere_arm.merge_breakpoints()

        if len(path_centromere_arm.segments) == 1:
            path.path_chr = path_centromere_arm.segments[0].chr_name
            if not path_centromere_arm.segments[0].direction():
                # flip if centromere is backward
                # print(path.path_name + " is flipped")
                rotate_path(path)
        elif len(path_centromere_arm.segments) == 0:
            # print(path)
            chr_bin = get_highest_represented_chr(path_centromere + path_telomeres)
            if chr_bin is not None:
                # bin by cen + telo1 + telo2 %content
                path.path_chr = chr_bin
            else:
                # bin by the whole chr %centent
                path.path_chr = get_highest_represented_chr(path.linear_path.segments)
            path.path_chr += "-no centromere"
        else:
            # print(path)
            chr_bin = get_highest_represented_chr(path_centromere + path_telomeres)
            if chr_bin is not None:
                # bin by cen + telo1 + telo2 %content
                path.path_chr = chr_bin
            else:
                # bin by the whole chr %centent
                path.path_chr = get_highest_represented_chr(path.linear_path.segments)
            path.path_chr += "-multiple centromere: "
            for centromere_itr in path_centromere:
                path.path_chr += centromere_itr.chr_name + " "


def report_centromere_anomaly(path_list):
    for path in path_list:
        if "no centromere" in path.path_chr or "multiple centromere" in path.path_chr:
            print(path.get_path_notes())


def get_highest_represented_chr(segment_list):
    tally = {}
    for segment in segment_list:
        if segment.chr_name in tally:
            tally[segment.chr_name] += len(segment)
        else:
            tally[segment.chr_name] = len(segment)

    max_count = -1
    max_count_chr = None
    for key in tally:
        if tally[key] > max_count:
            max_count = tally[key]
            max_count_chr = key
    return max_count_chr


def rotate_path(input_path):
    segment_list = input_path.linear_path.segments
    segment_list.reverse()
    for segment_itr in segment_list:
        segment_itr.invert()
    input_path.linear_path.segments = segment_list

# def cmd_centromere_anomaly():
#     # TODO: verify errors before restoring
#     import argparse
#     parser = argparse.ArgumentParser(description="dicentromeric and acentromeric checker")
#     parser.add_argument("--file", type=str, dest='omkar_file', help="file path to OMKar's solved path")
#     args = parser.parse_args()
#
#     path_list = read_OMKar_output(args.omkar_file)
#     for path in path_list:
#         path.linear_path.merge_breakpoints()
#     path_list = rotate_and_bin_path(path_list, "Metadata/merged_forbidden_regions_unique.bed")
#     report_centromere_anomaly(path_list)


### pending updates
def get_segments_by_type(forbidden_region_file, segment_type):
    masking_arm = read_forbidden_regions(forbidden_region_file)
    centromere_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == segment_type:
            centromere_segments.append(segment)
    return Arm(centromere_segments, segment_type)


def bin_path_by_chr_content(input_path):
    tally = {}
    for segment in input_path.linear_path.segments:
        if segment.chr_name in tally:
            tally[segment.chr_name] += len(segment)
        else:
            tally[segment.chr_name] = len(segment)

    max_count = -1
    max_count_chr = None
    for key in tally:
        if tally[key] > max_count:
            max_count = tally[key]
            max_count_chr = key
    return max_count_chr


def test():
    path_list = read_OMKar_output("/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final/1q21-1_recurrent_microdeletion_r1.1/1q21-1_recurrent_microdeletion_r1.1.txt")
    path_list = rotate_and_bin_path(path_list, "../Metadata/merged_forbidden_regions_unique.bed")
    # report_centromere_anomaly(path_list)
    for path in path_list:
        print(path)


def test_read_OMKar_output():
    path_list = read_OMKar_output("sample_input/23Y_C ri_du_Chat_r1.1.txt")
    for path in path_list:
        print(path)


def test_read_OMKar_to_path():
    path_list = read_OMKar_output_to_path("sample_input/23Y_Cri_du_Chat_r1.1.txt", "Metadata/acrocentric_telo_cen.bed")
    for path in path_list:
        print(path)

if __name__ == "__main__":
    test_read_OMKar_to_path()
