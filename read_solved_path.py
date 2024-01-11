from Structures import *
from read_masking_regions import read_masking_regions


def read_solved_path(file):
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
                segment_dict[int(line[1])] = Segment(chr_name, start, end, "solved_path")
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
                    new_segment.kt_index += direction
                    if direction == "+":
                        path_segments.append(new_segment)
                    elif direction == "-":
                        new_segment.invert()
                        path_segments.append(new_segment)
                    else:
                        print(direction)
                        raise ValueError("direction must be + or -")
                path_list.append(Path(Arm(path_segments, "solved_path"), path_name))

    return path_list


def get_segments_by_type(masking_file, segment_type):
    masking_arm = read_masking_regions(masking_file)
    centromere_segments = []
    for segment in masking_arm.segments:
        if segment.segment_type == segment_type:
            centromere_segments.append(segment)
    return Arm(centromere_segments, segment_type)


def rotate_and_bin_path(path_list, masking_file):
    """
    only works if each path contains exactly one centromere
    :param masking_file:
    :param path_list:
    :return: path_list
    """
    centromere_arm = get_segments_by_type(masking_file, 'centromere')
    centromere_path = Path(centromere_arm)
    telomere1_arm = get_segments_by_type(masking_file, 'telomere1')
    telomere1_path = Path(telomere1_arm)
    telomere2_arm = get_segments_by_type(masking_file, 'telomere2')
    telomere2_path = Path(telomere2_arm)

    # isolate centromere
    for path in path_list:
        path.generate_mutual_breakpoints(other_path=centromere_path, mutual=False)
        path.generate_mutual_breakpoints(other_path=telomere1_path, mutual=False)
        path.generate_mutual_breakpoints(other_path=telomere2_path, mutual=False)

    # get centromere, rotate if backward, and bin path
    for path in path_list:
        path_centromere = []
        for centromere_segment in centromere_path.linear_path.segments:
            for segment_itr in path.linear_path.segments:
                if segment_itr.same_segment_ignore_dir(centromere_segment):
                    path_centromere.append(segment_itr)

        if len(path_centromere) == 1:
            flip_path_direction(path, centromere_path)
        elif len(path_centromere) == 0:
            path.path_chr = "no centromere"
            print(path.get_path_notes())
            if not flip_path_direction(path, telomere2_path):
                if not flip_path_direction(path, telomere1_path):
                    print("DEBUG: this path cannot be flipped: " + path.get_path_notes())
        else:
            path.path_chr = "multiple centromere: "
            for centromere_itr in path_centromere:
                path.path_chr += centromere_itr.chr_name + " "
            print(path.get_path_notes())
            if not flip_path_direction(path, telomere2_path):
                if not flip_path_direction(path, telomere1_path):
                    print("DEBUG: this path cannot be flipped: " + path.get_path_notes())

        path.path_chr = bin_path_by_chr_content(path)  # always bin by %chr content

    return path_list


def flip_path_direction(input_path, path_with_direction_anchor):
    for segment_itr1 in input_path.linear_path.segments:
        for segment_itr2 in path_with_direction_anchor.linear_path.segments:
            if segment_itr1.same_segment_ignore_dir(segment_itr2):
                if not segment_itr1.direction():
                    reversed_segment_list = []
                    for reversed_segment_itr in reversed(input_path.linear_path.segments):
                        new_segment = reversed_segment_itr.duplicate()
                        new_segment.invert()
                        reversed_segment_list.append(new_segment)
                    input_path.linear_path.segments = reversed_segment_list
                return True
    return False


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


def report_centromere_anomaly(path_list):
    print(len(path_list))
    for path in path_list:
        if path.path_chr.startswith("no centromere") or path.path_chr.startswith("multiple centromere"):
            print(path.get_path_notes())


def test():
    path_list = read_solved_path("/media/zhaoyang-new/workspace/KarSim/OMKar_outputs/simulation_final/1q21-1_recurrent_microdeletion_r1.1/1q21-1_recurrent_microdeletion_r1.1.txt")
    path_list = rotate_and_bin_path(path_list, "../Metadata/merged_masking_unique.bed")
    # report_centromere_anomaly(path_list)
    for path in path_list:
        print(path)


def cmd():
    import argparse
    parser = argparse.ArgumentParser(description="dicentromeric and acentromeric checker")
    parser.add_argument("--file", type=str, dest='omkar_file', help="file path to OMKar's solved path")
    args = parser.parse_args()

    path_list = read_solved_path(args.omkar_file)
    for path in path_list:
        path.linear_path.merge_breakpoints()
    path_list = rotate_and_bin_path(path_list, "Metadata/merged_masking_unique.bed")
    report_centromere_anomaly(path_list)


if __name__ == "__main__":
    cmd()
