from Structures import *


def generate_genome_from_KT(input_file: str) -> Genome:
    """
    directly generating Genome from a formatted KT file
    :param input_file: path to input file
    :return: Genome Object
    """
    full_KT = {}
    motherboard_segments = []
    centromere_segments = []
    chromosome_of_interest = set()
    segment_dict = {}

    # read-in segment indexing
    with open(input_file) as fp_read:
        line = fp_read.readline()  # column headers
        while True:
            line = fp_read.readline()
            if line[0] == '-':
                break
            info = line.replace('\n', '').split('\t')
            new_segment = Segment(info[1], int(info[2]), int(info[3]))
            segment_dict[info[0]] = new_segment
            chromosome_of_interest.add(info[1])
            if info[0][0] == "C":
                centromere_segments.append(new_segment)
            else:
                motherboard_segments.append(new_segment)

        line = fp_read.readline()  # column headers

        # construct slots
        for chromosome in chromosome_of_interest:
            full_KT[chromosome] = []

        # read-in the segments in each chromosome
        while True:
            line = fp_read.readline().replace('\n', '')
            if line[0] == '-':
                break
            info = line.split('\t')

            if info[1] == 'deleted':
                new_chromosome = Chromosome(info[0], Arm([], 'deleted'), Arm([], 'deleted'), int(info[2]),
                                            int(info[3]), Arm([], 'deleted'), True)
                full_KT[info[0][:-1]].append(new_chromosome)
                continue

            p_arm_segments = []
            q_arm_segments = []
            current_centromere_segments = []
            segment_indices = info[1].split(',')
            # p_arm
            current_index = 0
            while segment_indices[current_index][0] != 'C':
                index_value = int(segment_indices[current_index][:-1])
                index_direction = segment_indices[current_index][-1]
                current_segment = segment_dict[str(index_value)].duplicate()
                if index_direction == "-":
                    current_segment.invert()
                p_arm_segments.append(current_segment)
                current_index += 1
            # centromere
            index_direction = segment_indices[current_index][-1]
            current_segment = segment_dict[segment_indices[current_index][:-1]].duplicate()
            if index_direction == "-":
                current_segment.invert()
            current_centromere_segments.append(current_segment)
            current_index += 1
            # q_arm
            while current_index < len(segment_indices):
                index_value = int(segment_indices[current_index][:-1])
                index_direction = segment_indices[current_index][-1]
                current_segment = segment_dict[str(index_value)].duplicate()
                if index_direction == "-":
                    current_segment.invert()
                q_arm_segments.append(current_segment)
                current_index += 1

            new_chromosome = Chromosome(info[0], Arm(p_arm_segments, 'p'), Arm(q_arm_segments, 'q'),
                                        int(info[2]), int(info[3]),
                                        Arm(current_centromere_segments, 'centromere'))
            # append chromosome to slot in genome, ignore the last char in chromosome's unique ID
            full_KT[info[0][:-1]].append(new_chromosome)

    return Genome(full_KT, motherboard_segments, centromere_segments)


def get_segment_location(input_current_segment: Segment, input_segment_ordinal: int, input_chr: Chromosome):
    # print(input_current_segment)
    # print(input_segment_ordinal)
    output_segment_location = -1
    finder_ptr = 0
    if len(input_chr.p_arm) == 0:
        p_arm_exhausted = True
    else:
        p_arm_exhausted = False
    for occurrence in range(0, input_segment_ordinal):
        segment_not_matched = True
        while segment_not_matched:
            if not p_arm_exhausted:
                if input_chr.p_arm.segments[finder_ptr].same_segment_ignore_dir(input_current_segment):
                    if occurrence + 1 == input_segment_ordinal:
                        output_segment_location = finder_ptr
                    else:
                        finder_ptr += 1
                        if finder_ptr >= len(input_chr.p_arm.segments):
                            p_arm_exhausted = True
                            finder_ptr = 0
                    segment_not_matched = False
                else:
                    finder_ptr += 1
                    if finder_ptr >= len(input_chr.p_arm.segments):
                        p_arm_exhausted = True
                        finder_ptr = 0
            else:
                if input_chr.q_arm.segments[finder_ptr].same_segment_ignore_dir(input_current_segment):
                    if occurrence + 1 == input_segment_ordinal:
                        output_segment_location = finder_ptr
                    else:
                        finder_ptr += 1
                        if finder_ptr >= len(input_chr.q_arm.segments):
                            raise RuntimeError('segment not found with the correct ordinal in given chr')
                    segment_not_matched = False
                else:
                    finder_ptr += 1
                    if finder_ptr >= len(input_chr.q_arm.segments):
                        print(input_current_segment)
                        print(input_segment_ordinal)
                        print(input_chr)
                        raise RuntimeError('segment not found with the correct ordinal in given chr')
    return p_arm_exhausted, output_segment_location


def test():
    genome = generate_genome_from_KT('sample_input/23X_Cri_du_Chat_r1.kt.txt')
    genome.output_KT('sample_output/23X_Cri_du_Chat_r1.kt.txt')


if __name__ == "__main__":
    test()
