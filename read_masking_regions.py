from Structures import Arm, Segment


def read_masking_regions(masking_file) -> Arm:
    segment_list = []
    with open(masking_file) as fp_read:
        fp_read.readline()  # skip index line
        for line in fp_read:
            line = line.replace('\n', '').split('\t')
            new_segment = Segment(str(line[0]), int(line[1]), int(line[2]), str(line[3]))
            segment_list.append(new_segment)
    return Arm(segment_list, 'masking_regions')


def test():
    print(read_masking_regions('../Metadata/merged_masking_unique.bed'))


if __name__ == "__main__":
    test()
