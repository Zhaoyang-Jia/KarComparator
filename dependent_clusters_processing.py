from Structures import *


def genome_wide_mutual_breaking(karsim_path_list, omkar_path_list):
    """
    in a complete bipartite graph fashion perform mutual breaking
    :param karsim_path_list:
    :param omkar_path_list:
    :return:
    """
    for karsim_path in karsim_path_list:
        for omkar_path in omkar_path_list:
            karsim_path.generate_mutual_breakpoints(omkar_path, mutual=True)


def form_dependent_clusters(karsim_path_list, omkar_path_list, output_dir):
    # check if the two path_list have mutually broken segments, only runs if this is ensured
    check_all_segments = []
    for path in karsim_path_list:
        check_all_segments += path.linear_path.segments
    for path in omkar_path_list:
        check_all_segments += path.linear_path.segments
    check_path = Path(Arm(check_all_segments, 'verifier'))
    if not check_path.is_disjoint():
        raise RuntimeError('input is not disjoint')

    karsim_origins = []
    omkar_origins = []

    for path in karsim_path_list:
        karsim_origins.append(path.get_origins())
    for path in omkar_path_list:
        omkar_origins.append(path.get_origins())

    dependent_clusters = form_least_disjoint_supergroups(omkar_origins + karsim_origins)
    clustered_origin_chr = {}
    clustered_karsim_path = {}
    clustered_omkar_path = {}
    for cluster_ind in range(len(dependent_clusters)):
        clustered_origin_chr[cluster_ind] = dependent_clusters[cluster_ind]

        new_karsim_cluster = []
        for path in karsim_path_list:
            path_origin_chr = path.path_chr.split('-')[0]
            if path_origin_chr in dependent_clusters[cluster_ind]:
                new_karsim_cluster.append(path)
        clustered_karsim_path[cluster_ind] = new_karsim_cluster

        new_omkar_cluster = []
        for path in omkar_path_list:
            path_origin_chr = path.path_chr.split('-')[0]
            if path_origin_chr in dependent_clusters[cluster_ind]:
                new_omkar_cluster.append(path)
        clustered_omkar_path[cluster_ind] = new_omkar_cluster

    for cluster_ind in range(len(clustered_origin_chr)):
        all_segments = set()
        for path in clustered_karsim_path[cluster_ind]:
            for segment in path.linear_path.segments:
                new_segment = segment.duplicate()
                if not new_segment.direction():
                    new_segment.invert()
                all_segments.add(new_segment)
        for path in clustered_omkar_path[cluster_ind]:
            for segment in path.linear_path.segments:
                new_segment = segment.duplicate()
                if not new_segment.direction():
                    new_segment.invert()
                all_segments.add(new_segment)
        all_segments = sorted(list(all_segments))

        output_str = ""
        # cluster info + segment indexing
        output_str += "cluster<{}>, origin_chr<{}>, n_karsim<{}>, n_omkar<{}>\n".format(cluster_ind,
                                                                                        clustered_origin_chr[cluster_ind],
                                                                                        len(clustered_karsim_path[cluster_ind]),
                                                                                        len(clustered_omkar_path[cluster_ind]))
        segment_dict = {}  # also form dict for later parsing
        for segment_ind in range(len(all_segments)):
            output_str += "{}\t{}\t{}\t{}\t{}\n".format(segment_ind,
                                                        all_segments[segment_ind].chr_name,
                                                        all_segments[segment_ind].start,
                                                        all_segments[segment_ind].end,
                                                        len(all_segments[segment_ind]))
            segment_dict[all_segments[segment_ind]] = str(segment_ind)
        output_str += "---\n"

        # karsim path
        for path in clustered_karsim_path[cluster_ind]:
            output_str += path.tostring_path_by_index(segment_dict) + '\n'
        output_str += "---\n"

        # omkar path
        for path in clustered_omkar_path[cluster_ind]:
            output_str += path.tostring_path_by_index(segment_dict) + '\n'
        output_str += "---\n"

        with open("{}/cluster_{}.txt".format(output_dir, str(cluster_ind)), 'w') as fp_write:
            fp_write.write(output_str)


def form_least_disjoint_supergroups(all_origins):
    """
    Given a list of values, form supergroups that are the smallest size but disjoint in the initial grouping
    :param all_origins:
    :return:
    """
    class Bin:
        def __init__(self, chr_components: set):
            self.chr_components = set(chr_components)

        def __len__(self):
            return len(self.chr_components)

        def includes_chr(self, input_chr):
            return input_chr in self.chr_components

        def chr_add(self, input_chr):
            self.chr_components.add(input_chr)

        def merge_bins(self, other):
            new_components = self.chr_components
            for new_chr in other.chr_components:
                new_components.add(new_chr)
            new_bin = Bin(new_components)
            return new_bin

        def get_sorted_chr_list(self):
            chr_list = list(self.chr_components)

            def custom_sort_key(chromosome):
                numerical_value = chromosome[3:]
                if numerical_value.upper() == 'X':
                    numerical_value = 23
                elif numerical_value.upper() == 'Y':
                    numerical_value = 24
                else:
                    numerical_value = int(numerical_value)
                return numerical_value

            chr_list = sorted(chr_list, key=custom_sort_key)
            return chr_list

    bin_list = []
    for origin in all_origins:
        bin_to_construct = Bin(set())
        for origin_chr in origin:
            # if an intersection found in a bin, add all members of the group into that bin
            for bin_itr in bin_list:
                if bin_itr.includes_chr(origin_chr):
                    bin_to_construct = bin_to_construct.merge_bins(bin_itr)

        # search rest of the chr in the merged bin, so we can merge their corresponding bins also
        bin_to_remove = []
        for bin_ind in range(len(bin_list)):
            current_bin_itr = bin_list[bin_ind]
            if current_bin_itr.chr_components.intersection(bin_to_construct.chr_components):
                bin_to_remove.append(bin_ind)
                bin_to_construct = bin_to_construct.merge_bins(current_bin_itr)

        if len(bin_to_construct) == 0:
            # none of the origin chrs is in a bin yet
            bin_to_construct = Bin(origin)
        else:
            new_bin_list = [value for index, value in enumerate(bin_list) if index not in bin_to_remove]
            bin_list = new_bin_list

        bin_list.append(bin_to_construct)

    output = []
    for bin_itr in bin_list:
        output.append(bin_itr.get_sorted_chr_list())
    return output


def test_form_least_disjoint_supergroups():
    test_case1 = [[1, 2, 3],
                  [4, 5],
                  [6],
                  [7],
                  [8],
                  [1, 2],
                  [3, 4],
                  [5, 6],
                  [7, 8]]
    print(form_least_disjoint_supergroups(test_case1))

    test_case2 = [[1, 2, 3],
                  [4, 5],
                  [6],
                  [7],
                  [1],
                  [2, 3, 4, 5, 6],
                  [7]]
    print(form_least_disjoint_supergroups(test_case2))

    test_case3 = [[1, 2, 3],
                  [4, 5],
                  [6],
                  [7],
                  [1],
                  [2, 3],
                  [4, 5],
                  [6, 7]]
    print(form_least_disjoint_supergroups(test_case3))


if __name__ == "__main__":
    pass
