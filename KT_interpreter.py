class Aligned_Haplotype:
    id: int  # for debugging
    alignment_score: int
    mt_aligned: [str]
    wt_aligned: [str]
    block_indices: {int: (int, int)}  # [block_id: (a, b)] where event block are [a, b)
    mt_blocks: {int: (str,)}  # insertion discordant blocks
    wt_blocks: {int: (str,)}  # deletion discordant blocks
    concordant_blocks: {int: (str,)}
    discordant_block_assignment: {int: str}  # block index: event name, exactly one event can be assigned to one block; only for discordant blocks
    size_dict: {str: int}  # seg_name: bp size

    def __init__(self, mt_aligned, wt_aligned, alignment_score, size_dict, id):
        self.id = id
        self.mt_aligned = mt_aligned
        self.wt_aligned = wt_aligned
        self.alignment_score = alignment_score
        self.size_dict = size_dict

        ## get event_blocks
        event_block = []
        block_types = []  # [str], a block is either ins/del
        # we know len(wt_hap) == len(mt_hap) as they are global alignments
        hap_len = len(self.mt_aligned)
        seg_idx = 0
        while seg_idx < hap_len:
            mt_seg = self.mt_aligned[seg_idx]
            wt_seg = self.wt_aligned[seg_idx]
            if mt_seg == wt_seg:
                seg_idx += 1
                continue
            else:
                # this will be the first seg of (potentially) a section, that is indel
                if wt_seg == '-':
                    # insertion
                    # block ends whenever continuous section OR insertion section ends
                    continuous_end = continuous_extension(self.mt_aligned, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if self.wt_aligned[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('ins')
                elif mt_seg == '-':
                    # deletion
                    continuous_end = continuous_extension(self.wt_aligned, seg_idx)
                    ins_end = seg_idx + 1
                    while ins_end < hap_len:
                        if self.mt_aligned[ins_end] == '-':
                            ins_end += 1
                        else:
                            break
                    block_end = min(continuous_end, ins_end)
                    block_types.append('del')
                else:
                    raise ValueError('mismatch not allowed')
                event_block.append((seg_idx, block_end))
                seg_idx = block_end

        ## extract concordant blocks, mt/wt blocks
        discordant_blocks = {event_block_boundaries: block_types[event_block_idx] for event_block_idx, event_block_boundaries in enumerate(event_block)}
        # all gaps in discordant blocks are concordant blocks
        self.block_indices = {}
        self.wt_blocks = {}
        self.mt_blocks = {}
        self.concordant_blocks = {}
        self.discordant_block_assignment = {}
        block_id = 0
        p_endpoint = 0
        for discordant_block in discordant_blocks:
            c_startpoint = int(discordant_block[0])
            if p_endpoint < c_startpoint:
                # fill gap with concordant blocks
                # between two discordant blocks, the region belongs to at most 1 concordant block because the wt is always continuous
                # (i.e. will not split into 2 blocks)
                self.block_indices[block_id] = (p_endpoint, c_startpoint)
                self.concordant_blocks[block_id] = wt_aligned[p_endpoint: c_startpoint]
                block_id += 1
            self.block_indices[block_id] = discordant_block
            self.discordant_block_assignment[block_id] = ''
            if discordant_blocks[discordant_block] == 'del':
                self.wt_blocks[block_id] = wt_aligned[discordant_block[0]: discordant_block[1]]
            elif discordant_blocks[discordant_block] == 'ins':
                self.mt_blocks[block_id] = mt_aligned[discordant_block[0]: discordant_block[1]]
            block_id += 1
            p_endpoint = discordant_block[1]
        # add the last concordant block, if present
        if p_endpoint < len(self.wt_aligned):
            self.block_indices[block_id] = (p_endpoint, len(self.wt_aligned))
            self.concordant_blocks[block_id] = wt_aligned[p_endpoint: len(self.wt_aligned)]

    def __str__(self):
        return "ID<{}>".format(self.id)

    def get_block_segs_and_block_type(self, block_idx):
        """
        given block_idx, get block segs and block type
        :param block_idx:
        :return:
        """
        if block_idx in self.wt_blocks:
            return self.wt_blocks[block_idx], self.wt_blocks
        elif block_idx in self.mt_blocks:
            return self.mt_blocks[block_idx], self.mt_blocks
        elif block_idx in self.concordant_blocks:
            return self.concordant_blocks[block_idx], self.concordant_blocks
        else:
            raise RuntimeError('block idx not found in the Haplotype')

    def report_SV(self, event_blocks, event_types):
        """
        :param event_blocks: cumulative dict, where value is a list of event_block, this helps to group paired-events that were on different chromosomes
        :param event_types: cumulative dict, where value is a list of event_type, same reason as above
        :return:
        """
        def get_block_boundaries(block_idx):
            if block_idx == 0:
                # first block
                left = 'p-ter'
            else:
                previous_block, _ = self.get_block_segs_and_block_type(block_idx - 1)
                left = previous_block[-1]
            if block_idx == len(self.block_indices) - 1:
                # last block
                right = 'q-ter'
            else:
                next_block, _ = self.get_block_segs_and_block_type(block_idx + 1)
                right = next_block[0]
            return left, right

        # TODO: conglomerate event boundaries, (i.e. report one pair of boundaries and event segment for inv, as it is always consecutive)
        for mt_block_idx, mt_block in self.mt_blocks.items():
            mt_block_str = ','.join(list(mt_block))
            event_id = int(self.discordant_block_assignment[mt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[mt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(mt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.mt({}).{}.{}'.format(self.id, mt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.mt({}).{}.{}'.format(self.id, mt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]
        for wt_block_idx, wt_block in self.wt_blocks.items():
            wt_block_str = ','.join(list(wt_block))
            event_id = int(self.discordant_block_assignment[wt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[wt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(wt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.wt({}).{}.{}'.format(self.id, wt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.wt({}).{}.{}'.format(self.id, wt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]
        return event_blocks, event_types

    def search_seed(self, query_section, seed_type, d, eps):
        """
        :param query_section:
        :param seed_type: 'ins' OR 'del', which type of seed to be searched of
        :param d:
        :param eps:
        :return: block_idx if found, -1 if not
        """
        if seed_type == 'del':
            query_block = self.wt_blocks
        elif seed_type == 'ins':
            query_block = self.mt_blocks
        else:
            raise ValueError('invalid seed type')
        for c_block_idx, c_block in query_block.items():
            if len(self.discordant_block_assignment[c_block_idx]) != 0:
                # already has assignment
                # TODO: can improve this by bipartite matching, instead of greedy matching
                continue
            seed_start, seed_end = is_seeded(c_block, query_section, self.size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                return c_block_idx
        return -1


def interpret_haplotypes(mt_hap_list: [[str]], wt_hap_list: [[str]], segment_size_dict: {str: int}, d=200000, eps=200000):
    """
    Input a haplotype and a WT, report SVs; hap from mt list must correspond to hap from wt list
    :param segment_size_dict: {str: int} segment mapped to size of the segment (eg. 11: 12,345)
    :param mt_hap_list:
    :param wt_hap_list:
    :return:
    """
    event_id = 1  # static variable, keep different events separate, keep paired events under same id (eg. balanced trans)
    hap_id = 1  # used for naming haplotypes, preserves the order they come-in in mt_hap/wt_hap lists
    aligned_haplotypes = []
    for idx, mt_hap in enumerate(mt_hap_list):
        wt_hap = wt_hap_list[idx]
        a1, a2, a3 = lcs(mt_hap, wt_hap, segment_size_dict)
        aligned_haplotypes.append(Aligned_Haplotype(a2, a3, a1, segment_size_dict, hap_id))
        hap_id += 1

    def genomewide_seed_search(query_section, query_block_idx, query_hap_idx, query_type):
        """
        :param query_section: (segs) to be searched for
        :param query_block_idx: needed for writing discordant_block_assignment if seed found
        :param query_hap_idx: needed to prioritize searching inter-chr first, then intra-chr if inter-chr not found
        :param query_type:
        :return:
        """
        seed_found = False
        for aligned_hap_idx1, aligned_hap1 in enumerate(aligned_haplotypes):
            if aligned_hap_idx1 == query_hap_idx:
                # first get inter-chr event
                continue
            seeded_block_idx = aligned_hap1.search_seed(query_section, query_type, d=d, eps=eps)
            if seeded_block_idx != -1:
                aligned_hap1.discordant_block_assignment[seeded_block_idx] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[query_block_idx] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
                break
        if not seed_found:
            # at last, check intra-chr event
            seeded_block_idx = aligned_haplotypes[query_hap_idx].search_seed(query_section, query_type, d=d, eps=eps)
            if seeded_block_idx != -1:
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[seeded_block_idx] = 'balanced_translocation,{}'.format(event_id)
                aligned_haplotypes[query_hap_idx].discordant_block_assignment[query_block_idx] = 'balanced_translocation,{}'.format(event_id)
                seed_found = True
        if seed_found:
            return True
        else:
            return False

    ### Order of resolving: balanced-translocation, inv, dup-inv, tandem-dup, del, unbalanced-translocation
    ## resolve all balanced translocations
    # TODO: categorize into 1/2/3/4 breakpoint balanced translocation
    for aligned_hap_idx, aligned_hap in enumerate(aligned_haplotypes):
        for c_wt_block_idx, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_idx]) > 0:
                # already has assignment
                continue
            aligned_hap.discordant_block_assignment[c_wt_block_idx] = 'under investigation'  # block it from being matched in intra-chr event; maybe not necessary
            balanced_translocation_found = genomewide_seed_search(c_wt_block, c_wt_block_idx, aligned_hap_idx, 'ins')
            if balanced_translocation_found:
                event_id += 1
            else:
                # reset block assignment for query
                aligned_hap.discordant_block_assignment[c_wt_block_idx] = ''
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                # already has assignment
                continue
            aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'under investigation'  # block it from being matched in intra-chr event
            balanced_translocation_found = genomewide_seed_search(c_mt_block, c_mt_block_idx, aligned_hap_idx, 'del')
            if balanced_translocation_found:
                event_id += 1
            else:
                # reset block assignment for query
                aligned_hap.discordant_block_assignment[c_mt_block_idx] = ''

    ## inversion: mt{- k-} wt{k+ -} OR mt{k- -} wt{- k+}
    for aligned_hap in aligned_haplotypes:
        # case: mt{k- -} wt{- k+}, ins of inverted-seg + del of uninverted-seg
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if c_mt_block[0][-1] == "+":
                # not inverted
                continue
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            if c_mt_block_idx + 1 >= len(aligned_hap.block_indices) or \
                    (c_mt_block_idx + 1) not in aligned_hap.wt_blocks or \
                    len(aligned_hap.discordant_block_assignment[c_mt_block_idx + 1]) > 0:
                # this is the last block, next block is not a del-block, OR next del-block has assignment
                continue
            uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
            next_del_block = aligned_hap.wt_blocks[c_mt_block_idx + 1]
            seed_start, seed_end = is_seeded(next_del_block, uninverted_segs, segment_size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'inversion,{}'.format(event_id)
                aligned_hap.discordant_block_assignment[c_mt_block_idx + 1] = 'inversion,{}'.format(event_id)
                event_id += 1
        # case: mt{- k-} wt{k+ -}, del of uninverted-seg + ins of inverted-seg
        for c_wt_block_idx, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_idx]) > 0:
                continue
            if c_wt_block_idx + 1 >= len(aligned_hap.block_indices) or \
                    (c_wt_block_idx + 1) not in aligned_hap.mt_blocks or \
                    len(aligned_hap.discordant_block_assignment[c_wt_block_idx + 1]) > 0:
                # this is the last block, next block is not a ins-block, OR next ins-block has assignment
                continue
            inverted_segs = [seg[:-1] + '-' for seg in c_wt_block[::-1]]
            next_ins_block = aligned_hap.mt_blocks[c_wt_block_idx + 1]
            seed_start, seed_end = is_seeded(next_ins_block, inverted_segs, segment_size_dict, indel_direction='both', d=d, eps=eps)
            if seed_start != -1:
                aligned_hap.discordant_block_assignment[c_wt_block_idx] = 'inversion,{}'.format(event_id)
                aligned_hap.discordant_block_assignment[c_wt_block_idx + 1] = 'inversion,{}'.format(event_id)
                event_id += 1

    ## duplication inversion
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            # left-dup-inv: ins of inverted-seg + concordant block w/ uninverted-seg as seed
            if c_mt_block[0][-1] == "+":
                # not inverted
                continue
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            uninverted_segs = [seg[:-1] + '+' for seg in c_mt_block[::-1]]
            # this is not the last block, AND next block is concordant
            if c_mt_block_idx + 1 < len(aligned_hap.block_indices) and \
                    (c_mt_block_idx + 1) in aligned_hap.concordant_blocks:
                next_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx + 1]
                seed_start, seed_end = is_seeded(next_concordant_block, uninverted_segs, segment_size_dict, indel_direction='left', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'left_duplication_inversion,{}'.format(event_id)
                    event_id += 1
                    continue
        # right-dup-inv: concordant block w/ uninverted-seg as seed + ins of inverted-seg
            # this is not the first block, AND previous block is concordant
            if c_mt_block_idx - 1 >= 0 and \
                    (c_mt_block_idx - 1) in aligned_hap.concordant_blocks:
                previous_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx - 1]
                seed_start, seed_end = is_seeded(previous_concordant_block, uninverted_segs, segment_size_dict, indel_direction='right', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'right_duplication_inversion,{}'.format(event_id)
                    event_id += 1

    ## tandem-dup: mt{k+ k+}, wt{- k+} OR wt{k+ -}
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            block_segs = list(c_mt_block)
            # case: mt{k+ k+}, wt{- k+}, ins seg + concordant block w/ ins-seg as seed
            # this is not the last block, AND next block is concordant
            if c_mt_block_idx + 1 < len(aligned_hap.block_indices) and \
                    (c_mt_block_idx + 1) in aligned_hap.concordant_blocks:
                next_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx + 1]
                seed_start, seed_end = is_seeded(next_concordant_block, block_segs, segment_size_dict, indel_direction='left', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1
                    continue
            # case: mt{k+ k+}, wt{k+ -}, concordant block w/ ins-seg as seed + ins seg
            # this is not the first block, AND previous block is concordant
            if c_mt_block_idx - 1 >= 0 and \
                    (c_mt_block_idx - 1) in aligned_hap.concordant_blocks:
                previous_concordant_block = aligned_hap.concordant_blocks[c_mt_block_idx - 1]
                seed_start, seed_end = is_seeded(previous_concordant_block, block_segs, segment_size_dict, indel_direction='right', d=d, eps=eps)
                if seed_start != -1:
                    aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'tandem_duplication,{}'.format(event_id)
                    event_id += 1

    ## deletion: all remaining wt_blocks are deletions
    for aligned_hap in aligned_haplotypes:
        for c_wt_block_idx, c_wt_block in aligned_hap.wt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_wt_block_idx]) > 0:
                continue
            else:
                aligned_hap.discordant_block_assignment[c_wt_block_idx] = 'deletion,{}'.format(event_id)
                event_id += 1

    ## unbalanced translocation: all remaining mt_blocks are insertions (i.e. unbalanced translocations)
    for aligned_hap in aligned_haplotypes:
        for c_mt_block_idx, c_mt_block in aligned_hap.mt_blocks.items():
            if len(aligned_hap.discordant_block_assignment[c_mt_block_idx]) > 0:
                continue
            else:
                aligned_hap.discordant_block_assignment[c_mt_block_idx] = 'insertion,{}'.format(event_id)
                event_id += 1

    ## congregate events for report
    event_blocks = {}
    event_types = {}
    for aligned_hap in aligned_haplotypes:
        event_blocks, event_types = aligned_hap.report_SV(event_blocks, event_types)
    sorted_event_id = sorted(list(event_blocks.keys()))

    # if for the same event ID, we have different event types, we have an issue, otherwise, name the event as the singular name
    conglomerated_event_types = {}
    for c_id, event_type_list in event_types.items():
        c_event_type = event_type_list[0]
        for event_type in event_type_list:
            if event_type != c_event_type:
                raise ValueError('same ID, multiple event types')
        conglomerated_event_types[c_id] = c_event_type

    output_list = []
    for event_id in sorted_event_id:
        print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
        output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    return output_list


def continuous_extension(input_hap, idx_ptr):
    """
    start from the idx_ptr location, moving leftward until no longer can form a continuous section
    :param idx_ptr:
    :param input_hap:
    :return: final_ptr, where [idx_ptr, final_ptr) is a continuous section
    """
    final_ptr = idx_ptr + 1
    hap_len = len(input_hap)
    if idx_ptr >= hap_len:
        raise ValueError('idx_ptr out of bound')

    current_seg = input_hap[idx_ptr]
    section_sign = current_seg[-1]
    while True:
        if final_ptr == hap_len:
            break

        next_seg = input_hap[final_ptr]
        if next_seg == '-':
            break
        if section_sign == '+':
            if int(next_seg[:-1]) != int(current_seg[:-1]) + 1:
                break
        elif section_sign == '-':
            if int(next_seg[:-1]) != int(current_seg[:-1]) - 1:
                break
        else:
            raise ValueError('sign error')
        final_ptr += 1

    return final_ptr


def is_seeded(supergroup_section, cont_section, size_dict, indel_direction, d, eps):
    """
    whether cont_section is seeded in input_hap, with significant size (>d)
    :param supergroup_section:
    :param cont_section:
    :param size_dict:
    :param indel_direction: 'both', 'left', OR 'right': which direction to compute indels size on
    :param d: required seed overlap size
    :param eps: epsilon, allowed 1/2 indel size (i.e. 2 * eps >= indel size)
    :return: [start, end) indices in the supergroup
    """
    def sublist_idx(list1, list2):
        # list 1 is small, list 2 is large (presumed superlist)
        for i in range(len(list2) - len(list1) + 1):
            if list2[i:i + len(list1)] == list1:
                return i
        return None

    all_sublists = [cont_section[i:j + 1] for i in range(len(cont_section)) for j in range(i, len(cont_section))]
    sublist_sizes = [section_size(sublist, size_dict) for sublist in all_sublists]

    sublists_found = []
    seed_start_location_in_supergroup = []
    max_size = -1
    max_size_sublist_idx = -1
    for idx, sublist in enumerate(all_sublists):
        current_size = sublist_sizes[idx]
        if current_size < d:
            seed_start_location_in_supergroup.append(-1)
            continue
        current_start_location = sublist_idx(sublist, supergroup_section)
        if current_start_location is not None:
            sublists_found.append(sublist)
            seed_start_location_in_supergroup.append(current_start_location)
            if current_size > max_size:
                max_size = current_size
                max_size_sublist_idx = idx
        else:
            seed_start_location_in_supergroup.append(-1)

    if max_size_sublist_idx == -1:
        # not found
        return -1, -1

    max_size_sublist = all_sublists[max_size_sublist_idx]
    max_size_sublist_start_location = seed_start_location_in_supergroup[max_size_sublist_idx]
    max_size_sublist_end_location = seed_start_location_in_supergroup[max_size_sublist_idx] + len(max_size_sublist)

    del_size = section_size(cont_section, size_dict) - section_size(max_size_sublist, size_dict)
    left_ins_size = section_size(supergroup_section[:max_size_sublist_start_location], size_dict)
    right_ins_size = section_size(supergroup_section[max_size_sublist_end_location:], size_dict)

    if indel_direction == 'both':
        indel_size = del_size + left_ins_size + right_ins_size
        if indel_size <= 2 * eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    elif indel_direction == 'left':
        indel_size = del_size + left_ins_size
        if indel_size <= eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    elif indel_direction == 'right':
        indel_size = del_size + right_ins_size
        if indel_size <= eps:
            return max_size_sublist_start_location, max_size_sublist_end_location
    else:
        raise ValueError('indel_direction parameter input illegal')

    # indel size too large
    return -1, -1


def section_size(input_section, size_dict):
    tot_size = 0
    for seg in input_section:
        tot_size += size_dict[seg[:-1]]
    return tot_size


def lcs(list1, list2, size_dict):
    """
    longest common subsequence finding
    :param list1:
    :param list2:
    :param size_dict:
    :return:
    """
    indel_penalty_per_nt = -1
    alignment_1 = []
    alignment_2 = []
    scoring_matrix = [[0 for i in range(0, len(list2) + 1)] for j in range(0, len(list1) + 1)]
    backtrack_matrix = [["" for i in range(0, len(list2) + 1)] for j in range(0, len(list1) + 1)]

    # initialize starting grid
    scoring_matrix[0][0] = 0
    for row_index in range(1, len(list1) + 1):
        current_segment = list1[row_index - 1]
        len_current_segment = size_dict[current_segment[:-1]]  # remove sign
        scoring_matrix[row_index][0] = scoring_matrix[row_index - 1][0] + len_current_segment * indel_penalty_per_nt
        backtrack_matrix[row_index][0] = "down"

    for col_index in range(1, len(list2) + 1):
        current_segment = list2[col_index - 1]
        len_current_segment = size_dict[current_segment[:-1]]  # remove sign
        scoring_matrix[0][col_index] = scoring_matrix[0][col_index - 1] + len_current_segment * indel_penalty_per_nt
        backtrack_matrix[0][col_index] = "rigt"

    # DP recursion
    for row_index in range(1, len(list1) + 1):
        for col_index in range(1, len(list2) + 1):
            len_down_segment = size_dict[list1[row_index - 1][:-1]]
            down_value = scoring_matrix[row_index - 1][col_index] + len_down_segment * indel_penalty_per_nt

            len_right_segment = size_dict[list2[col_index - 1][:-1]]
            right_value = scoring_matrix[row_index][col_index - 1] + len_right_segment * indel_penalty_per_nt

            if list1[row_index - 1] != list2[col_index - 1]:
                # mismatch: not allowed
                diagonal_value = float('-inf')
            else:
                # match
                diagonal_value = scoring_matrix[row_index - 1][col_index - 1]

            if diagonal_value >= down_value and diagonal_value >= right_value:
                scoring_matrix[row_index][col_index] = diagonal_value
                backtrack_matrix[row_index][col_index] = "diag"
            elif down_value >= right_value:
                scoring_matrix[row_index][col_index] = down_value
                backtrack_matrix[row_index][col_index] = "down"
            else:
                scoring_matrix[row_index][col_index] = right_value
                backtrack_matrix[row_index][col_index] = "rigt"

    # backtracking
    final_score = scoring_matrix[len(list1)][len(list2)]
    current_row = len(list1)
    current_col = len(list2)

    while True:
        if current_row == 0 and current_col == 0:
            break
        if backtrack_matrix[current_row][current_col] == "diag":
            alignment_1.insert(0, list1[current_row - 1])
            alignment_2.insert(0, list2[current_col - 1])
            current_col -= 1
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "down":
            alignment_1.insert(0, list1[current_row - 1])
            alignment_2.insert(0, "-")
            current_row -= 1
        elif backtrack_matrix[current_row][current_col] == "rigt":
            alignment_1.insert(0, "-")
            alignment_2.insert(0, list2[current_col - 1])
            current_col -= 1
        else:
            raise RuntimeError("error in backtrack matrix")

    return final_score, alignment_1, alignment_2


if __name__ == "__main__":
    i_mt_hap1 = ['1+', '1+', '2+', '3-', '4+', '9+', '10+']
    i_mt_hap2 = ['7+', '8+', '5+', '6+']
    i_mt_hap3 = ['11+', '12+', '11-']
    i_mt_hap4 = ['14-', '13+','14+']
    i_mt_hap5 = ['15+', '16+', '17-', '16-', '18+']
    i_wt_hap1 = ['1+', '2+', '3+', '4+', '5+', '6+']
    i_wt_hap2 = ['7+', '8+', '9+', '10+']
    i_wt_hap3 = ['11+', '12+']
    i_wt_hap4 = ['13+', '14+']
    i_wt_hap5 = ['15+', '16+', '17+', '18+']
    i_mt_list = [i_mt_hap1, i_mt_hap2, i_mt_hap3, i_mt_hap4, i_mt_hap5]
    i_wt_list = [i_wt_hap1, i_wt_hap2, i_wt_hap3, i_wt_hap4, i_wt_hap5]
    i_size_dict = {str(i): 1 for i in range(19)}
    # print(lcs(i_mt_hap, i_wt_hap, i_size_dict))
    out = interpret_haplotypes(i_mt_list, i_wt_list, i_size_dict, 1, 1)
    print(out)
