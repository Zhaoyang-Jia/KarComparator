from Report_Genes import *


class Aligned_Haplotype:
    id: int  # for debugging
    chrom: str
    alignment_score: int
    mt_aligned: [str]
    wt_aligned: [str]
    block_indices: {int: (int, int)}  # [block_id: (a, b)] where event block are [a, b)
    mt_blocks: {int: (str,)}  # insertion discordant blocks
    wt_blocks: {int: (str,)}  # deletion discordant blocks
    concordant_blocks: {int: (str,)}
    discordant_block_assignment: {int: str}  # block index: event name, exactly one event can be assigned to one block; only for discordant blocks
    size_dict: {str: int}  # seg_name: bp size

    def __init__(self, mt_aligned, wt_aligned, alignment_score, size_dict, id, chrom):
        self.id = id
        self.chrom = chrom
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

        for mt_block_idx, mt_block in self.mt_blocks.items():
            mt_block_str = ','.join(list(mt_block))
            event_id = int(self.discordant_block_assignment[mt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[mt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(mt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.{}.mt({}).{}.{}'.format(self.id, mt_block_idx, mt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.{}.mt({}).{}.{}'.format(self.id, mt_block_idx, mt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]
        for wt_block_idx, wt_block in self.wt_blocks.items():
            wt_block_str = ','.join(list(wt_block))
            event_id = int(self.discordant_block_assignment[wt_block_idx].split(',')[-1])
            event_type = self.discordant_block_assignment[wt_block_idx].split(',')[0]
            left_boundary_seg, right_boundary_seg = get_block_boundaries(wt_block_idx)
            if event_id in event_blocks:
                event_blocks[event_id].append('{}.{}.wt({}).{}.{}'.format(self.id, wt_block_idx, wt_block_str, left_boundary_seg, right_boundary_seg))
                event_types[event_id].append(event_type)
            else:
                event_blocks[event_id] = ['{}.{}.wt({}).{}.{}'.format(self.id, wt_block_idx, wt_block_str, left_boundary_seg, right_boundary_seg)]
                event_types[event_id] = [event_type]

        # sort each event_blocks by chr, then by block_idx
        for event_id, event_block_list in event_blocks.items():
            event_blocks[event_id] = sorted(event_block_list, key=lambda x: (int(x.split('.')[0]), int(x.split('.')[1])), reverse=False)

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


# def segs_union(seg_list1, seg_list2):
#     # Find the index of the overlapping region
#     overlap_index = len(seg_list1) - len(set(seg_list1) - set(seg_list2))
#     # Extend the union list to include elements from both lists
#     union_list = seg_list1[:overlap_index] + seg_list2
#     return union_list


def interpret_haplotypes(mt_hap_list: [[str]], wt_hap_list: [[str]], chrom_identities: [str], segment_size_dict: {str: int}, d=200000, eps=200000):
    """
    Input a haplotype and a WT, report SVs; hap from mt list must correspond to hap from wt list
    :param segment_size_dict: {str: int} segment mapped to size of the segment (eg. 11: 12,345)
    :param mt_hap_list:
    :param wt_hap_list:
    :param chrom_identities:
    :return:
    """
    event_id = 1  # static variable, keep different events separate, keep paired events under same id (eg. balanced trans)
    hap_id = 1  # used for naming haplotypes, preserves the order they come-in in mt_hap/wt_hap lists
    aligned_haplotypes = []
    for idx, mt_hap in enumerate(mt_hap_list):
        wt_hap = wt_hap_list[idx]
        c_chrom = chrom_identities[idx]
        a1, a2, a3 = lcs(mt_hap, wt_hap, segment_size_dict)
        aligned_haplotypes.append(Aligned_Haplotype(a2, a3, a1, segment_size_dict, hap_id, c_chrom))
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

    # consolidate neighboring blocks for events with paired blocks
    for event_id, event_type in event_types.items():
        event_type = event_type[0]
        event_block_list = event_blocks[event_id]
        if event_type == 'inversion':
            if len(event_block_list) != 2:
                raise RuntimeError('inversion needs to have 2 consecutive discordant blocks')
            block1_segs = event_block_list[0].split('.')[2].split('(')[1].split(')')[0].split(',')
            block2_segs = event_block_list[1].split('.')[2].split('(')[1].split(')')[0].split(',')
            block1_size = section_size(block1_segs, segment_size_dict)
            block2_size = section_size(block2_segs, segment_size_dict)
            bp1 = event_block_list[0].split('.')[3]
            bp2 = event_block_list[1].split('.')[4]
            if block1_segs[0][-1] == '-':
                # earlier is ins of inverted-segs
                blocks1_uninverted_segs = [seg[:-1] + '+' for seg in block1_segs[::-1]]
                block2_uninverted_segs = block2_segs
            elif block2_segs[0][-1] == '-':
                # later is ins of inverted-segs
                blocks1_uninverted_segs = block1_segs
                block2_uninverted_segs = [seg[:-1] + '+' for seg in block2_segs[::-1]]
            else:
                raise RuntimeError('inversion not found')
            if block1_size >= block2_size:
                c_block_segs = blocks1_uninverted_segs
            else:
                c_block_segs = block2_uninverted_segs
            c_block_str = '({})'.format(','.join(c_block_segs))
            new_block_list = ['{}.{}.{}.{}.{}'.format(event_block_list[0].split('.')[0],
                                                      event_block_list[0].split('.')[1] + ',' + event_block_list[1].split('.')[1],
                                                      c_block_str + ',' + event_block_list[0].split('.')[2] + ',' + event_block_list[1].split('.')[2],
                                                      bp1,
                                                      bp2)]
            event_blocks[event_id] = new_block_list
        elif event_type == 'balanced_translocation':
            # TODO: current solution only accounts for 2-break balanced translocation
            pass

    output_list = []
    for event_id in sorted_event_id:
        print('event<{}>,type<{}>,blocks<{}>'.format(event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
        output_list.append((event_id, conglomerated_event_types[event_id], event_blocks[event_id]))
    return output_list, aligned_haplotypes


def format_report(interpreted_events, aligned_haplotypes, index_to_segment_dict):
    main_bullets = []
    sub_bullets = []
    for event in interpreted_events:
        event_type = event[1]
        if event_type != 'balanced_translocation':
            # then only 1 block for each event
            if len(event[2]) != 1:
                raise RuntimeError('not exactly 1 block notation')
            path_idx = int(event[2][0].split('.')[0]) - 1
            path_chr = aligned_haplotypes[path_idx].chrom[3:]
            event_segs = event[2][0].split('.')[2].split(')')[0].split('(')[1].split(',')
            left_event_seg = index_to_segment_dict[int(event_segs[0][:-1])]
            right_event_seg = index_to_segment_dict[int(event_segs[-1][:-1])]
            left_event_seg_dir = event_segs[0][-1]
            right_event_seg_dir = event_segs[-1][-1]
            if left_event_seg_dir == '-':
                # we assume the event segments are continuous
                if right_event_seg_dir != '-':
                    raise RuntimeError('event segs not in the same direction')
                # only maintain the directionality if it is an INS()
                if event_type != 'insertion':
                    bp2 = right_event_seg.start
                    bp3 = left_event_seg.end
                else:
                    bp2 = left_event_seg.end
                    bp3 = right_event_seg.start
            else:
                bp2 = left_event_seg.start
                bp3 = right_event_seg.end
            if event[2][0].split('.')[3] == 'p-ter':
                bp1 = None
            else:
                if event[2][0].split('.')[3][-1] == '+':
                    bp1 = index_to_segment_dict[int(event[2][0].split('.')[3][:-1])].start
                else:
                    bp1 = index_to_segment_dict[int(event[2][0].split('.')[3][:-1])].end
            if event[2][0].split('.')[4] == 'q-ter':
                bp4 = None
            else:
                if event[2][0].split('.')[4][-1] == '+':
                    bp4 = index_to_segment_dict[int(event[2][0].split('.')[4][:-1])].start
                else:
                    bp4 = index_to_segment_dict[int(event[2][0].split('.')[4][:-1])].end
            if bp1 is not None:
                bp1_band = get_band_location('chr' + path_chr, bp1)
            else:
                bp1_band = 'pter'
            bp2_band = get_band_location('chr' + path_chr, bp2)
            bp3_band = get_band_location('chr' + path_chr, bp3)
            if bp4 is not None:
                bp4_band = get_band_location('chr' + path_chr, bp4)
            else:
                bp4_band = 'qter'

            def get_genes_near_breakpoints(breakpoints, proximity):
                pass

            if event_type == 'deletion':
                if bp2_band != bp3_band:
                    main_str = 'del({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'del({})({})'.format(path_chr, bp2_band)
                sub_str1 = 'deletion on Chr{}: {}({}) - {}({})'.format(path_chr, bp2, bp2_band, bp3, bp3_band)
                genes_near_bp = get_genes_in_region('chr' + path_chr, )
                sub_str2 = 'all genes near breakpoints: '

            main_bullets.append(main_str)
        else:
            # balanced trans
            main_bullets.append('balanced trans skipped')
            sub_bullets.append(('x', 'x', 'x'))



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


def test_interpreter():
    i_mt_hap1 = ['1+', '1+', '2+', '3-', '4+', '9+', '10+']
    i_mt_hap2 = ['7+', '8+', '5+', '6+']
    i_mt_hap3 = ['11+', '12+', '11-']
    i_mt_hap4 = ['14-', '13+', '14+']
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


# def test_segs_union():
#     l1 = ['7+', '8+']
#     l2 = ['7+']
#     l3 = segs_union(l1, l2)
#     print(l3)


if __name__ == "__main__":
    test_interpreter()
    # test_segs_union()
