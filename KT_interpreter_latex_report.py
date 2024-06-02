from KT_interpreter import *
from forbidden_region_processing import *
from KT_visualizer import *
import re


def format_report(interpreted_events, aligned_haplotypes, index_to_segment_dict):
    iscn_events = []
    gene_reports = []
    associated_event_already_reported = []
    for event in interpreted_events:
        main_str = ''
        iscn_interpretation = ''

        event_id = event[0]
        event_type = event[1]
        if event_id in associated_event_already_reported:
            continue
        cn_signature = 0
        cn_changed_genes = []
        cn_changed_genes_highlight = []
        if not event_type.startswith('balanced_translocation'):
            # then only 1 block for each event
            if len(event[2]) != 1:
                raise RuntimeError('not exactly 1 block notation')
            path_idx = int(event[2][0].split('.')[0])
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
                    bp2_chr = right_event_seg.chr_name
                    bp3_chr = left_event_seg.chr_name
                else:
                    bp2 = left_event_seg.end
                    bp3 = right_event_seg.start
                    bp2_chr = left_event_seg.chr_name
                    bp3_chr = right_event_seg.chr_name
            else:
                bp2 = left_event_seg.start
                bp3 = right_event_seg.end
                bp2_chr = left_event_seg.chr_name
                bp3_chr = right_event_seg.chr_name
            if event[2][0].split('.')[3] == 'p-ter':
                bp1 = None
                bp1_chr = None
            else:
                bp1_seg = index_to_segment_dict[int(event[2][0].split('.')[3][:-1])]
                bp1_chr = bp1_seg.chr_name
                if event[2][0].split('.')[3][-1] == '+':
                    bp1 = bp1_seg.start
                else:
                    bp1 = bp1_seg.end
            if event[2][0].split('.')[4] == 'q-ter':
                bp4 = None
                bp4_chr = None
            else:
                bp4_seg = index_to_segment_dict[int(event[2][0].split('.')[4][:-1])]
                bp4_chr = bp4_seg.chr_name
                if event[2][0].split('.')[4][-1] == '+':
                    bp4 = bp4_seg.start
                else:
                    bp4 = bp4_seg.end
            if bp1 is not None:
                bp1_band = get_band_location(bp1_chr, bp1)
            else:
                bp1_band = 'pter'
            bp2_band = get_band_location(bp2_chr, bp2)
            bp3_band = get_band_location(bp3_chr, bp3)
            if bp4 is not None:
                bp4_band = get_band_location(bp4_chr, bp4)
            else:
                bp4_band = 'qter'

            if event_type == 'deletion':
                if bp2_band != bp3_band:
                    main_str = 'del({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'del({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'deletion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = -1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'inversion':
                if bp2_band != bp3_band:
                    main_str = 'inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
            elif event_type == 'tandem_duplication':
                if bp2_band != bp3_band:
                    main_str = 'dup({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'dup({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'tandem duplication on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'left_duplication_inversion':
                if bp2_band != bp3_band:
                    main_str = 'left-dup-inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'left-dup-inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'left duplication inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'right_duplication_inversion':
                if bp2_band != bp3_band:
                    main_str = 'right-dup-inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'right-dup-inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                iscn_interpretation = 'right duplication inversion on Chr{}: {}'.format(path_chr, chr_range)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(path_chr, bp2, bp3)
            elif event_type == 'insertion':
                # different report format if insertion is from different chr
                if 'Chr' + path_chr == bp2_chr:
                    # TODO: check ISCN syntax if bp2_band == bp3_band
                    main_str = 'ins({})({}{}{})'.format(path_chr, bp1_band, bp2_band, bp3_band)
                else:
                    main_str = 'ins({};{})({};{}{})'.format(path_chr, bp2_chr, bp1_band, bp2_band, bp3_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                if bp1 is None:
                    bp1_text = bp1_band  # TODO: use 0/len(Chr) for pter/qter
                else:
                    bp1_text = format(bp1, ',d')
                iscn_interpretation = 'duplicated-insertion of Chr{}: {} into Chr{}: {} ({})'.\
                    format(bp2_chr[3:], chr_range, path_chr, bp1_text, bp1_band)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                cn_signature = 1
                cn_changed_genes, cn_changed_genes_highlight = report_cnv_genes_on_region(bp2_chr[3:], bp2, bp3)
            else:
                # continue
                raise RuntimeError('event not in allowed list')
        elif event_type.startswith("balanced_translocation_associated"):
            # TODO: only works with 2-break reciprocal balanced translocation
            o_event_id = int(event_type.split('<')[1].split('>')[0])
            # get extract the other event
            co_event_idx = -1
            for o_event_idx, event_itr in enumerate(interpreted_events):
                co_event_id = event_itr[0]
                if co_event_id == o_event_id:
                    co_event_idx = o_event_idx
                    break
            # check if o_event associate back
            o_event = interpreted_events[co_event_idx]
            if int(o_event[1].split('<')[1].split('>')[0]) != event_id:
                raise RuntimeError('more than 2-breaks detected')
            # get breakpoints and determine the swaps by number of qter/pter
            c_event_info = event[2]
            o_event_info = o_event[2]
            event_bps = c_event_info[0].split('.')[3:5] + c_event_info[1].split('.')[3:5] + o_event_info[0].split('.')[3:5] + o_event_info[1].split('.')[3:5]
            pter_idx = []
            qter_idx = []
            for event_bp_idx, event_bp_itr in enumerate(event_bps):
                if event_bp_itr == 'p-ter':
                    pter_idx.append(event_bp_idx)
                elif event_bp_itr == 'q-ter':
                    qter_idx.append(event_bp_idx)
            if len(pter_idx) + len(qter_idx) != 2 or len(pter_idx) == len(qter_idx):
                pass
                # raise RuntimeError('non-terminal 2-break reciprocal translocation detected')
            # locate endpoint of event segments, if p-side, choose right, if q-side, choose left
            indexed_event_segs1 = c_event_info[0].split('.')[2].split(')')[0].split('(')[1].split(',')
            indexed_event_segs2 = o_event_info[0].split('.')[2].split(')')[0].split('(')[1].split(',')
            typed_event_segs1 = []
            typed_event_segs2 = []
            for indexed_event_seg in indexed_event_segs1:
                typed_event_segs1.append(index_to_segment_dict[int(indexed_event_seg[:-1])])
            for indexed_event_seg in indexed_event_segs2:
                typed_event_segs2.append(index_to_segment_dict[int(indexed_event_seg[:-1])])
            seg1_left_bp = typed_event_segs1[0].start
            seg1_right_bp = typed_event_segs1[-1].end
            seg2_left_bp = typed_event_segs2[0].start
            seg2_right_bp = typed_event_segs2[-1].end
            seg1_left_band = get_band_location(typed_event_segs1[0].chr_name, seg1_left_bp)
            seg1_right_band = get_band_location(typed_event_segs1[-1].chr_name, seg1_right_bp)
            seg2_left_band = get_band_location(typed_event_segs2[0].chr_name, seg2_left_bp)
            seg2_right_band = get_band_location(typed_event_segs2[-1].chr_name, seg2_right_bp)
            if seg1_left_band[0] == 'p':
                is_pside = True
            elif seg1_left_band[0] == 'q':
                is_pside = False
            else:
                raise RuntimeError('illegal band name')
            chr1 = typed_event_segs1[0].chr_name[3:]  # assumes the segments have the same chr
            chr2 = typed_event_segs2[0].chr_name[3:]
            if is_pside:
                bp1 = typed_event_segs1[-1].end
                bp1_band = seg1_right_band
                bp2 = typed_event_segs2[-1].end
                bp2_band = seg2_right_band
            else:
                bp1 = typed_event_segs1[-1].start
                bp1_band = seg1_left_band
                bp2 = typed_event_segs2[-1].start
                bp2_band = seg2_left_band

            # if there is a sex chr, place it first
            flip_order = False
            if chr2.lower() == 'x':
                if chr1.lower() != 'x':
                    # temp_chr1, temp_bp1, temp_bp1_band = chr1, bp1, bp1_band
                    # chr1, bp1, bp1_band = chr2, bp2, bp2_band
                    # chr2, bp2, bp2_band = temp_chr1, temp_bp1, temp_bp1_band
                    flip_order = True
            elif chr2.lower() == 'y':
                if chr1.lower() != 'x' and chr1.lower() != 'y':
                    # temp_chr1, temp_bp1, temp_bp1_band = chr1, bp1, bp1_band
                    # chr1, bp1, bp1_band = chr2, bp2, bp2_band
                    # chr2, bp2, bp2_band = temp_chr1, temp_bp1, temp_bp1_band
                    flip_order = True
            elif chr1.lower() not in ['x', 'y'] and int(chr1) > int(chr2):
                flip_order = True
            if not flip_order:
                main_str = 't({};{})({};{})'.format(chr1, chr2, bp1_band, bp2_band)
                chr_range1 = chr_range_tostr(seg1_left_bp, seg1_right_bp, seg1_left_band, seg1_right_band)
                chr_range2 = chr_range_tostr(seg2_left_bp, seg2_right_bp, seg2_left_band, seg2_right_band)
                iscn_interpretation = 'balanced translocation between Chr{} and Chr{}, ' \
                                      'between segments Chr{}: {} and Chr{}: {}'. \
                    format(chr1, chr2,
                           chr1, chr_range1,
                           chr2, chr_range2)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(chr1, seg1_left_bp), (chr1, seg1_right_bp),
                                                                                     (chr2, seg2_left_bp), (chr2, seg2_right_bp)])
            else:
                main_str = 't({};{})({};{})'.format(chr2, chr1, bp2_band, bp1_band)
                chr_range1 = chr_range_tostr(seg2_left_bp, seg2_right_bp, seg2_left_band, seg2_right_band)
                chr_range2 = chr_range_tostr(seg1_left_bp, seg1_right_bp, seg1_left_band, seg1_right_band)
                iscn_interpretation = 'balanced translocation between Chr{} and Chr{}, ' \
                                      'between segments Chr{}: {} and Chr{}: {}'. \
                    format(chr2, chr1,
                           chr2, chr_range1,
                           chr1, chr_range2)
                bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(chr2, seg2_left_bp), (chr2, seg2_right_bp),
                                                                                     (chr1, seg1_left_bp), (chr1, seg1_right_bp)])
            associated_event_already_reported.append(o_event_id)
        elif event_type.startswith("balanced_translocation_unassociated"):
            event_info = event[2]
            if event_info[0].split('.')[2].startswith('mt'):
                del_idx = 0
                ins_idx = 1
            else:
                del_idx = 1
                ins_idx = 0
            ins_path_idx = int(event_info[ins_idx].split('.')[0])
            ins_chr = aligned_haplotypes[ins_path_idx].chrom[3:]
            indexed_event_segs = event_info[0].split('.')[2].split(')')[0].split('(')[1].split(',')
            typed_event_segs = []
            for indexed_event_seg in indexed_event_segs:
                typed_event_segs.append(index_to_segment_dict[int(indexed_event_seg[:-1])])
            event_seg_left_bp = typed_event_segs[0].start
            event_seg_right_bp = typed_event_segs[-1].end
            event_seg_left_band = get_band_location(typed_event_segs[0].chr_name, event_seg_left_bp)
            event_seg_right_band = get_band_location(typed_event_segs[-1].chr_name, event_seg_right_bp)
            if typed_event_segs[0].chr_name != typed_event_segs[-1].chr_name:
                raise RuntimeError('diff chr in event segs')
            else:
                event_seg_chr = typed_event_segs[0].chr_name[3:]
            ins_site_left_seg = event_info[ins_idx].split('.')[3]
            if ins_site_left_seg == 'p-ter':
                ins_site_left_bp = 0
                ins_site_left_band = 'p-ter'
            else:
                ins_site_left_bp = index_to_segment_dict[int(ins_site_left_seg[:-1])].start
                ins_site_left_band = get_band_location('chr' + ins_chr, ins_site_left_bp)
            main_str = 'ins-t({};{})({};{}{})'.format(ins_chr, event_seg_chr, ins_site_left_band, event_seg_left_band, event_seg_right_band)
            chr_range = chr_range_tostr(event_seg_left_bp, event_seg_right_bp, event_seg_left_band, event_seg_right_band)
            iscn_interpretation = 'balanced non-reciprocal translocation of Chr{}: {} into Chr{}: {}({})' \
                .format(event_seg_chr, chr_range, ins_chr, ins_site_left_bp, ins_site_left_band)
            bp_genes, bp_genes_highlight = report_on_genes_based_on_breakpoints([(event_seg_chr, event_seg_left_bp),
                                                                                 (event_seg_chr, event_seg_left_bp),
                                                                                 (ins_chr, ins_site_left_bp)])
        else:
            raise RuntimeError('illegal type assigned')
        if len(main_str) == 0 or len(iscn_interpretation) == 0:
            raise RuntimeError('missed interpretation')
        iscn_events.append([main_str, iscn_interpretation])
        gene_reports.append({'bp_genes': bp_genes,
                             'bp_genes_highlight': bp_genes_highlight,
                             'cnv': cn_signature,
                             'cnv_genes': cn_changed_genes,
                             'cnv_genes_highlight': cn_changed_genes_highlight})
    if len(iscn_events) != len(gene_reports):
        raise RuntimeError('unmatched reports')
    return iscn_events, gene_reports


def generate_latex_frontpage(title,
                             genefile_name,
                             breakpoint_reporting_proximity,
                             interpretation_insertion_threshold,
                             interpretation_deletion_threshold,
                             omkar_version='xxx'):
    output_str = ''

    def append_action(input_str, o):
        o = o + input_str + "\n"
        return o

    output_str += "\\catcode`\_=12\n"
    output_str += "\\documentclass[12pt]{article}\n"
    output_str += "\\input{macros}\n"
    output_str += "\\usepackage[letterpaper, margin=0.75in]{geometry}\n"
    output_str += "\\setcounter{secnumdepth}{0}\n"
    output_str += "\\usepackage{graphicx}\n"
    output_str += "\\usepackage{setspace}\n"
    output_str += "\\usepackage{titling}\n"
    output_str += "\\usepackage{enumitem}\n"
    # output_str += "\\usepackage[T1]{fontenc}\n"
    output_str += "\\renewcommand\\maketitlehookc{\\vspace{-10ex}}\n"
    # output_str += "\\graphicspath{ {./images/} }\n"

    output_str += "\n"
    output_str += "\\begin{document}\n"
    output_str += "\n"
    output_str += "\\title{" + title + "}\n"
    today = datetime.today()
    formatted_date = today.strftime("%b.%dth, %Y")
    output_str += "\\date{" + str(formatted_date) + "}\n"
    # output_str += "\\maketitle\n"
    output_str += "\n"
    # output_str += "\\textbf{{Samples Included: }} "
    # output_str += ', '.join(sample_names) + "\n"
    output_str += "\n"

    output_str = append_action('\\paragraph{Parameters}', output_str)
    output_str = append_action('\\begin{packed_itemize}', output_str)
    output_str = append_action('\\item OMKar version: {}'.format(omkar_version), output_str)
    output_str = append_action('\\item Genome assembly: hg38', output_str)
    output_str = append_action('\\item Genes: protein coding', output_str)
    output_str = append_action('\\item Breakpoint Gene Reporting Proximity: {}'.format(breakpoint_reporting_proximity), output_str)
    output_str = append_action('\\item Threshold for event insertion size: {}'.format(interpretation_insertion_threshold), output_str)
    output_str = append_action('\\item Threshold for event deletion size: {}'.format(interpretation_deletion_threshold), output_str)
    output_str = append_action('\\item Supported SV types:', output_str)
    output_str = append_action('\\begin{packed_itemize}', output_str)
    output_str = append_action('\\item Deletion', output_str)
    output_str = append_action('\\item Inversion', output_str)
    output_str = append_action('\\item Single/repeated Tandem-duplication', output_str)
    output_str = append_action('\\item Left/right Duplication-inversion', output_str)
    output_str = append_action('\\item 2/multi-break Reciprocal-balanced-translocation', output_str)
    output_str = append_action('\\item Nonreciprocal-balanced-translocation', output_str)
    output_str = append_action('\\item Duplicated-insertion', output_str)
    output_str = append_action('\\end{packed_itemize}', output_str)
    output_str = append_action('\\end{packed_itemize}', output_str)
    output_str += "\n"
    output_str += "\\newpage\n"

    return output_str


def batch_generate_latex_case_str(omkar_output_dir, image_dir, compile_image=False):
    """
    :param omkar_output_dir: assume files are named as (int).txt
    :param image_dir: relative path to the latex DIR, assume exists an image file with the same name, but (int).pdf
    :return:
    """
    image_suffix = '.pdf'

    final_str = ""
    cases_with_events = []
    files = [file for file in os.listdir(omkar_output_dir)]
    # files = sorted(files, key=int_file_keys)  # TODO: turn back on for files with int ID

    # highlight_files = ['23X_1q21_recurrent_microduplication_r1',
    #                    '23X_22q11_duplication_r2',
    #                    '23X_Angelman_r1',
    #                    '23Y_Cri_du_Chat_r1',
    #                    '23Y_WAGR_11p13_deletion_r2']
    highlight_files = []
    highlight_files = [i + '.1.txt' for i in highlight_files]
    files = [file for file in files if file not in highlight_files]
    files = highlight_files + files

    for file in files:
        # if True:
        # if file in ['3.txt', '39.txt', '49.txt', '12.txt', '45.txt']:
        if file == 'CMT1A_example':
            filename = file.split('.')[0]
            file_path = omkar_output_dir + file
            print(file)
            mt_indexed_lists, mt_path_chrs, segment_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
            mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
            wt_path_dict = generate_wt_from_OMKar_output(segment_dict)
            wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
            events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
            if len(events) == 0:
                continue
            dependent_clusters, cluster_events = form_dependent_clusters(events, aligned_haplotypes)
            print(dependent_clusters)
            final_str += "\\section{{Sample Id: {}}}\n".format(filename)
            final_str += "\n"
            ## iterate over all clusters
            n_clusters = len(dependent_clusters)
            for image_cluster_idx, (c_cluster, c_events) in enumerate(zip(dependent_clusters, cluster_events)):
                ## include all homologues
                event_chr = set()
                for cluster_idx in c_cluster:
                    event_chr.add(aligned_haplotypes[cluster_idx].chrom)
                hap_idx_to_plot = []
                for hap_idx, hap in enumerate(aligned_haplotypes):
                    if hap.chrom in event_chr:
                        hap_idx_to_plot.append(hap_idx)

                # if len(hap_idx_to_plot) > 4:
                #     raise RuntimeError('more than 4 chrom selected')
                c_aligned_haplotypes = [aligned_haplotypes[i] for i in hap_idx_to_plot]

                ## generate report text
                c_events = sort_events(c_events)
                iscn_events, genes_report = format_report(c_events, aligned_haplotypes, reverse_dict(segment_dict))
                ## generate image
                c_vis_input = generate_visualizer_input(c_events, c_aligned_haplotypes, segment_dict)

                def vis_key(input_vis):
                    chr_val = input_vis['chr'][3:]
                    if chr_val == "X":
                        return_val = 23.0
                    elif chr_val == "Y":
                        return_val = 24.0
                    else:
                        return_val = float(chr_val)
                    if input_vis['highlight']:
                        return_val += 0.5  # highlight always later
                    return return_val

                c_vis_input = sorted(c_vis_input, key=vis_key)

                image_prefix = "{}/{}_imagecluster{}".format(image_dir, filename, image_cluster_idx)
                image_path = image_prefix + '_rotated.png'
                overleaf_relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
                pycharm_relative_image_path = image_path
                if compile_image:
                    if len(c_vis_input) <= 4:
                        make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                    else:
                        make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)

                final_str += "\\subsection{{Event Cluster {} (of {})}}\n".format(image_cluster_idx + 1, n_clusters)
                final_str += "\n"
                # final_str += "\\noindent\n"
                # final_str += "\\begin{wrapfigure}{r}{0.5\\textwidth}\n"

                # final_str += "\\end{wrapfigure}\n"

                # Iterate through all events
                final_str += "\\begin{minipage}[t][][t]{0.5\\textwidth}\n"
                final_str += "\\vspace*{0pt}\n"
                final_str += "\\paragraph{SVs}\n"
                final_str += "\\medskip\n"
                final_str += "\\begin{flushleft}\n"
                SV_counter = 1
                for bullet_idx, (main_str, iscn_interpretation) in enumerate(iscn_events):
                    iscn_interpretation = hyperlink_iscn_interpretation(iscn_interpretation)
                    # format space-symbol better
                    iscn_interpretation = iscn_interpretation.replace(' on ', '&^&^&^&')
                    iscn_interpretation = iscn_interpretation.replace(' between ', '!@!@!@!@!')
                    iscn_interpretation = iscn_interpretation.replace(' ', '\\,')
                    iscn_interpretation = iscn_interpretation.replace('&^&^&^&', '\\,on ')
                    iscn_interpretation = iscn_interpretation.replace('!@!@!@!@!', ' between ')

                    # final_str += "\\item \\textbf{{{}}}.\\,{}\n".format(main_str, iscn_interpretation)
                    final_str += "\\textbf{{{}.\\;{}}}. {}\n\n".format(SV_counter, main_str, iscn_interpretation)
                    SV_counter += 1
                final_str += "\\end{flushleft}\n"
                final_str += "\n"
                final_str += "\\paragraph{Impacted genes in DDG2P}\\;\n"
                final_str += latex_gene_table(genes_report)
                final_str += "\\end{minipage}%\n"
                final_str += "\\hfill\n"
                final_str += "\\begin{minipage}[t][][t]{0.5\\textwidth}\n"
                final_str += "\\vspace*{0pt}\n"
                final_str += "\\centering\n"
                # final_str += "\\fbox{{\\includegraphics[width=\\linewidth]{{{}}}}}\n".format(overleaf_relative_image_path)
                final_str += "\\includegraphics[width=\\linewidth]{{{}}}\n".format(overleaf_relative_image_path)
                final_str += "\\captionof{{figure}}{{sample {}, event cluster {}}}\n".format(filename, image_cluster_idx + 1)
                final_str += "\\end{minipage}\n"
                final_str += "\\newpage\n\n"

    final_str += "\n"
    final_str += "\\end{document}\n"

    #         image_path = "{}/{}".format(image_dir, str(filename).zfill(3) + image_suffix)
    #         # image_path = "{}/{}".format(image_dir, str(filename) + image_suffix)
    #
    #         ## insert image
    #         if os.path.exists('latex_reports/' + image_path):
    #             # make sure the image file exists
    #             final_str += "\\section{{Sample Id: {}}}\n".format(filename)
    #             final_str += "\\begin{figure}[h!]\n"
    #             final_str += "\\centering\n"
    #             final_str += "\\includegraphics[width=3in]{{{}}}\n".format(image_path)
    #             final_str += "\\caption{{\\footnotesize Chromosomes with aberrant karyotypes}}\n"
    #             final_str += "\\label{{fig:karyotype_id{}}}\n".format(filename)
    #             final_str += "\\end{figure}\n"
    #             print('latex_reports/' + image_path, 'found')
    #         else:
    #             final_str += "\\section{{Sample Id: {}}}\n".format(filename)
    #             final_str += "\n"
    #             print('latex_reports/' + image_path, 'not found')
    #
    #         ## Iterate through all events
    #         final_str += "\\paragraph{Events}\n"
    #         final_str += "\\begin{packed_enum}\n"
    #         for bullet_idx, (main_str, iscn_interpretation) in enumerate(iscn_events):
    #             iscn_interpretation = hyperlink_iscn_interpretation(iscn_interpretation)
    #             final_str += "\\item {{\\bf {}}}. {}\n".format(main_str, iscn_interpretation)
    #         final_str += "\\end{packed_enum}\n"
    #
    #         final_str += "\n"
    #         final_str += "\\paragraph{Impacted genes in DDG2P}$\\;$\\\\\\\\\n"
    #         final_str += latex_gene_table(genes_report)
    #
    #         final_str += "\n"
    #         final_str += "\\newpage\n"
    # final_str += "\n"
    # final_str += "\\end{document}\n"
    return final_str, cases_with_events


def generate_latex_report(output_filename_prefix, front_page_str, batch_case_str):
    directory_path = os.path.dirname(output_filename_prefix) + '/'
    latex_path = output_filename_prefix + '.tex'
    with open(latex_path, 'w') as fp_write:
        fp_write.write(front_page_str)
        fp_write.write(batch_case_str)
    os.chdir('/'.join(output_filename_prefix.split('/')[:-1]))
    relative_latex_path = latex_path.split('/')[-1]
    subprocess.run(['pdflatex', relative_latex_path])


def chr_range_tostr(bpa, bpb, bpa_band, bpb_band):
    return "{}-{} ({} - {})".format(format(bpa, ',d'), format(bpb, ',d'), bpa_band, bpb_band)


def int_file_keys(f):
    return int(f.split('.')[0])


def latex_gene_table(genes_report):
    ## format genes to report
    genes_to_report = []
    for event_idx, report_dict in enumerate(genes_report):
        one_based_idx = event_idx + 1
        for entry in report_dict['bp_genes_highlight']:
            gene_entry = entry[0]
            disease_entry = entry[1]
            gene_name = gene_entry[0]
            gene_omim = gene_entry[1]
            disease_names = []
            disease_omim = []
            for disease in disease_entry:
                disease_names.append(disease[0])
                disease_omim.append(disease[1])
            genes_to_report.append({'SV': one_based_idx,
                                    'rationale': 'breakpoint proximal',
                                    'gene name': gene_name,
                                    'gene omim': gene_omim,
                                    'diseases': disease_names,
                                    'disease omims': disease_omim})
        for entry in report_dict['cnv_genes_highlight']:
            gene_entry = entry[0]
            disease_entry = entry[1]
            gene_name = gene_entry[0]
            gene_omim = gene_entry[1]
            disease_names = []
            disease_omim = []
            for disease in disease_entry:
                disease_names.append(disease[0])
                disease_omim.append(disease[1])
            if report_dict['cnv'] > 0:
                cnv_str = "+" + str(report_dict['cnv'])
            else:
                cnv_str = str(report_dict['cnv'])
            genes_to_report.append({'SV': one_based_idx,
                                    'rationale': 'CN' + cnv_str,
                                    'gene name': gene_name,
                                    'gene omim': gene_omim,
                                    'diseases': disease_names,
                                    'disease omims': disease_omim})

    ## form latex table
    if len(genes_to_report) == 0:
        return '\\\\\\quad None\n\n'
    return_str = "\\medskip\n"
    return_str += "{\\\\\\scriptsize\n"
    # return_str += "\\begin{flushleft}\n"
    return_str += "\\begin{tabular}{|llll|}\\hline\n"
    return_str += "SV & Rationale & Gene Name & Gene Omim  \\\\\\hline\n"
    for entry_idx, entry in enumerate(genes_to_report):
        new_line = '{SV} & {Rationale} & {Gene_Name} & {Gene_Omim} \\\\'. \
            format(SV=entry['SV'],
                   Rationale=entry['rationale'],
                   Gene_Name=entry['gene name'],
                   Gene_Omim=entry['gene omim'])
        if entry_idx == len(genes_to_report) - 1:
            return_str += new_line + "\\hline\n"
        else:
            return_str += new_line + "\n"
    return_str += "\\end{tabular}\n"
    return_str += "}\n"
    # return_str += "\\end{flushleft}\n"
    return return_str


def latex_hyperlink_coordinates(input_str, proximity=50000):
    return_dict = {}  # {replacement_string: hyperlinked_string}

    pattern = r'Chr(\d+|X|Y): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        start_pos = int(match.group(2).replace(',', ''))
        end_pos = int(match.group(3).replace(',', ''))
        ucsc_url = get_ucsc_url('chr' + match.group(1), start_pos, end_pos)
        hyperlinked_str = '\\href{{{}}}{{{}}}'.format(ucsc_url, replacement_str)
        return_dict[replacement_str] = hyperlinked_str

    pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        chrom = 'chr' + match.group(1)
        pos = int(match.group(2).replace(',', ''))
        c_chr_length = get_chr_length_from_forbidden_file(chrom)
        ucsc_url = get_ucsc_url(chrom, max(0, pos - proximity), min(c_chr_length, pos + proximity))
        hyperlinked_str = '\\href{{{}}}{{{}}}'.format(ucsc_url, replacement_str)
        return_dict[replacement_str] = hyperlinked_str

    return return_dict


def hyperlink_iscn_interpretation(input_str):
    hyperlinked_mapping = latex_hyperlink_coordinates(input_str)
    for replacement_str, hyperlinked_str in hyperlinked_mapping.items():
        input_str = input_str.replace(replacement_str, hyperlinked_str)
    return input_str


def test_latex(output_name, compile_image):
    output_path = 'latex_reports/{}'.format(output_name)
    batch_case_str, cases_in_report = batch_generate_latex_case_str(data_dir, image_dir, compile_image=compile_image)
    front_str = generate_latex_frontpage('{} Data'.format(' '.join(output_name.split('_'))),
                                         'hg38 all coding genes',
                                         50,
                                         200,
                                         200)
    generate_latex_report(output_path, front_str, batch_case_str)


if __name__ == "__main__":
    forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"
    # test_interpreter()
    # test_segs_union()
    # test_reciprocal_trans()
    # c_output_name, data_dir, image_dir = 'Dremsek', 'real_case_data/dremsek_OMKar_output_paths/', 'latex_reports/paul_dremsek_plots_new/'
    # c_output_name, data_dir, image_dir = 'Keyhole', 'real_case_data/keyhole_OMKar_output_paths/', 'latex_reports/keyhole_plots_new/'
    # c_output_name, data_dir, image_dir = 'Sunnyside', 'real_case_data/sunnyside_OMKar_output_paths/', 'latex_reports/sunnyside_plots_new/'
    c_output_name, data_dir, image_dir = 'ACC_Simulation', 'omkar_analyses_pipeline/builds/b14/omkar_paths/', 'latex_reports/ACC_simulation_plots/'
    # batch_case_str = batch_generate_latex_case_str(data_dir, image_dir)
    os.makedirs(image_dir, exist_ok=True)
    test_latex(c_output_name, compile_image=True)
