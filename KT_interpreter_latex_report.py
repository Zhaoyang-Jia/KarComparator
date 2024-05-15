from KT_interpreter import *

def format_report(interpreted_events, aligned_haplotypes, index_to_segment_dict):
    main_bullets = []
    sub_bullets = []
    associated_event_already_reported = []
    for event in interpreted_events:
        event_id = event[0]
        event_type = event[1]
        if event_id in associated_event_already_reported:
            continue
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
                sub_str1 = 'deletion on Chr{}: {}'.format(path_chr, chr_range)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                sub_str4, sub_str5 = report_cnv_genes_on_region(path_chr, bp2, bp3, '-1 CN')
                sub_bullets.append((sub_str1, sub_str2, sub_str3, sub_str4, sub_str5))
            elif event_type == 'inversion':
                if bp2_band != bp3_band:
                    main_str = 'inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                sub_str1 = 'inversion on Chr{}: {}'.format(path_chr, chr_range)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                sub_bullets.append((sub_str1, sub_str2, sub_str3))
            elif event_type == 'tandem_duplication':
                if bp2_band != bp3_band:
                    main_str = 'dup({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'dup({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                sub_str1 = 'tandem duplication on Chr{}: {}'.format(path_chr, chr_range)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                sub_str4, sub_str5 = report_cnv_genes_on_region(path_chr, bp2, bp3, '+1 CN')
                sub_bullets.append((sub_str1, sub_str2, sub_str3, sub_str4, sub_str5))
            elif event_type == 'left_duplication_inversion':
                if bp2_band != bp3_band:
                    main_str = 'left-dup-inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'left-dup-inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                sub_str1 = 'left duplication inversion on Chr{}: {}'.format(path_chr, chr_range)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                sub_str4, sub_str5 = report_cnv_genes_on_region(path_chr, bp2, bp3, '+1 CN')
                sub_bullets.append((sub_str1, sub_str2, sub_str3, sub_str4, sub_str5))
            elif event_type == 'right_duplication_inversion':
                if bp2_band != bp3_band:
                    main_str = 'right-dup-inv({})({}{})'.format(path_chr, bp2_band, bp3_band)
                else:
                    main_str = 'right-dup-inv({})({})'.format(path_chr, bp2_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                sub_str1 = 'right duplication inversion on Chr{}: {}'.format(path_chr, chr_range)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                sub_str4, sub_str5 = report_cnv_genes_on_region(path_chr, bp2, bp3, '+1 CN')
                sub_bullets.append((sub_str1, sub_str2, sub_str3, sub_str4, sub_str5))
            elif event_type == 'insertion':
                # different report format if insertion is from different chr
                if 'Chr' + path_chr == bp2_chr:
                    # TODO: check ISCN syntax if bp2_band == bp3_band
                    main_str = 'ins({})({}{}{})'.format(path_chr, bp1_band, bp2_band, bp3_band)
                else:
                    main_str = 'ins({};{})({};{}{})'.format(path_chr, bp2_chr, bp1_band, bp2_band, bp3_band)
                chr_range = chr_range_tostr(bp2, bp3, bp2_band, bp3_band)
                sub_str1 = 'duplicated-insertion of Chr{}: {} into Chr{}: {}({})'.format(bp2_chr[3:], chr_range, path_chr, bp1, bp1_band)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(bp1_chr, bp1), (bp2_chr, bp2), (bp3_chr, bp3), (bp4_chr, bp4)])
                sub_str4, sub_str5 = report_cnv_genes_on_region(bp2_chr[3:], bp2, bp3, '+1 CN')
                sub_bullets.append((sub_str1, sub_str2, sub_str3, sub_str4, sub_str5))
            else:
                # continue
                raise RuntimeError('event not in allowed list')
            main_bullets.append(main_str)
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
                raise RuntimeError('non-terminal 2-break reciprocal translocation detected')
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
                main_bullets.append('t({};{})({};{})'.format(chr1, chr2, bp1_band, bp2_band))
                chr_range1 = chr_range_tostr(seg1_left_bp, seg1_right_bp, seg1_left_band, seg1_right_band)
                chr_range2 = chr_range_tostr(seg2_left_bp, seg2_right_bp, seg2_left_band, seg2_right_band)
                sub_str1 = 'balanced translocation between Chr{} and Chr{}, ' \
                           'between segments Chr{}: {} and Chr{}: {}'. \
                    format(chr1, chr2,
                           chr1, chr_range1,
                           chr2, chr_range2)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(chr1, seg1_left_bp), (chr1, seg1_right_bp),
                                                                           (chr2, seg2_left_bp), (chr2, seg2_right_bp)])
            else:
                main_bullets.append('t({};{})({};{})'.format(chr2, chr1, bp2_band, bp1_band))
                chr_range1 = chr_range_tostr(seg2_left_bp, seg2_right_bp, seg2_left_band, seg2_right_band)
                chr_range2 = chr_range_tostr(seg1_left_bp, seg1_right_bp, seg1_left_band, seg1_right_band)
                sub_str1 = 'balanced translocation between Chr{} and Chr{}, ' \
                           'between segments Chr{}: {} and Chr{}: {}'. \
                    format(chr2, chr1,
                           chr2, chr_range1,
                           chr1, chr_range2)
                sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(chr2, seg2_left_bp), (chr2, seg2_right_bp),
                                                                           (chr1, seg1_left_bp), (chr1, seg1_right_bp)])
            sub_bullets.append((sub_str1, sub_str2, sub_str3))
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
            main_bullets.append('ins-t({};{})({};{}{})'.format(ins_chr, event_seg_chr,
                                                               ins_site_left_band, event_seg_left_band, event_seg_right_band))
            chr_range = chr_range_tostr(event_seg_left_bp, event_seg_right_bp, event_seg_left_band, event_seg_right_band)
            sub_str1 = 'balanced non-reciprocal translocation of Chr{}: {} into Chr{}: {}({})'.\
                format(event_seg_chr, chr_range,
                       ins_chr, ins_site_left_bp, ins_site_left_band)
            sub_str2, sub_str3 = report_on_genes_based_on_breakpoints([(event_seg_chr, event_seg_left_bp),
                                                                       (event_seg_chr, event_seg_left_bp),
                                                                       (ins_chr, ins_site_left_bp)])
            sub_bullets.append((sub_str1, sub_str2, sub_str3))
        else:
            raise RuntimeError('illegal type assigned')
    return main_bullets, sub_bullets


def generate_latex_frontpage(title,
                             sample_names,
                             genefile_name,
                             breakpoint_reporting_proximity,
                             interpretation_insertion_threshold,
                             interpretation_deletion_threshold):
    output_str = "\\documentclass[12pt]{article}\n"
    output_str += "\\usepackage[letterpaper, margin=0.75in]{geometry}\n"
    output_str += "\\setcounter{secnumdepth}{0}\n"
    output_str += "\\usepackage{graphicx}\n"
    output_str += "\\usepackage{setspace}\n"
    output_str += "\\usepackage{titling}\n"
    output_str += "\\usepackage{enumitem}\n"
    output_str += "\\usepackage[T1]{fontenc}\n"
    output_str += "\\catcode`\_=12\n"
    output_str += "\\renewcommand\\maketitlehookc{\\vspace{-10ex}}\n"
    output_str += "\\usepackage{lipsum}\n"
    output_str += "\\graphicspath{ {./images/} }\n"
    output_str += "\\begin{document}\n"
    output_str += "\n"
    output_str += "\\title{" + title + "}\n"
    today = datetime.today()
    formatted_date = today.strftime("%b.%dth, %Y")
    output_str += "\\date{" + str(formatted_date) + "}\n"
    output_str += "\\maketitle\n"
    output_str += "\n"
    output_str += "\\textbf{{Samples Included: }} "
    output_str += ', '.join(sample_names) + "\n"
    output_str += "\n"
    output_str += "\\hfill\n"
    output_str += "\n"
    output_str += "\\textbf{Gene file used: } " + genefile_name + "\n"
    output_str += "\n"
    output_str += "\\hfill\n"
    output_str += "\n"
    output_str += "\\textbf{Breakpoint Gene Reporting Proximity: } " + str(breakpoint_reporting_proximity) + "kbp\n"
    output_str += "\n"
    output_str += "\\textbf{Threashold for event insertion size: } " + str(interpretation_insertion_threshold) + "kbp\n"
    output_str += "\n"
    output_str += "\\textbf{Threashold for event deletion size: } " + str(interpretation_deletion_threshold) + "kbp\n"
    output_str += "\n"
    output_str += "\\textbf{Supported SV types for interpretation: } \n"
    output_str += "\\begin{itemize}[leftmargin=3.5em,labelsep=0.5em]\n"
    output_str += "\\itemsep0em\n"
    output_str += "\\item deletion\n"
    output_str += "\\item inversion\n"
    output_str += "\\item single/repeated tandem duplication\n"
    output_str += "\\item left/right duplication inversion\n"
    output_str += "\\item 2/multiple-break reciprocal balanced translocation\n"
    output_str += "\\item nonreciprocal balanced translocation\n"
    output_str += "\\item duplicated insertion\n"
    output_str += "\\end{itemize}\n"
    output_str += "\n"
    output_str += "\\newpage\n"
    # output_str += "\\end{document}\n"

    return output_str

def batch_generate_latex_case_str(omkar_output_dir, image_dir):
    """
    :param omkar_output_dir: assume files are named as (int).txt
    :param image_dir: relative path to the latex DIR, assume exists an image file with the same name, but (int).pdf
    :return:
    """
    image_suffix = '.pdf'

    final_str = ""
    cases_with_events = []
    files = [file for file in os.listdir(omkar_output_dir)]
    files = sorted(files, key=int_file_keys)

    for file in files:
        # if file == '3.txt':
        if True:
            # if file in ['45.txt']:
            #     continue
            # if file in ['145.txt', '450.txt']:
            #     continue
            if file in ['54.txt', '205.txt']:
                continue
            filename = file.split('.')[0]
            file_path = omkar_output_dir + file
            print(file)
            mt_indexed_lists, mt_path_chrs, segment_dict, segment_size_dict = read_OMKar_to_indexed_list(file_path, forbidden_region_file)
            mt_path_chrs = [info.split(': ')[-1] for info in mt_path_chrs]
            wt_path_dict = generate_wt_from_OMKar_output(segment_dict)
            wt_indexed_lists = populate_wt_indexed_lists(mt_path_chrs, wt_path_dict)
            events, aligned_haplotypes = interpret_haplotypes(mt_indexed_lists, wt_indexed_lists, mt_path_chrs, segment_size_dict)
            main_bullets, sub_bullets = format_report(events, aligned_haplotypes, reverse_dict(segment_dict))

            if len(main_bullets) == 0:
                continue
            else:
                cases_with_events.append(filename)
            image_path = "{}/{}".format(image_dir, str(filename).zfill(3) + image_suffix)
            # image_path = "{}/{}".format(image_dir, str(filename) + image_suffix)
            if os.path.exists('latex_reports/' + image_path):
                # make sure the image file exists
                final_str += "\\section{{Sample Id: {}}}\n".format(filename)
                final_str += "\\begin{figure}[h!]\n"
                final_str += "\\centering\n"
                final_str += "\\includegraphics[width=3in]{{{}}}\n".format(image_path)
                final_str += "\\caption{{\\footnotesize Chromosomes with aberrant karyotypes}}\n"
                final_str += "\\label{{fig:karyotype_id{}}}\n".format(filename)
                final_str += "\\end{figure}\n"
                print('latex_reports/' + image_path, 'found')
            else:
                final_str += "\\section{{Sample Id: {}}}\n".format(filename)
                final_str += "\n"
                print('latex_reports/' + image_path, 'not found')
            final_str += "\\subsection{Events}\n"

            ## Iterate through all events
            c_event_id = 1  # 1-indexed



            for bullet_idx, main_bullet in enumerate(main_bullets):
                sub_bullet_list = sub_bullets[bullet_idx]
                final_str += "\\textbf{{Event {}: {}}}\n".format(bullet_idx, main_bullet)
                final_str += "\\begin{itemize}[leftmargin=3.5em,labelsep=0.5em]\n"
                for sub_bullet in sub_bullet_list:
                    if 'DDG2P' in sub_bullet and len(sub_bullet.split('\n\t')) > 1:
                        final_str += "\\begin{itemize}\n"
                        for ddg2p_str in sub_bullet.split('\n\t')[1:]:
                            final_str += "\\item {}\n".format(ddg2p_str)
                        final_str += "\\end{itemize}\n"
                    else:
                        final_str += "\\item {}\n".format(sub_bullet)
                final_str += "\\end{itemize}\n"
                final_str += "\n"
                final_str += "\\hfill"
                final_str += "\n"
                final_str += "\n"

            final_str += "\n"
            final_str += "\\newpage\n"
    final_str += "\n"
    final_str += "\\end{document}\n"
    return final_str, cases_with_events


def generate_latex_report(output_filename_prefix, front_page_str, batch_case_str):
    directory_path = os.path.dirname(output_filename_prefix) + '/'
    latex_path = output_filename_prefix + '.tex'
    with open(latex_path, 'w') as fp_write:
        fp_write.write(front_page_str)
        fp_write.write(batch_case_str)
    subprocess.run(['pdflatex', '-output-directory=' + directory_path, latex_path])


def chr_range_tostr(bpa, bpb, bpa_band, bpb_band):
    return "{}-{} ({} - {})".format(format(bpa, ',d'), format(bpb, ',d'), bpa_band, bpb_band)


def int_file_keys(f):
    return int(f.split('.')[0])


def test_latex(output_name):
    output_path = 'latex_reports/{}'.format(output_name)
    batch_case_str, cases_in_report = batch_generate_latex_case_str(data_dir, image_dir)
    front_str = generate_latex_frontpage('{} Data'.format(' '.join(output_name.split('_'))),
                                         cases_in_report,
                                         'hg38 all coding genes',
                                         50,
                                         200,
                                         200)
    generate_latex_report(output_path, front_str, batch_case_str)


if __name__ == "__main__":
    # test_interpreter()
    # test_segs_union()
    # test_reciprocal_trans()
    c_output_name = 'Dremsek'
    data_dir = '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_omkar/'
    # data_dir = '/media/zhaoyang-new/workspace/paul_dremsek/omkar_output/'
    # data_dir = '/media/zhaoyang-new/workspace/sunnyside/OMKar_output_paths/'
    # data_dir = '/media/zhaoyang-new/workspace/keyhole/OMKar_output_paths/'
    forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"
    # image_dir = '/media/zhaoyang-new/workspace/KarSim/KarComparator/latex_reports/paul_dremsek_plots/'
    # image_dir = '/media/zhaoyang-new/workspace/KarSim/KarComparator/latex_reports/sunnyside_plots/'
    image_dir = 'paul_dremsek_plots/'
    # batch_case_str = batch_generate_latex_case_str(data_dir, 'dremsek_images_2')
    test_latex(c_output_name)
