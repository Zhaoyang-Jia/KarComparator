from KT_visualizer import *
from jinja2 import Environment, FileSystemLoader
import os
import base64
import re

def image_to_base64(image_path):
    try:
        with open(image_path, "rb") as img_file:
            return base64.b64encode(img_file.read()).decode('utf-8')
    except FileNotFoundError:
        print(f"Error: File {image_path} not found.")
        return ""


def batch_populate_contents(omkar_output_dir, image_dir, file_of_interest=None, compile_image=False):
    cases_with_events = []
    image_paths = []
    iscn_reports = []
    genes_reports = []
    files = [file for file in os.listdir(omkar_output_dir)]
    for file in files:
        if file_of_interest is not None:
            if file not in file_of_interest:
                continue
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
        else:
            cases_with_events.append(filename)
        dependent_clusters, cluster_events = form_dependent_clusters(events, aligned_haplotypes)
        print(dependent_clusters)
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
            relative_image_path = image_dir.replace('latex_reports/', '') + image_path.split('/')[-1]
            if compile_image:
                if len(c_vis_input) <= 4:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_VERTICAL_SPLIT)
                else:
                    make_image(c_vis_input, max_chr_length(c_vis_input), image_prefix, IMG_LENGTH_SCALE_HORIZONTAL_SPLIT)

            image_paths.append(relative_image_path)
            iscn_reports.append(iscn_events)
            genes_reports.append(genes_report)

    return cases_with_events, image_paths, iscn_reports, genes_reports


def html_hyperlink_coordinates(input_str, proximity=50000):
    return_dict = {}  # {replacement_string: hyperlinked_string}

    pattern = r'Chr(\d+|X|Y): (\d{1,3}(?:,\d{3})*)-(\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        start_pos = int(match.group(2).replace(',', ''))
        end_pos = int(match.group(3).replace(',', ''))
        ucsc_url = get_ucsc_url('chr' + match.group(1), start_pos, end_pos)
        hyperlinked_str = f'<a href="{ucsc_url}">{replacement_str}</a>'
        return_dict[replacement_str] = hyperlinked_str

    pattern = r'Chr(\d+): (\d{1,3}(?:,\d{3})*) \(.*?\)'
    matches_itr = re.finditer(pattern, input_str)
    for match in matches_itr:
        replacement_str = input_str[match.start(): match.end()]
        chrom = 'chr' + match.group(1)
        pos = int(match.group(2).replace(',', ''))
        c_chr_length = get_chr_length_from_forbidden_file(chrom)
        ucsc_url = get_ucsc_url(chrom, max(0, pos - proximity), min(c_chr_length, pos + proximity))
        hyperlinked_str = f'<a href="{ucsc_url}">{replacement_str}</a>'
        return_dict[replacement_str] = hyperlinked_str

    return return_dict


def hyperlink_iscn_interpretation(input_str):
    hyperlinked_mapping = html_hyperlink_coordinates(input_str)
    for replacement_str, hyperlinked_str in hyperlinked_mapping.items():
        input_str = input_str.replace(replacement_str, hyperlinked_str)
    return input_str


def test(compile_image, cases_of_interest):
    omkar_output_dir = 'real_case_data/sunnyside_OMKar_output_paths/'
    image_output_dir = 'html_reports/test_HTML_plots/'
    os.makedirs(image_output_dir, exist_ok=True)

    title = 'Sunnyside'
    cases_with_events, image_paths, iscn_reports, genes_reports = batch_populate_contents(omkar_output_dir, image_output_dir,
                                                                                          file_of_interest=cases_of_interest, compile_image=compile_image)
    images_base64 = [image_to_base64(img) for img in image_paths]

    formatted_genes_reports = [format_genes_report(genes_report) for genes_report in genes_reports]
    columns_order = ['SV', 'rationale', 'gene name', 'gene omim']

    ## hyperlinking
    for iscn_report in iscn_reports:
        for iscn_report_idx, (iscn, sv_interpretation) in enumerate(iscn_report):
            hyperlinked_sv_interpretation = hyperlink_iscn_interpretation(sv_interpretation)
            iscn_report[iscn_report_idx][1] = hyperlinked_sv_interpretation

    content = [(text, image, table_content) for text, image, table_content in zip(iscn_reports, images_base64, formatted_genes_reports)]

    env = Environment(loader=FileSystemLoader('html_reports/'))
    template = env.get_template('template.html')
    rendered_html = template.render(title=title, content=content, columns_order=columns_order)

    output_file = 'html_reports/test.html'
    with open(output_file, 'w') as f:
        f.write(rendered_html)
    print(f"HTML file generated: {os.path.abspath(output_file)}")


def manual_test():
    # Define the data
    title = "My Text and Images"
    texts = ['text1',
            'text2',
            'text3']

    ## ZJ: image paths need to be relative path
    images = ['/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster0_rotated.png',
              '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster1_rotated.png',
              '/Users/zhaoyangjia/PyCharm_Repos/KarComparator/latex_reports/paul_dremsek_plots_new/3_imagecluster2_rotated.png']
    images_base64 = []
    for img in images:
        images_base64.append(image_to_base64(img))


    table_contents = [
        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]],

        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]],

        [["Row1-Col1", "Row1-Col2"], ["Row2-Col1", "Row2-Col2"], ["Row3-Col1", "Row3-Col2"], ["Row4-Col1", "Row4-Col2"],
         ["Row5-Col1", "Row5-Col2"], ["Row6-Col1", "Row6-Col2"], ["Row7-Col1", "Row7-Col2"], ["Row8-Col1", "Row8-Col2"],
         ["Row9-Col1", "Row9-Col2"], ["Row10-Col1", "Row10-Col2"], ["Row11-Col1", "Row11-Col2"], ["Row12-Col1", "Row12-Col2"]]
    ]

    content = [(text, image, table_content) for text, image, table_content in zip(texts, images_base64, table_contents)]

    # Create an environment for Jinja2
    env = Environment(loader=FileSystemLoader('html_reports/'))
    template = env.get_template('template.html')

    # Render the template with the data
    rendered_html = template.render(title=title, content=content)

    # Write the rendered HTML to a file
    output_file = 'html_reports/test.html'
    with open(output_file, 'w') as f:
        f.write(rendered_html)

    print(f"HTML file generated: {os.path.abspath(output_file)}")


if __name__ == "__main__":
    forbidden_region_file = "Metadata/acrocentric_telo_cen.bed"
    test(True, ['130.txt', '205.txt'])
