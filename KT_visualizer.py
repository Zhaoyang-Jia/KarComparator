import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as mcolors
import colorsys


def reduce_saturation(color, factor):
    """Reduce the saturation of a color by a given factor (0-1)."""
    rgb = mcolors.to_rgb(color)
    h, l, s = colorsys.rgb_to_hls(*rgb)
    s *= factor
    return colorsys.hls_to_rgb(h, l, s)


def get_text_color(bg_color):
    """Determine if black or white text would be more readable on the given background color."""
    rgb = mcolors.to_rgb(bg_color)
    luminance = 0.2126 * rgb[0] + 0.7152 * rgb[1] + 0.0722 * rgb[2]
    return 'black' if luminance > 0.5 else 'white'


# constant parameters, x and y based on the unrotated orientation
WHOLE_CHR_Y_OFFSET = 2
CHR_HEADER_Y_OFFSET = -0.2 + WHOLE_CHR_Y_OFFSET
CHR_BAND_Y_OFFSET = WHOLE_CHR_Y_OFFSET
CHR_HEADER_HIGHLIGHT_Y_OFFSET = -0.63 + WHOLE_CHR_Y_OFFSET
ORIGIN_Y_OFFSET = 1.3 + WHOLE_CHR_Y_OFFSET
TICK_Y_OFFSET = -0.05 + WHOLE_CHR_Y_OFFSET
TICK_MARKING_Y_OFFSET = -0.3 + WHOLE_CHR_Y_OFFSET
LABEL_BAR_Y_OFFSET = -0.05 + WHOLE_CHR_Y_OFFSET
LABEL_MARK_Y_OFFSET = -2.05 + WHOLE_CHR_Y_OFFSET

CHR_HEADER_X_OFFSET = 5
CHR_HEADER_HIGHLIGHT_COLOR = 'red'
BAND_WIDTH = 1
ORIGIN_WIDTH = 0.5

BAND_SATURATION = 0.7
BAND_ALPHA = 0.5
BAND_TEXT_WEIGHT = 'normal'
ORIGIN_ALPHA = 0.7

TICK_MARKING_X_OFFSET = 0.25  # helps to center the tickmarkings
TICK_LEN = 0.2
TICK_THICKNESS = 0.9
TICK_ALPHA = 0.65
TICK_MARKING_ALPHA = 0.85
SUBTICK_LEN = 0.065
SUBTICK_THICKNESS = 0.9
SUBTICK_ALPHA = 0.5

LABEL_BAR_LEN = 0.4
LABEL_BAR_THICKNESS = 1.5
LABEL_BAR_ALPHA = 1
LABEL_MARK_ALPHA = 1

BAND_FONTSIZE = 8
LABEL_MARK_FONTSIZE = 8
CHR_HEADER_HIGHLIGHT_FONTSIZE = 18
CHR_HEADER_FONTSIZE = 11
SCALE_MARKING_FONTSIZE = 7
LABEL_MARK_COLOR = 'red'

# constant parameters: do not adjust, dependent to above
TICK_END_Y_OFFSET = TICK_Y_OFFSET - TICK_LEN
SUBTICK_END_Y_OFFSET = TICK_Y_OFFSET - SUBTICK_LEN
LABEL_BAR_END_Y_OFFSET = LABEL_BAR_Y_OFFSET - LABEL_BAR_LEN


color_mapping = {
    'gneg': 'white',
    'gpos25': '#C0C0C0',  # light grey
    'gpos50': '#808080',  # medium grey
    'gpos75': '#404040',  # dark grey
    'gpos100': 'black',   # full black
    'acen': 'red',        # centromere
    'gvar': 'blue',       # variable region
    'stalk': '#87CEEB'    # light blue (skyblue)
}
chr_color_mapping = {
    '1': '#73a9a3',
    '2': '#ffffd9',
    '3': '#a5a4c6',
    '4': '#fb8072',
    '5': '#73a9c3',
    '6': '#d6a780',
    '7': '#a8d8b4',
    '8': '#fcebf5',
    '9': '#d9d9d9',
    '10': '#b780bd',
    '11': '#c5e7c5',
    '12': '#ffed9f',
    '13': '#ff6a6a',
    '14': '#73a8d3',
    '15': '#78af78',
    '16': '#a94ea3',
    '17': '#ff7f40',
    '18': '#ffff99',
    '19': '#a56e40',
    '20': '#f7a1bf',
    '21': '#999999',
    '22': '#66c2a5',
    'X': '#fc8d62',
    'Y': '#8da0cb'
}


# Reduce saturation of the colors in the color mapping
reduced_saturation_mapping = {k: reduce_saturation(v, BAND_SATURATION) for k, v in color_mapping.items()}


def plot_chromosome(ax, chromosome_data, y_offset, x_offset):
    ## Chrom header
    ax.text(x_offset, y_offset + CHR_HEADER_Y_OFFSET, chromosome_data['name'],
            va='bottom', fontsize=CHR_HEADER_FONTSIZE, rotation=90, weight='bold')

    ## Bands
    for band in chromosome_data['bands']:
        start = band['start']
        end = band['end']
        name = band['name']
        stain = band['stain']
        color = reduced_saturation_mapping[stain]
        text_color = get_text_color(color)

        chrom_bands = patches.Rectangle((x_offset + start + CHR_HEADER_X_OFFSET, y_offset + CHR_BAND_Y_OFFSET), end - start, BAND_WIDTH,
                                        linewidth=1, edgecolor='black', facecolor=color, alpha=BAND_ALPHA)
        ax.add_patch(chrom_bands)
        ax.text(x_offset + (start + end) / 2 + CHR_HEADER_X_OFFSET, y_offset + BAND_WIDTH / 2 + CHR_BAND_Y_OFFSET, name,
                ha='center', va='center', fontsize=BAND_FONTSIZE, color=text_color, rotation=90, weight=BAND_TEXT_WEIGHT)

    ## Origins
    for origin in chromosome_data['origins']:
        start = origin['start']
        end = origin['end']
        name = origin['name']
        color = chr_color_mapping[name]
        text_color = get_text_color(color)
        origin_bands = patches.Rectangle((x_offset + start + CHR_HEADER_X_OFFSET, y_offset + ORIGIN_Y_OFFSET), end - start, ORIGIN_WIDTH,
                                         linewidth=1, edgecolor='black', facecolor=color, alpha=ORIGIN_ALPHA)
        ax.add_patch(origin_bands)
        ax.text(x_offset + (start + end) / 2 + CHR_HEADER_X_OFFSET, y_offset + ORIGIN_Y_OFFSET + ORIGIN_WIDTH / 2, name,
                ha='center', va='center', fontsize=8, color=text_color, rotation=90, weight=BAND_TEXT_WEIGHT)

    ## Modified Chrom's header
    if chromosome_data['highlight']:
        ax.text(x_offset, y_offset + CHR_HEADER_HIGHLIGHT_Y_OFFSET, "*",
                va='bottom', ha='left',
                fontsize=CHR_HEADER_HIGHLIGHT_FONTSIZE, rotation=90, weight='bold', color=CHR_HEADER_HIGHLIGHT_COLOR)

    ## Add sub-scale ticks
    for i in range(0, chromosome_data['length'] + 1, 2):
        subtick_x_loc = x_offset + i + CHR_HEADER_X_OFFSET
        subtick_y_start_loc = y_offset + TICK_Y_OFFSET
        subtick_y_end_loc = y_offset + SUBTICK_END_Y_OFFSET
        ax.plot([subtick_x_loc, subtick_x_loc], [subtick_y_start_loc, subtick_y_end_loc],
                color='grey', linewidth=SUBTICK_THICKNESS, alpha=SUBTICK_ALPHA)
    ## Add scale ticks
    for i in range(0, chromosome_data['length'] + 1, 10):
        tick_x_loc = x_offset + i + CHR_HEADER_X_OFFSET
        tick_y_start_loc = y_offset + TICK_Y_OFFSET - SUBTICK_LEN  # to not overlap with the subticks
        tick_y_end_loc = y_offset + TICK_END_Y_OFFSET
        tickmark_x_loc = x_offset + i + CHR_HEADER_X_OFFSET + TICK_MARKING_X_OFFSET
        tickmark_y_loc = y_offset + TICK_MARKING_Y_OFFSET
        ax.plot([tick_x_loc, tick_x_loc], [tick_y_start_loc, tick_y_end_loc],
                color='red', linewidth=TICK_THICKNESS, alpha=TICK_ALPHA)
        ax.text(tickmark_x_loc, tickmark_y_loc, str(i),
                ha='center', va='top', fontsize=SCALE_MARKING_FONTSIZE, rotation=90, alpha=TICK_MARKING_ALPHA)

    ## sv_labels
    for sv_label in chromosome_data['sv_labels']:
        pos = sv_label['pos']
        label = sv_label['label']
        label_bar_x_loc = x_offset + pos + CHR_HEADER_X_OFFSET
        label_bar_y_start_loc = y_offset + LABEL_BAR_Y_OFFSET
        label_bar_y_end_loc = y_offset + LABEL_BAR_END_Y_OFFSET - SUBTICK_LEN
        label_mark_y_loc = label_bar_y_end_loc + LABEL_MARK_Y_OFFSET
        ax.plot([label_bar_x_loc, label_bar_x_loc], [label_bar_y_start_loc, label_bar_y_end_loc],
                color=LABEL_MARK_COLOR, linewidth=LABEL_BAR_THICKNESS, alpha=LABEL_BAR_ALPHA)
        ax.text(label_bar_x_loc, label_mark_y_loc, label,
                ha='center', va='top', fontsize=LABEL_MARK_FONTSIZE, rotation=90,
                alpha=LABEL_MARK_ALPHA, color=LABEL_MARK_COLOR, weight='normal')

    ## Limit chrom plot size
    # ax.set_xlim(0, chromosome_data['length']+CHR_HEADER_X_OFFSET+10)
    ax.set_xlim(0, chromosome_data['length'] + CHR_HEADER_X_OFFSET + 6)
    ax.set_ylim(0, len(chromosomes_data) * 4)
    ax.axis('off')


# Example data for multiple chromosomes## Chrom header
#     ax.text(x_offset, y_offset + CHR_HEADER_Y_OFFSET, chromosome_data['name'], va='bottom', fontsize=CHR_HEADER_FONTSIZE, rotation=90)
chromosomes_data = [
    {
        'name': 'Chr1',
        'length': 125,
        'bands': [
            {'start': 0, 'end': 15, 'name': 'p15', 'stain': 'gneg'},
            {'start': 15, 'end': 25, 'name': 'p14', 'stain': 'gpos25'},
            {'start': 25, 'end': 35, 'name': 'p13', 'stain': 'gpos50'},
            {'start': 35, 'end': 45, 'name': 'p12', 'stain': 'gpos75'},
            {'start': 45, 'end': 55, 'name': 'p11.2', 'stain': 'acen'},
            {'start': 55, 'end': 65, 'name': 'q11.2', 'stain': 'acen'},
            {'start': 65, 'end': 75, 'name': 'q21', 'stain': 'gpos100'},
            {'start': 75, 'end': 85, 'name': 'q22', 'stain': 'gneg'},
            {'start': 85, 'end': 95, 'name': 'q23', 'stain': 'gvar'},
            {'start': 95, 'end': 105, 'name': 'q24', 'stain': 'gneg'},
            {'start': 105, 'end': 115, 'name': 'q25', 'stain': 'stalk'},
            {'start': 115, 'end': 125, 'name': 'q26', 'stain': 'gneg'}
        ],
        'origins': [
            {'start': 0, 'end': 95, 'name': '1'},
            {'start': 95, 'end': 125, 'name': '2'}
        ],
        'highlight': False,
        'sv_labels': [
            {'pos': 24, 'label': '[2]INS'}
        ]
    },
    {
        'name': 'Chr2',
        'length': 120,
        'bands': [
            {'start': 0, 'end': 10, 'name': 'p15', 'stain': 'gneg'},
            {'start': 10, 'end': 20, 'name': 'p14', 'stain': 'gpos25'},
            {'start': 20, 'end': 30, 'name': 'p13', 'stain': 'gpos50'},
            {'start': 30, 'end': 40, 'name': 'p12', 'stain': 'gpos75'},
            {'start': 40, 'end': 50, 'name': 'p11.2', 'stain': 'acen'},
            {'start': 50, 'end': 60, 'name': 'q11.2', 'stain': 'acen'},
            {'start': 60, 'end': 70, 'name': 'q21', 'stain': 'gpos100'},
            {'start': 70, 'end': 80, 'name': 'q22', 'stain': 'gneg'},
            {'start': 80, 'end': 90, 'name': 'q23', 'stain': 'gvar'},
            {'start': 90, 'end': 100, 'name': 'q24', 'stain': 'gneg'},
            {'start': 100, 'end': 110, 'name': 'q25', 'stain': 'stalk'},
            {'start': 110, 'end': 120, 'name': 'q26', 'stain': 'gneg'}
        ],
        'origins': [
            {'start': 0, 'end': 60, 'name': '3'},
            {'start': 60, 'end': 100, 'name': '5'},
            {'start': 100, 'end': 120, 'name': '1'}
        ],
        'highlight': True,
        'sv_labels': [
            {'pos': 58, 'label': '[2]DEL'},
            {'pos': 79, 'label': '[2]DUPINV'}
        ]
    },
    {
        'name': 'Chr2',
        'length': 120,
        'bands': [
            {'start': 0, 'end': 10, 'name': 'p15', 'stain': 'gneg'},
            {'start': 10, 'end': 20, 'name': 'p14', 'stain': 'gpos25'},
            {'start': 20, 'end': 30, 'name': 'p13', 'stain': 'gpos50'},
            {'start': 30, 'end': 40, 'name': 'p12', 'stain': 'gpos75'},
            {'start': 40, 'end': 50, 'name': 'p11.2', 'stain': 'acen'},
            {'start': 50, 'end': 60, 'name': 'q11.2', 'stain': 'acen'},
            {'start': 60, 'end': 70, 'name': 'q21', 'stain': 'gpos100'},
            {'start': 70, 'end': 80, 'name': 'q22', 'stain': 'gneg'},
            {'start': 80, 'end': 90, 'name': 'q23', 'stain': 'gvar'},
            {'start': 90, 'end': 100, 'name': 'q24', 'stain': 'gneg'},
            {'start': 100, 'end': 110, 'name': 'q25', 'stain': 'stalk'},
            {'start': 110, 'end': 120, 'name': 'q26', 'stain': 'gneg'}
        ],
        'origins': [
            {'start': 0, 'end': 60, 'name': '3'},
            {'start': 60, 'end': 100, 'name': '5'},
            {'start': 100, 'end': 120, 'name': '1'}
        ],
        'highlight': True,
        'sv_labels': [
            {'pos': 58, 'label': '[2]DEL'},
            {'pos': 79, 'label': '[2]T'}
        ]
    },
{
        'name': 'Chr2',
        'length': 120,
        'bands': [
            {'start': 0, 'end': 10, 'name': 'p15', 'stain': 'gneg'},
            {'start': 10, 'end': 20, 'name': 'p14', 'stain': 'gpos25'},
            {'start': 20, 'end': 30, 'name': 'p13', 'stain': 'gpos50'},
            {'start': 30, 'end': 40, 'name': 'p12', 'stain': 'gpos75'},
            {'start': 40, 'end': 50, 'name': 'p11.2', 'stain': 'acen'},
            {'start': 50, 'end': 60, 'name': 'q11.2', 'stain': 'acen'},
            {'start': 60, 'end': 70, 'name': 'q21', 'stain': 'gpos100'},
            {'start': 70, 'end': 80, 'name': 'q22', 'stain': 'gneg'},
            {'start': 80, 'end': 90, 'name': 'q23', 'stain': 'gvar'},
            {'start': 90, 'end': 100, 'name': 'q24', 'stain': 'gneg'},
            {'start': 100, 'end': 110, 'name': 'q25', 'stain': 'stalk'},
            {'start': 110, 'end': 120, 'name': 'q26', 'stain': 'gneg'}
        ],
        'origins': [
            {'start': 0, 'end': 60, 'name': '3'},
            {'start': 60, 'end': 100, 'name': '5'},
            {'start': 100, 'end': 120, 'name': '1'}
        ],
        'highlight': True,
        'sv_labels': [
            {'pos': 58, 'label': '[2]DEL'},
            {'pos': 79, 'label': '[2]DUP'}
        ]
    }
]


if __name__ == '__main__':
    plt.rcParams['figure.dpi'] = 250
    n_chrom = len(chromosomes_data)
    n_row = n_chrom // 4 + 1
    fig, i_ax = plt.subplots(figsize=(5 * n_row, 2 * min(4, n_chrom)))

    for chrom_idx, i_chromosome_data in enumerate(chromosomes_data):
        row = chrom_idx // 4
        col = chrom_idx % 4
        plot_chromosome(i_ax, i_chromosome_data, y_offset=col * 4, x_offset=row * 28)

    # plt.show(bbox_inches='tight')
    plt.savefig('test_fig.png', bbox_inches='tight')

    from PIL import Image

    def rotate_image(input_image_path, output_image_path):
        # Open an image file
        with Image.open(input_image_path) as img:
            # Rotate the image by 90 degrees
            rotated_img = img.rotate(270, expand=True)
            # Save the rotated image
            rotated_img.save(output_image_path)


    rotate_image('test_fig.png', 'test_fig_rotated.png')

