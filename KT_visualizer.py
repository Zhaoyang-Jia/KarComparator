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


# constant parameters, x and y based on the unrotated orientation
CHR_HEADER_X_OFFSET = 13
CHR_WIDTH = 0.75
ORIGIN_Y_OFFSET = 1
ORIGIN_WIDTH = 0.4

BAND_SATURATION = 0.7
BAND_ALPHA = 0.5

TICK_Y_OFFSET = -0.05
TICK_MARKING_Y_OFFSET = -0.3
TICK_MARKING_X_OFFSET = 0.25  # helps to center the tickmarkings
TICK_LEN = 0.2
TICK_THICKNESS = 0.9
TICK_ALPHA = 0.75
TICK_MARKING_ALPHA = 0.85
SUBTICK_LEN = 0.065
SUBTICK_THICKNESS = 0.9
SUBTICK_ALPHA = 0.75

# constant parameters: do not adjust, dependent to above
TICK_END_Y_OFFSET = TICK_Y_OFFSET - TICK_LEN
SUBTICK_END_Y_OFFSET = TICK_Y_OFFSET - SUBTICK_LEN


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

# Reduce saturation of the colors in the color mapping
reduced_saturation_mapping = {k: reduce_saturation(v, BAND_SATURATION) for k, v in color_mapping.items()}


def plot_chromosome(ax, chromosome_data, y_offset=0):
    ## Bands
    for band in chromosome_data['bands']:
        start = band['start']
        end = band['end']
        name = band['name']
        stain = band['stain']
        color = reduced_saturation_mapping[stain]

        chrom_bands = patches.Rectangle((start + CHR_HEADER_X_OFFSET, y_offset), end - start, CHR_WIDTH,
                                        linewidth=1, edgecolor='black', facecolor=color, alpha=BAND_ALPHA)
        ax.add_patch(chrom_bands)
        ax.text((start + end) / 2 + CHR_HEADER_X_OFFSET, y_offset + CHR_WIDTH / 2, name,
                ha='center', va='center', fontsize=8, color='blue', rotation=90)

        origin_bands = patches.Rectangle((start + CHR_HEADER_X_OFFSET, y_offset + ORIGIN_Y_OFFSET), end - start, ORIGIN_WIDTH,
                                         linewidth=1, edgecolor='black', facecolor=color)
        ax.add_patch(origin_bands)

    ## Chrom header
    ax.text(0, y_offset, chromosome_data['name'], va='center', fontsize=12, rotation=90)

    ## Add sub-scale ticks
    for i in range(0, chromosome_data['length'] + 1, 1):
        subtick_x_loc = i + CHR_HEADER_X_OFFSET
        subtick_y_start_loc = y_offset + TICK_Y_OFFSET
        subtick_y_end_loc = y_offset + SUBTICK_END_Y_OFFSET
        ax.plot([subtick_x_loc, subtick_x_loc], [subtick_y_start_loc, subtick_y_end_loc],
                color='grey', linewidth=SUBTICK_THICKNESS, alpha=SUBTICK_ALPHA)
    ## Add scale ticks
    for i in range(0, chromosome_data['length'] + 1, 10):
        tick_x_loc = i + CHR_HEADER_X_OFFSET
        tick_y_start_loc = y_offset + TICK_Y_OFFSET - SUBTICK_LEN  # to not overlap with the subticks
        tick_y_end_loc = y_offset + TICK_END_Y_OFFSET
        tickmark_x_loc = i + CHR_HEADER_X_OFFSET + TICK_MARKING_X_OFFSET
        tickmark_y_loc = y_offset + TICK_MARKING_Y_OFFSET
        ax.plot([tick_x_loc, tick_x_loc], [tick_y_start_loc, tick_y_end_loc],
                color='red', linewidth=TICK_THICKNESS, alpha=TICK_ALPHA)
        ax.text(tickmark_x_loc, tickmark_y_loc, str(i),
                ha='center', va='top', fontsize=8, rotation=90, alpha=TICK_MARKING_ALPHA)

    ## Limit chrom plot size
    ax.set_xlim(0, chromosome_data['length'] + CHR_HEADER_X_OFFSET+10)  # TODO: scale x based on max-length chrom
    ax.set_ylim(-2, len(chromosomes_data) * 4)
    ax.axis('off')


# Example data for multiple chromosomes
chromosomes_data = [
    {
        'name': '1',
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
        ]
    },
    {
        'name': '2',
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
        ]
    }
]


if __name__ == '__main__':
    fig, i_ax = plt.subplots(figsize=(10, 10))  # TODO: scale y based on n-chrom, scale x based on max-length chrom

    for chrom_idx, i_chromosome_data in enumerate(chromosomes_data):
        plot_chromosome(i_ax, i_chromosome_data, y_offset=chrom_idx * 4)

    plt.show()
