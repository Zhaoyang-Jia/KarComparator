import matplotlib.pyplot as plt
import matplotlib.patches as patches

# Example data for one chromosome
chromosome_data = {
    'name': 'chr1',
    'bands': [
        {'start': 0, 'end': 15, 'name': 'p15', 'color': 'lightgrey'},
        {'start': 15, 'end': 25, 'name': 'p14', 'color': 'grey'},
        {'start': 25, 'end': 35, 'name': 'p13', 'color': 'lightgrey'},
        {'start': 35, 'end': 45, 'name': 'p12', 'color': 'grey'},
        {'start': 45, 'end': 55, 'name': 'p11.2', 'color': 'red'},
        {'start': 55, 'end': 65, 'name': 'q11.2', 'color': 'blue'},
        {'start': 65, 'end': 75, 'name': 'q21', 'color': 'black'},
        {'start': 75, 'end': 85, 'name': 'q22', 'color': 'white'},
        {'start': 85, 'end': 95, 'name': 'q23', 'color': 'grey'},
        {'start': 95, 'end': 105, 'name': 'q24', 'color': 'white'},
        {'start': 105, 'end': 115, 'name': 'q25', 'color': 'grey'},
        {'start': 115, 'end': 125, 'name': 'q26', 'color': 'white'}
    ]
}


def plot_chromosome(ax, chromosome_data, y_offset=0):
    for band in chromosome_data['bands']:
        start = band['start']
        end = band['end']
        name = band['name']
        color = band['color']

        rect = patches.Rectangle((start, y_offset), end - start, 1, linewidth=1, edgecolor='black', facecolor=color)
        ax.add_patch(rect)
        ax.text((start + end) / 2, y_offset + 0.5, name, ha='center', va='center', fontsize=8, color='blue')

    ax.text(130, y_offset + 0.5, chromosome_data['name'], va='center', fontsize=12)
    ax.set_xlim(0, 140)
    ax.set_ylim(-1, 2)
    ax.axis('off')


fig, ax = plt.subplots(figsize=(10, 2))

plot_chromosome(ax, chromosome_data)

# Rotate the entire plot 90 degrees counterclockwise
plt.gca().set_aspect('equal', adjustable='box')
plt.gca().invert_yaxis()
plt.gca().invert_xaxis()
plt.gca().set_xlim(ax.get_xlim()[::-1])
plt.gca().set_ylim(ax.get_ylim()[::-1])
plt.gca().transData.rotate_deg(90)

plt.show()
