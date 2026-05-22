import pandas as pd
from Bio import Phylo
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.lines as mlines
import numpy as np
import os
import matplotlib as mpl

work_dir = '/mnt/d/hanjia/2.10/fig_tree_global_shenzhen'
tree_path = os.path.join(work_dir, 'iqtree_out.treefile')
meta_path = os.path.join(work_dir, 'metadata.csv')

# 1. Load tree
tree = Phylo.read(tree_path, "newick")
outgroup_node = next(tree.find_clades(name='NC_063383.1'))
tree.root_with_outgroup(outgroup_node)
tree.ladderize()

# 2. Load metadata
df = pd.read_csv(meta_path)
df['Country'] = df['Country'].replace('Dominican_Republic', 'Dominican')
meta_dict = df.set_index('seqName').to_dict('index')

# Calculate lineage frequencies to sort legend and assign colors
lineage_counts = df['Lineage'].value_counts()
sorted_lineages = sorted([l for l in lineage_counts.keys() if l != 'Unknown'])
top_lineages = sorted_lineages

# Assign colors for Country and Lineage
# Lineage colors
lineage_colors = {
    "C.1.1": "#96c37d",
    "E.3": "#e79a90",
    "C.1": "#efbc91",
    "F.2": "#fae5b8",
    "B.1": "#c82423",
    "E.1": "#d87659",
    "C.1.3": "#80b1d3",
    "E.4": "#bc80bd",
    "B.1.3": "#fccde5"
}

# Use a colormap but exclude grey colors (14 and 15 in tab20)
cmap_lin = plt.get_cmap('tab20')
allowed_colors = [mpl.colors.to_hex(cmap_lin(i)) for i in range(20) if i not in [14, 15]]
# Add tab20b or other colors if we need more
cmap_lin2 = plt.get_cmap('tab20b')
allowed_colors += [mpl.colors.to_hex(cmap_lin2(i)) for i in range(20) if i not in [14, 15]]

color_idx = 0
for lin in top_lineages:
    if lin not in lineage_colors:
        lineage_colors[lin] = allowed_colors[color_idx % len(allowed_colors)]
        color_idx += 1

# Country colors
country_counts = df['Country'].value_counts()
sorted_countries = sorted([c for c in country_counts.keys() if c != 'Unknown'])
top_countries = sorted_countries

country_colors = {}
# Use tab20/Dark2 manually picked colors to avoid anything too light/bright (like yellow/light green)
allowed_c_colors = [
    '#1f77b4', '#ff7f0e', '#2ca02c', '#9467bd', '#8c564b',
    '#e377c2', '#7f7f7f', '#bcbd22', '#17becf', '#393b79',
    '#8c6d31', '#843c39', '#7b4173', '#5254a3', '#637939',
    '#000080', '#006400', '#8b0000', '#8b008b', '#ff8c00',
    '#483d8b', '#2f4f4f', '#008080', '#b8860b', '#a0522d',
    '#d2691e', '#cd5c5c', '#4682b4', '#556b2f', '#9932cc'
]

c_idx = 0
for c in top_countries:
    if c == 'China':
        country_colors[c] = '#d62728' # A deeper, clear red for China
    elif c == 'Reference':
        country_colors[c] = '#000000'
    else:
        country_colors[c] = allowed_c_colors[c_idx % len(allowed_c_colors)]
        c_idx += 1

# 3. Calculate coordinates (rectangular)
y_coords = {}
terminals = tree.get_terminals()
for i, term in enumerate(terminals):
    y_coords[term] = i

for clade in tree.get_nonterminals(order='postorder'):
    if clade.clades:
        y_coords[clade] = sum(y_coords[child] for child in clade.clades) / len(clade.clades)
    else:
        y_coords[clade] = 0

def get_depths(tree):
    depths = {tree.root: 0}
    for clade in tree.find_clades(order='preorder'):
        for child in clade:
            branch_len = child.branch_length if child.branch_length is not None else 0
            depths[child] = depths[clade] + branch_len
    return depths

x_coords = get_depths(tree)

# Convert to circular coordinates
def to_polar(x, y, y_max, max_radius=10):
    # theta from 0 to 2*pi
    # leave a small gap so it doesn't overlap completely
    theta = 2 * np.pi * (y / y_max)
    r = x
    return theta, r

y_max = max(y_coords.values())
max_x = max(x_coords.values())

fig = plt.figure(figsize=(24, 24))
ax = fig.add_subplot(111, polar=True)
ax.set_theta_zero_location("N")
ax.set_theta_direction(-1)

# Draw tree
def draw_clade_polar(clade):
    theta_parent, r_parent = to_polar(x_coords[clade], y_coords[clade], y_max)
    if clade.clades:
        # vertical arc
        y_min = min(y_coords[c] for c in clade.clades)
        y_max_c = max(y_coords[c] for c in clade.clades)
        theta_min, _ = to_polar(x_coords[clade], y_min, y_max)
        theta_max, _ = to_polar(x_coords[clade], y_max_c, y_max)
        
        thetas = np.linspace(theta_min, theta_max, 100)
        rs = np.full_like(thetas, r_parent)
        # Increase line width for internal branches
        ax.plot(thetas, rs, color='black', lw=1.0)
        
        for child in clade.clades:
            theta_child, r_child = to_polar(x_coords[child], y_coords[child], y_max)
            # horizontal line
            ax.plot([theta_child, theta_child], [r_parent, r_child], color='black', lw=1.0)
            draw_clade_polar(child)

draw_clade_polar(tree.root)

# Draw rings and tips
# Bring the rings closer to the tree to reduce the gap
ring1_r = max_x * 1.01
ring2_r = max_x * 1.05
ring_width = max_x * 0.035

for term in terminals:
    name = term.name
    theta, r = to_polar(x_coords[term], y_coords[term], y_max)
    
    info = meta_dict.get(name, {'Country': 'Unknown', 'Lineage': 'Unknown', 'is_shenzhen': False})
    
    # Shenzhen node in red
    if info['is_shenzhen']:
        ax.plot(theta, r, marker='o', color='red', markersize=4, zorder=10)
    
    # Ring 1: Country
    c_color = country_colors.get(info['Country'], '#cccccc')
    ax.bar(theta, ring_width, bottom=ring1_r, width=2*np.pi/y_max, color=c_color, edgecolor='none', align='center')
    
    # Ring 2: Lineage
    l_color = lineage_colors.get(info['Lineage'], '#cccccc')
    ax.bar(theta, ring_width, bottom=ring2_r, width=2*np.pi/y_max, color=l_color, edgecolor='none', align='center')

# Remove axes
ax.axis('off')

# Configure font sizes
LEGEND_TITLE_FONT_SIZE = 30
LEGEND_LABEL_FONT_SIZE = 22

# Add legends
# Lineage legend
# Show all lineages in legend, formatted in multiple columns
handles_lin = [mlines.Line2D([], [], color='w', marker='s', markerfacecolor=lineage_colors[l], markersize=15, label=f"{l}") 
               for l in top_lineages]
leg1 = ax.legend(handles=handles_lin, title="Lineage", loc='center', bbox_to_anchor=(0.38, 0.5), frameon=False, fontsize=LEGEND_LABEL_FONT_SIZE, title_fontsize=LEGEND_TITLE_FONT_SIZE, ncol=3, columnspacing=1.0, handletextpad=0.5)
ax.add_artist(leg1)

# Country legend
handles_c = [mlines.Line2D([], [], color='w', marker='s', markerfacecolor=country_colors[c], markersize=15, label=f"{c}") 
               for c in top_countries]
leg2 = ax.legend(handles=handles_c, title="Country", loc='center', bbox_to_anchor=(0.64, 0.5), frameon=False, fontsize=LEGEND_LABEL_FONT_SIZE, title_fontsize=LEGEND_TITLE_FONT_SIZE, ncol=2, columnspacing=1.0, handletextpad=0.5)
ax.add_artist(leg2)

plt.tight_layout()
plt.savefig(os.path.join(work_dir, 'combined_tree.pdf'), bbox_inches='tight')
plt.savefig(os.path.join(work_dir, 'combined_tree.svg'), bbox_inches='tight')
plt.savefig(os.path.join(work_dir, 'combined_tree.png'), bbox_inches='tight', dpi=500)
plt.savefig(os.path.join(work_dir, 'combined_tree.jpg'), bbox_inches='tight', dpi=500)
print("Tree successfully drawn and saved.")
