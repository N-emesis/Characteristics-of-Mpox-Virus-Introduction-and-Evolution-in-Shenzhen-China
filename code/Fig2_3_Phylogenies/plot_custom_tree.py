import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.dates as mdates
from datetime import datetime
import pandas as pd
from Bio import Phylo
import os
import numpy as np

# Files
rename_tsv = "/mnt/d/hanjia/2.10/rename/197_shenzhen_rename.tsv"
nextclade_csv = "/mnt/d/hanjia/2.10/dnds/data/nextclade_shenzhen/nextclade.csv"
date_config = "/mnt/d/hanjia/2.10/tree_time_new/itol_text_config.txt"
nwk_path = "/mnt/d/hanjia/2.10/tree_time_new/tree.nwk"
output_dir = os.path.dirname(nwk_path)

# 1. Parse Nextclade Lineages
rename_df = pd.read_csv(rename_tsv, sep='\t', header=0)
# Original name -> New name
orig_to_new = dict(zip(rename_df.iloc[:, 0], rename_df.iloc[:, 1]))

nc_df = pd.read_csv(nextclade_csv, sep=';', low_memory=False)
# 'seqName' is col 1, 'Nextclade_pango' is col 3
seq_col = nc_df.columns[1]
lin_col = nc_df.columns[3]
new_to_lineage = dict(zip(nc_df[seq_col], nc_df[lin_col]))

# Combine: orig name -> lineage
tip_to_lineage = {}
for orig, new in orig_to_new.items():
    if new in new_to_lineage:
        tip_to_lineage[orig] = new_to_lineage[new]
    else:
        tip_to_lineage[orig] = "Unknown"

# Lineage color mapping (using your previous palette)
lineage_colors = {
    "C.1.1": "#96c37d",
    "E.3": "#e79a90",
    "C.1": "#efbc91",
    "F.2": "#fae5b8",
    "B.1": "#c82423",
    "E.1": "#d87659",
    "C.1.3": "#80b1d3",
    "E.4": "#bc80bd",
    "B.1.3": "#fccde5",
    "Unknown": "#cccccc"
}

# Add any new lineages found in nextclade that weren't in the original config
found_lineages = set(tip_to_lineage.values())
for lin in found_lineages:
    if lin not in lineage_colors:
        lineage_colors[lin] = "#333333" # Fallback dark gray

# 2. Parse Exact Dates
dates = {}
with open(date_config, 'r') as f:
    in_data = False
    for line in f:
        line = line.strip()
        if line == "DATA":
            in_data = True
            continue
        if in_data and line and not line.startswith("#"):
            parts = line.split(',')
            if len(parts) >= 2:
                tip_name = parts[0]
                date_str = parts[1]
                if date_str != 'nan':
                    try:
                        dt = datetime.strptime(date_str, '%Y-%m-%d')
                        dates[tip_name] = mdates.date2num(dt)
                    except ValueError:
                        pass

# 3. Read Tree and Calculate Depths
tree = Phylo.read(nwk_path, "newick")

def get_depths(tree):
    depths = {tree.root: 0}
    for clade in tree.find_clades(order='preorder'):
        for child in clade:
            branch_len = child.branch_length if child.branch_length is not None else 0
            depths[child] = depths[clade] + branch_len
    return depths

depths = get_depths(tree)

anchor_depths = []
anchor_dates = []
for term in tree.get_terminals():
    if term.name in dates:
        anchor_depths.append(depths[term])
        anchor_dates.append(dates[term.name])

x_coords = {}
if len(anchor_depths) > 1:
    slope, intercept = np.polyfit(anchor_depths, anchor_dates, 1)
    for node in tree.find_clades():
        x_coords[node] = depths[node] * slope + intercept
else:
    for node in tree.find_clades():
        x_coords[node] = depths[node]

y_coords = {}
terminals = tree.get_terminals()
for i, term in enumerate(terminals):
    y_coords[term] = i

for clade in tree.get_nonterminals(order='postorder'):
    if clade.clades:
        y_coords[clade] = sum(y_coords[child] for child in clade.clades) / len(clade.clades)
    else:
        y_coords[clade] = 0

# 4. Set up plot
fig, ax = plt.subplots(figsize=(18, 9))

def draw_clade(clade):
    x = x_coords[clade]
    y = y_coords[clade]
    if clade.clades:
        y_min = min(y_coords[c] for c in clade.clades)
        y_max = max(y_coords[c] for c in clade.clades)
        ax.plot([x, x], [y_min, y_max], color='black', lw=1.2)
        
        for child in clade.clades:
            cx = x_coords[child]
            cy = y_coords[child]
            ax.plot([x, cx], [cy, cy], color='black', lw=1.2)
            draw_clade(child)

draw_clade(tree.root)

# Draw Tips with NEXTCLADE Lineage Colors
for term in terminals:
    x = x_coords[term]
    y = y_coords[term]
    name = term.name
    
    lin = tip_to_lineage.get(name, "Unknown")
    color = lineage_colors.get(lin, "#ffffff")
        
    ax.scatter(x, y, s=80, facecolor=color, edgecolor='black', linewidth=1, zorder=10)

ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.get_yaxis().set_visible(False)

locator = mdates.MonthLocator(interval=2)
formatter = mdates.DateFormatter('%b %Y')
ax.xaxis.set_major_locator(locator)
ax.xaxis.set_major_formatter(formatter)
ax.tick_params(axis='x', labelsize=12, rotation=0)

ax.xaxis.grid(True, linestyle='--', color='gray', alpha=0.5)

min_date = min(x_coords.values())
max_date = max(x_coords.values())
margin = (max_date - min_date) * 0.05
ax.set_xlim(min_date - margin, max_date + margin)

# Add Legend (Only for lineages actually present in the plot)
present_lineages = set([tip_to_lineage.get(term.name, "Unknown") for term in terminals])
# Sort them according to a predefined order, others at the end
legend_order_pref = ["C.1.1", "C.1", "C.1.3", "E.3", "E.1", "E.4", "F.2", "B.1", "B.1.3", "Unknown"]
legend_order = []
for lin in legend_order_pref:
    if lin in present_lineages:
        legend_order.append((lineage_colors[lin], lin))
for lin in present_lineages:
    if lin not in legend_order_pref:
        legend_order.append((lineage_colors[lin], lin))

handles = []
for color, label in legend_order:
    h = mlines.Line2D([], [], color='w', marker='o', markerfacecolor=color, 
                      markeredgecolor='black', markersize=10, label=label)
    handles.append(h)

ax.legend(handles=handles, title="Nextclade Lineage", loc='upper left', frameon=False, 
          fontsize=12, title_fontsize=14, borderpad=1, labelspacing=0.8)

plt.tight_layout()

pdf_path = os.path.join(output_dir, "tree_time_axis.pdf")
svg_path = os.path.join(output_dir, "tree_time_axis.svg")
png_path = os.path.join(output_dir, "tree_time_axis.png")

plt.savefig(pdf_path, bbox_inches='tight')
plt.savefig(svg_path, bbox_inches='tight')
plt.savefig(png_path, bbox_inches='tight', dpi=300)

print("Success! Tree replotted using Nextclade lineages.")
