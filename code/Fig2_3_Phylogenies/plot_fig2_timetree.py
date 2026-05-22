import argparse
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import matplotlib.dates as mdates
import pandas as pd
from Bio import Phylo
import os
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser(description="Plot Figure 2 TimeTree")
    parser.add_argument("--rename_tsv", default="../../data/197_shenzhen_rename.tsv", help="Path to renaming TSV")
    parser.add_argument("--nextclade_csv", default="../../data/nextclade.csv", help="Path to Nextclade CSV")
    parser.add_argument("--timetree", default="../../results/treetime_out/timetree.nexus", help="Path to TreeTime nexus tree")
    parser.add_argument("--dates_tsv", default="../../results/treetime_out/dates.tsv", help="Path to TreeTime dates.tsv")
    parser.add_argument("--output_dir", default="../../results/", help="Directory to save the plots")
    return parser.parse_args()

def decimal_year_to_date(year_decimal):
    year_decimal = float(year_decimal)
    year = int(year_decimal)
    rem = year_decimal - year
    base = datetime(year, 1, 1)
    result = base + pd.Timedelta(seconds=(base.replace(year=base.year + 1) - base).total_seconds() * rem)
    return result

def get_depths_years(tree):
    depths = {tree.root: 0}
    for clade in tree.find_clades(order='preorder'):
        for child in clade:
            branch_len = child.branch_length if child.branch_length is not None else 0
            depths[child] = depths[clade] + branch_len
    return depths

def main():
    args = parse_args()
    
    os.makedirs(args.output_dir, exist_ok=True)

    # 1. Parse Nextclade Lineages
    rename_df = pd.read_csv(args.rename_tsv, sep='\t', header=0)
    orig_to_new = dict(zip(rename_df.iloc[:, 0], rename_df.iloc[:, 1]))

    nc_df = pd.read_csv(args.nextclade_csv, sep=';', low_memory=False)
    seq_col = nc_df.columns[1]
    lin_col = nc_df.columns[3]
    new_to_lineage = dict(zip(nc_df[seq_col], nc_df[lin_col]))

    tip_to_lineage = {}
    for orig, new in orig_to_new.items():
        tip_to_lineage[orig] = new_to_lineage.get(new, "Unknown")

    lineage_colors = {
        "C.1.1": "#96c37d", "E.3": "#e79a90", "C.1": "#efbc91",
        "F.2": "#fae5b8", "B.1": "#c82423", "E.1": "#d87659",
        "C.1.3": "#80b1d3", "E.4": "#bc80bd", "B.1.3": "#fccde5",
        "Unknown": "#cccccc"
    }

    # 2. Read Timetree & Dates
    tree = Phylo.read(args.timetree, "nexus")
    dates_df = pd.read_csv(args.dates_tsv, sep='\t')

    node_dates_numeric = {}
    for _, row in dates_df.iterrows():
        node_name = row['#node'] if '#node' in row else row.iloc[0]
        node_dates_numeric[node_name] = row['numeric date']

    x_coords = {}
    for node in tree.find_clades():
        name = node.name
        if name in node_dates_numeric:
            val = node_dates_numeric[name]
            try:
                dt = decimal_year_to_date(val)
                x_coords[node] = mdates.date2num(dt)
            except ValueError:
                x_coords[node] = 0
        else:
            x_coords[node] = 0

    # Backfill missing dates from root
    depths_yr = get_depths_years(tree)
    root_date_yr = None
    for term in tree.get_terminals():
        if term.name in node_dates_numeric:
            try:
                root_date_yr = float(node_dates_numeric[term.name]) - depths_yr[term]
                break
            except ValueError:
                pass

    if root_date_yr is not None:
        for node in tree.find_clades():
            if x_coords.get(node, 0) == 0:
                yr = root_date_yr + depths_yr[node]
                dt = decimal_year_to_date(yr)
                x_coords[node] = mdates.date2num(dt)

    y_coords = {}
    terminals = tree.get_terminals()
    for i, term in enumerate(terminals):
        y_coords[term] = i

    for clade in tree.get_nonterminals(order='postorder'):
        if clade.clades:
            y_coords[clade] = sum(y_coords[child] for child in clade.clades) / len(clade.clades)
        else:
            y_coords[clade] = 0

    # 3. Plotting
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

    # Draw Tips
    for term in terminals:
        x = x_coords[term]
        y = y_coords[term]
        lin = tip_to_lineage.get(term.name, "Unknown")
        color = lineage_colors.get(lin, "#333333")
        ax.scatter(x, y, s=80, facecolor=color, edgecolor='black', linewidth=1, zorder=10)

    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.get_yaxis().set_visible(False)

    min_date, max_date = min(x_coords.values()), max(x_coords.values())
    margin = (max_date - min_date) * 0.05
    ax.set_xlim(min_date - margin, max_date + margin)

    ax.xaxis.set_major_locator(mdates.MonthLocator(interval=6))
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%Y-%m'))
    ax.tick_params(axis='x', labelsize=18, rotation=0)
    ax.set_ylim(0, len(terminals) + 1)
    ax.xaxis.grid(True, linestyle='--', color='gray', alpha=0.5)

    # Legend
    present_lineages = set([tip_to_lineage.get(term.name, "Unknown") for term in terminals])
    legend_order_pref = ["C.1.1", "C.1", "C.1.3", "E.3", "E.1", "E.4", "F.2", "B.1", "B.1.3", "Unknown"]
    
    handles = []
    for lin in legend_order_pref + list(present_lineages - set(legend_order_pref)):
        if lin in present_lineages:
            color = lineage_colors.get(lin, "#333333")
            handles.append(mlines.Line2D([], [], color='w', marker='o', markerfacecolor=color, 
                                         markeredgecolor='black', markersize=10, label=lin))

    ax.legend(handles=handles, title="Lineage", loc='upper left', frameon=False, 
              fontsize=18, title_fontsize=18, borderpad=1, labelspacing=0.8)

    plt.tight_layout()

    # Save
    base_out = os.path.join(args.output_dir, "fig2_timetree")
    plt.savefig(f"{base_out}.pdf", bbox_inches='tight')
    plt.savefig(f"{base_out}.svg", bbox_inches='tight')
    plt.savefig(f"{base_out}.png", bbox_inches='tight', dpi=300)
    print(f"Saved Figure 2 TimeTree to {base_out}.[pdf|svg|png]")

if __name__ == "__main__":
    main()
