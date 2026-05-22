import argparse
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import os

def parse_args():
    parser = argparse.ArgumentParser(description="Plot Figure 2 Bar Charts")
    parser.add_argument("--metadata_csv", default="../../data/metadata.csv", help="Path to local metadata CSV")
    parser.add_argument("--rename_tsv", default="../../data/197_shenzhen_rename.tsv", help="Path to renaming TSV")
    parser.add_argument("--nextclade_csv", default="../../data/nextclade.csv", help="Path to Nextclade CSV")
    parser.add_argument("--global_tsv", default="../../data/global.tsv", help="Path to global metadata TSV")
    parser.add_argument("--output_svg", default="../../results/fig2_bar_chart.svg", help="Output SVG path")
    return parser.parse_args()

def main():
    args = parse_args()
    
    os.makedirs(os.path.dirname(args.output_svg), exist_ok=True)

    # 1. Load Local Data
    metadata_df = pd.read_csv(args.metadata_csv)
    metadata_df = metadata_df[~metadata_df['date'].str.contains('XX', na=False)]
    metadata_df['date'] = pd.to_datetime(metadata_df['date'])
    
    rename_df = pd.read_csv(args.rename_tsv, sep='\t', header=0)
    orig_to_new = dict(zip(rename_df.iloc[:, 0], rename_df.iloc[:, 1]))
    
    nc_df = pd.read_csv(args.nextclade_csv, sep=';', low_memory=False)
    new_to_lineage = dict(zip(nc_df.iloc[:, 1], nc_df.iloc[:, 3]))
    
    metadata_df['new_name'] = metadata_df['name'].map(orig_to_new)
    metadata_df['lineage'] = metadata_df['new_name'].map(new_to_lineage).fillna("Unknown")
    
    weekly_lineages = metadata_df.groupby([pd.Grouper(key='date', freq='W-MON'), 'lineage']).size().unstack(fill_value=0)

    # 2. Load Global Data
    global_df = pd.read_csv(args.global_tsv, sep="\t")
    if "date" not in global_df.columns or "strain" not in global_df.columns:
        raise ValueError("Missing 'date' or 'strain' column in global metadata file")
    
    global_df["date"] = pd.to_datetime(global_df["date"], errors='coerce')
    global_df = global_df.dropna(subset=["date"])
    weekly_metadata = global_df.resample("W-Mon", on="date")["strain"].size()

    # 3. Plotting
    plt.rcParams['font.family'] = 'Arial'
    plt.rcParams['font.size'] = 22

    fig, ax1 = plt.subplots(figsize=(12, 6))

    color_metadata = "#14517c"
    ax1.bar(
        weekly_metadata.index, 
        weekly_metadata, 
        color=color_metadata, 
        alpha=0.6, 
        label="Global", 
        width=6
    )

    lineage_colors = {
        "C.1.1": "#96c37d", "E.3": "#e79a90", "C.1": "#efbc91",
        "F.2": "#fae5b8", "B.1": "#c82423", "E.1": "#d87659",
        "C.1.3": "#80b1d3", "E.4": "#bc80bd", "B.1.3": "#fccde5",
        "Unknown": "#cccccc"
    }

    columns = list(weekly_lineages.columns)
    bottom = np.zeros(len(weekly_lineages))
    
    for lineage in columns:
        color = lineage_colors.get(lineage, '#999999')
        ax1.bar(
            weekly_lineages.index,
            weekly_lineages[lineage],
            bottom=bottom,
            color=color,
            alpha=1,
            label=lineage,
            width=6
        )
        bottom += weekly_lineages[lineage]
    
    ax1.tick_params(axis="y", labelcolor="#000000", labelsize=22)
    ax1.set_ylim(0, 36)
    ax1.yaxis.set_major_locator(plt.MultipleLocator(5))

    x_min, x_max = global_df["date"].min(), global_df["date"].max()
    ax1.set_xlim(x_min, x_max)

    xticks = pd.date_range(start=x_min, end=x_max, freq="2MS")
    ax1.set_xticks(xticks)
    ax1.set_xticklabels([date.strftime("%Y-%m") for date in xticks], rotation=30, ha="right")
    ax1.tick_params(axis="x", labelsize=22)

    fig.legend(loc="upper left", bbox_to_anchor=(0.6, 1), frameon=False, ncol=2)

    ax1.spines["top"].set_visible(False)
    ax1.spines["right"].set_visible(False)
    ax1.grid(False)

    plt.tight_layout()

    plt.savefig(args.output_svg, format="svg", bbox_inches='tight')
    plt.savefig(args.output_svg.replace(".svg", ".pdf"), bbox_inches='tight')
    plt.savefig(args.output_svg.replace(".svg", ".png"), bbox_inches='tight', dpi=300)

    print(f"Saved Figure 2 bar charts to {args.output_svg}")

if __name__ == "__main__":
    main()
