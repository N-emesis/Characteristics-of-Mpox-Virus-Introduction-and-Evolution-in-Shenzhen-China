import pandas as pd
import matplotlib.pyplot as plt
import math
from pathlib import Path

BASE_DIR = Path("/mnt/d/hanjia/2.10/dnds")
TSV_PATH = BASE_DIR / "output_lineage/dnds_by_lineage.tsv"
OUTPUT_DIR = BASE_DIR / "output_lineage"

def make_plot():
    # Read the data
    df = pd.read_csv(TSV_PATH, sep="\t")
    
    # Filter out Genome-wide as we want gene-level
    df = df[df["gene"] != "Genome-wide"].copy()
    
    # We only care about C.1 and E.3
    df = df[df["lineage"].isin(["C.1", "E.3"])]
    
    if df.empty:
        print("No C.1 or E.3 data found!")
        return

    # Extract the ordering of genes from the original dataset
    ordered_genes = df["gene"].unique()
    gene_positions = {gene: idx for idx, gene in enumerate(ordered_genes)}
    df["x"] = df["gene"].map(gene_positions)
    
    fig, axes = plt.subplots(2, 1, figsize=(34, 12), sharex=True, constrained_layout=True)
    
    lineages = [("C.1", axes[0], "A"), ("E.3", axes[1], "B")]
    colors = {"Global": "#3b82b8", "Shenzhen": "#d96b45"}
    
    for lineage, ax, panel_label in lineages:
        subset = df[df["lineage"] == lineage].sort_values("x")
        
        for dataset in ["Global", "Shenzhen"]:
            dataset_subset = subset[subset["dataset"] == dataset]
            if dataset_subset.empty:
                continue
                
            sample_n = dataset_subset["samples"].iloc[0]
            label = f"{dataset} (n={sample_n})"
            
            ax.scatter(
                dataset_subset["x"],
                dataset_subset["dnds_adj"],
                s=40,
                color=colors[dataset],
                edgecolor="black",
                linewidth=0.3,
                alpha=0.8,
                label=label
            )
            
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
        ax.set_ylabel("dN/dS", fontsize=13)
        ax.set_title(f"{panel_label}  Lineage {lineage} Comparison", loc="left", fontsize=15, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle=":", alpha=0.35)
        ax.legend(loc="upper right", fontsize=12)

        ymax = subset["dnds_adj"].replace([math.inf, -math.inf], math.nan).max()
        ymax = 1.2 if pd.isna(ymax) else max(1.2, min(ymax * 1.15, 8.0))
        ax.set_ylim(0, ymax)

    axes[-1].set_xlabel("MPXV genes / protein regions", fontsize=13)
    
    tick_step = 5
    tick_positions = [idx for idx in range(len(ordered_genes)) if idx % tick_step == 0]
    tick_labels = [ordered_genes[idx] for idx in tick_positions]
    
    axes[-1].set_xticks(tick_positions)
    axes[-1].set_xticklabels(tick_labels, rotation=90, fontsize=7)

    fig.suptitle(
        "Overlaid gene-level dN/dS comparison of C.1 and E.3 lineages between Global and Shenzhen MPXV",
        fontsize=16,
        y=1.02,
    )

    svg_path = OUTPUT_DIR / "dnds_C1_E3_overlay.svg"
    pdf_path = OUTPUT_DIR / "dnds_C1_E3_overlay.pdf"
    
    plt.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    
    print(f"Saved overlay plots to {svg_path} and {pdf_path}")

if __name__ == "__main__":
    make_plot()
