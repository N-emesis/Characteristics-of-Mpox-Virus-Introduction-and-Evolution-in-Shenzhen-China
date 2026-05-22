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
    
    fig, axes = plt.subplots(4, 1, figsize=(34, 24), sharex=True, constrained_layout=True)
    
    # We want 4 panels: Global C.1, Shenzhen C.1, Global E.3, Shenzhen E.3
    panel_specs = [
        ("Global", "C.1", "#3b82b8", axes[0], "A"),
        ("Shenzhen", "C.1", "#d96b45", axes[1], "B"),
        ("Global", "E.3", "#4caf50", axes[2], "C"),
        ("Shenzhen", "E.3", "#9c27b0", axes[3], "D"),
    ]
    
    for dataset, lineage, color, ax, panel_label in panel_specs:
        subset = df[(df["dataset"] == dataset) & (df["lineage"] == lineage)].sort_values("x")
        
        if subset.empty:
            ax.text(0.5, 0.5, f"No data for {dataset} {lineage}", ha='center', va='center', transform=ax.transAxes, fontsize=20)
            ax.set_title(f"{panel_label}  {dataset} - Lineage {lineage}", loc="left", fontsize=15, fontweight="bold")
            continue
            
        sample_n = subset["samples"].iloc[0]
        
        ax.scatter(
            subset["x"],
            subset["dnds_adj"],
            s=38,
            color=color,
            edgecolor="black",
            linewidth=0.3,
            alpha=0.9,
        )
        ax.axhline(1.0, color="gray", linestyle="--", linewidth=1)
        ax.set_ylabel("dN/dS", fontsize=13)
        ax.set_title(f"{panel_label}  {dataset} - Lineage {lineage} (n={sample_n})", loc="left", fontsize=15, fontweight="bold")
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.grid(axis="y", linestyle=":", alpha=0.35)

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
        "Gene-level dN/dS comparison of C.1 and E.3 lineages between Global and Shenzhen MPXV",
        fontsize=16,
        y=1.02,
    )

    svg_path = OUTPUT_DIR / "dnds_C1_E3_comparison.svg"
    pdf_path = OUTPUT_DIR / "dnds_C1_E3_comparison.pdf"
    
    plt.savefig(svg_path, format="svg", bbox_inches="tight")
    plt.savefig(pdf_path, format="pdf", bbox_inches="tight")
    plt.close(fig)
    
    print(f"Saved plots to {svg_path} and {pdf_path}")

if __name__ == "__main__":
    make_plot()
