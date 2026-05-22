import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# Load data
df = pd.read_csv('/mnt/d/hanjia/2.10/dnds/output_hamming_dnds/hamming_dnds_results.tsv', sep='\t')

# Filter valid values
df = df.replace([np.inf, -np.inf], np.nan).dropna(subset=['dnds_adj'])
df = df[df['gene'] != 'Genome-wide']

# Identify lineages and groups
df['lineage'] = df['group'].apply(lambda x: 'C.1.1' if 'C.1.1' in x else 'E.3')
# Fix sample_type matching based on actual group names (e.g. "C.1.1_sampled" and "C.1.1_unsampled")
df['sample_type'] = df['group'].apply(lambda x: 'Unsampled' if 'unsampled' in x else 'Sampled (100)')

# Drop any NaN or empty strings from hue column
df = df.dropna(subset=['sample_type'])
df = df[df['sample_type'] != '']

print(f"Data types in sample_type: {df['sample_type'].unique()}")

# Get ordered list of genes to keep x-axis consistent
# We'll just sort them alphabetically for now or use the order of appearance
ordered_genes = sorted(df['gene'].unique())

# Set up the plot style
sns.set_theme(style="ticks")
plt.rcParams['font.family'] = 'sans-serif'
plt.rcParams['font.sans-serif'] = ['Arial']

fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(16, 12), sharex=True)

# Colors matching the requested style (light blue/dark blue, light red/dark red)
palette_c11 = {'Sampled (100)': '#92c5de', 'Unsampled': '#2166ac'}
palette_e3 = {'Sampled (100)': '#f4a582', 'Unsampled': '#b2182b'}

# --- Plot A: C.1.1 ---
c11_df = df[df['lineage'] == 'C.1.1']

# We need to plot points for each gene. We'll use a categorical scatterplot with jitter=False
sns.pointplot(
    data=c11_df, 
    x='gene', 
    y='dnds_adj', 
    hue='sample_type',
    dodge=0.4,       # separation between the two groups
    join=False,      # don't connect points
    palette=palette_c11,
    markers='o',
    scale=1.5,
    errwidth=1.5,
    capsize=0.1,
    ax=ax1,
    order=ordered_genes
)

# Add a horizontal line at y=1 (neutral selection)
ax1.axhline(1.0, color='black', linestyle='--', linewidth=1, alpha=0.7)

# Customize Ax1
ax1.set_title('A\n\nC.1.1', loc='left', fontweight='bold', fontsize=16)
ax1.set_ylabel('dN/dS ratio', fontsize=14)
ax1.set_xlabel('')
ax1.legend(title='', loc='upper center', bbox_to_anchor=(0.5, -0.15), ncol=2, frameon=False, fontsize=12)
ax1.grid(axis='y', linestyle=':', alpha=0.6)
ax1.set_ylim(bottom=-0.5)

# --- Plot B: E.3 ---
e3_df = df[df['lineage'] == 'E.3']

sns.pointplot(
    data=e3_df, 
    x='gene', 
    y='dnds_adj', 
    hue='sample_type',
    dodge=0.4,
    join=False,
    palette=palette_e3,
    markers='o',
    scale=1.5,
    errwidth=1.5,
    capsize=0.1,
    ax=ax2,
    order=ordered_genes
)

ax2.axhline(1.0, color='black', linestyle='--', linewidth=1, alpha=0.7)

# Customize Ax2
ax2.set_title('B\n\nE.3', loc='left', fontweight='bold', fontsize=16)
ax2.set_ylabel('dN/dS ratio', fontsize=14)
ax2.set_xlabel('Genes', fontsize=14)
ax2.legend(title='', loc='upper center', bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=12)
ax2.grid(axis='y', linestyle=':', alpha=0.6)
ax2.set_ylim(bottom=-0.5)

# Rotate x-axis labels
plt.xticks(rotation=90, fontsize=10)

# Add background shading for visual grouping (optional, simulating the grey background in the reference)
# For example, shading every other gene
for ax in [ax1, ax2]:
    for i in range(1, len(ordered_genes), 2):
        ax.axvspan(i - 0.5, i + 0.5, color='lightgray', alpha=0.2, zorder=0)

plt.tight_layout()
plt.subplots_adjust(hspace=0.4)

output_pdf = '/mnt/d/hanjia/2.10/dnds/output_hamming_dnds/gene_dnds_comparison.pdf'
output_png = '/mnt/d/hanjia/2.10/dnds/output_hamming_dnds/gene_dnds_comparison.png'

plt.savefig(output_pdf, format='pdf', bbox_inches='tight')
plt.savefig(output_png, format='png', bbox_inches='tight', dpi=300)
print(f"Plots saved to {output_pdf} and {output_png}")
