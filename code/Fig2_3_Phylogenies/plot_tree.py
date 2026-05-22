import matplotlib.pyplot as plt
from Bio import Phylo
import os

nwk_path = "/mnt/d/hanjia/2.10/tree_time_new/tree.nwk"
output_dir = os.path.dirname(nwk_path)

try:
    print(f"Reading tree file: {nwk_path}")
    tree = Phylo.read(nwk_path, "newick")
    
    fig = plt.figure(figsize=(18, 9))
    ax = fig.add_subplot(1, 1, 1)
    
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    ax.spines['left'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.get_yaxis().set_visible(False)
    ax.get_xaxis().set_visible(False)
    
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda x: "")
    
    pdf_path = os.path.join(output_dir, "tree_no_labels.pdf")
    svg_path = os.path.join(output_dir, "tree_no_labels.svg")
    png_path = os.path.join(output_dir, "tree_no_labels.png")
    
    plt.savefig(pdf_path, bbox_inches='tight')
    plt.savefig(svg_path, bbox_inches='tight')
    plt.savefig(png_path, bbox_inches='tight', dpi=300)
    
    print(f"Tree plot (no labels, 18x9) successfully saved to:")
    print(f"- {pdf_path}")
    print(f"- {svg_path}")
    print(f"- {png_path}")
    
except Exception as e:
    print(f"Error occurred during plotting: {e}")
