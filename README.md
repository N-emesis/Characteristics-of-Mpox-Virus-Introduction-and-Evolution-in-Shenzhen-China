# Mpox Virus (hMpxV) Spatiotemporal Dynamics Analysis

This repository contains the analytical pipeline, scripts, Jupyter Notebooks, and intermediate data used for the study of Mpox virus spatiotemporal dynamics.

## Repository Structure

*   **code/**: Contains all scripts used for data processing and visualization.
    *   **Fig2_3_Phylogenies/**: Scripts for phylogenetic analysis and tree visualization.
    *   **Fig6_Introduction_Events/**: Pipeline for estimating viral importation events (ML tree inference, molecular clock, ancestral state reconstruction, and rarefaction analysis).
    *   **Fig7_Epidemiology/**: Scripts for calculating $R_0$ and plotting $R_e$.
    *   **dN_dS_and_Mutations/**: Visualization scripts for dN/dS ratios and mutation profiles.
    *   **Misc_Analysis/**: Additional analysis scripts for demographic data.
    *   **Jupyter_Notebooks/**: Interactive notebooks used to generate the final figures, including a `data/` subfolder with necessary intermediate datasets.
*   **results/**: Contains inferred ML trees, TimeTrees, and estimation results.
*   **xml_files/**: Contains the phyloXML files generated during the TreeTime ancestral state reconstruction. These files include time-calibrated trees with ancestral geographic states and serve as the source data for the importation events inferred in Figure 6.

## Dependencies

*   Software: IQ-TREE, TreeTime, R
*   Python 3.8+ Packages: `pandas`, `numpy`, `scipy`, `matplotlib`, `seaborn`, `biopython`, `jupyter`
