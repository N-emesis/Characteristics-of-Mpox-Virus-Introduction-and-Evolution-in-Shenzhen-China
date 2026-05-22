#!/bin/bash
set -e
export PATH="/root/miniforge3/envs/bio/bin:$PATH"

echo "1. Inferring tree... (Skipped, already done)"
# python3 code/infer_tree.py

echo "2. Inferring clock and mugration... (Skipped, already done)"
# python3 code/infer_clock_mugration.py

echo "3. Estimating importations... (Skipped, already done)"
# mkdir -p results/aligned.ufboot_tres
# 
# # Changed back to sequential loop because parallel execution (xargs -P 10) caused Out of Memory (OOM) error (signal 9)
# for i in $(seq 1 10); do 
#     echo "   Replicate $i"
#     python3 code/estimate_importations.py --bootstrap_replicate $i
# done

echo "4. Rarefaction..."
rm -rf results/rarefaction
python3 code/rarefaction.py --n_samples 197

echo "5. Plotting Shenzhen Rarefaction..."
python3 code/plot_rarefaction_shenzhen_only.py

echo "Done!"
