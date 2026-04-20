#!/bin/bash
#SBATCH --job-name=postqc
#SBATCH --output=logs/postqc_%j.out
#SBATCH --error=logs/postqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

echo "=========================================="
echo "Post-Cleaning QC Batch Job"
echo "=========================================="
date
echo ""

# Load modules
module load fastqc/0.12.1
source ~/.bashrc
conda activate qc-tools

# Navigate to directory
cd ~/gbba-biohack-2026-cellsquad/challenge1-qc-detective

# Run post-QC script
bash scripts/04_post_qc.sh

echo ""
echo "=========================================="
echo "Post-QC Complete!"
echo "=========================================="
date
