#!/bin/bash
#SBATCH --job-name=fastqc_ecoli
#SBATCH --output=logs/fastqc_%j.out
#SBATCH --error=logs/fastqc_%j.err
#SBATCH --time=01:00:00
#SBATCH --partition=short
#SBATCH --nodes=1
#SBATCH --cpus-per-task=2
#SBATCH --mem=8G

echo "=========================================="
echo "FastQC Batch Job Started"
echo "=========================================="
date
echo ""

# Load modules
module load fastqc/0.12.1
source ~/.bashrc
conda activate qc-tools

# Navigate to directory
cd ~/gbba-biohack-2026-cellsquad/challenge1-qc-detective

# Create output directories
mkdir -p results/fastqc_raw
mkdir -p results/multiqc_raw

# Run FastQC on each file sequentially
echo "Running FastQC..."
for file in data/raw/*.fastq.gz; do
    echo "Processing $(basename $file)..."
    fastqc $file -o results/fastqc_raw/ -q
    echo "  ✓ Complete"
done

echo ""
echo "Running MultiQC..."
multiqc results/fastqc_raw/ \
    -o results/multiqc_raw/ \
    -n multiqc_report \
    --force \
    -q

echo ""
echo "=========================================="
echo "QC Complete!"
echo "=========================================="
echo "Results: results/multiqc_raw/multiqc_report.html"
date
