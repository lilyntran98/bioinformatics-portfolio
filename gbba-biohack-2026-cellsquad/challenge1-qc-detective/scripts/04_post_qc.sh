#!/bin/bash

#############################################
# GBBA BioHack 2026 - Challenge 1
# Script 4: Post-Cleaning Quality Control
#############################################

set -e

echo "=========================================="
echo "Running QC on Cleaned Data"
echo "=========================================="
echo ""

# Configuration
DATADIR="data/cleaned"
OUTDIR="results/fastqc_cleaned"
MULTIQC_OUT="results/multiqc_cleaned"
THREADS=2

# Check if cleaned data exists
if [ ! -f "${DATADIR}/NaClrep1_clean_R1.fastq.gz" ]; then
    echo "Error: Cleaned data not found!"
    echo "Please run scripts/03_clean_reads.sh first"
    exit 1
fi

# Create output directories
mkdir -p ${OUTDIR}
mkdir -p ${MULTIQC_OUT}

# Count files
NUM_FILES=$(ls ${DATADIR}/*.fastq.gz | wc -l)
echo "Found ${NUM_FILES} cleaned FASTQ files to analyze"
echo ""

# Run FastQC on cleaned files
echo "Running FastQC on cleaned FASTQ files..."
echo "This will take ~5-10 minutes..."
echo ""

for file in ${DATADIR}/*.fastq.gz; do
    echo "Processing $(basename $file)..."
    fastqc $file -o ${OUTDIR} -t 1 -q
done

echo ""
echo "✓ FastQC complete!"
echo ""

# Run MultiQC
echo "Running MultiQC to aggregate results..."
multiqc ${OUTDIR} \
    -o ${MULTIQC_OUT} \
    -n multiqc_report \
    --force \
    -q

echo ""
echo "✓ MultiQC complete!"
echo ""
echo "=========================================="
echo "Post-Cleaning QC Complete!"
echo "=========================================="
echo ""
echo "Results saved to:"
echo "  - Individual reports: ${OUTDIR}/"
echo "  - Aggregated report:  ${MULTIQC_OUT}/multiqc_report.html"
echo ""
echo "Next: Compare before/after"
echo "  bash scripts/05_compare_qc.sh"
