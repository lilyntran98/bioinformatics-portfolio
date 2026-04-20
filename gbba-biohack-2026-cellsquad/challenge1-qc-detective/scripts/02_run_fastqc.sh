#!/bin/bash

#############################################
# GBBA BioHack 2026 - Challenge 1
# Script 2: Initial Quality Control
# Run FastQC and MultiQC on raw data
#############################################

set -e

echo "=========================================="
echo "Running Quality Control on Raw Data"
echo "=========================================="
echo ""

# Configuration
DATADIR="data/raw"
OUTDIR="results/fastqc_raw"
MULTIQC_OUT="results/multiqc_raw"
THREADS=4

# Check if data exists
if [ ! -f "${DATADIR}/NaClrep1_R1.fastq.gz" ]; then
    echo "Error: Raw data not found!"
    echo "Please run 01_download_data.sh first"
    exit 1
fi

# Create output directories
mkdir -p ${OUTDIR}
mkdir -p ${MULTIQC_OUT}

# Count files
NUM_FILES=$(ls ${DATADIR}/*.fastq.gz | wc -l)
echo "Found ${NUM_FILES} FASTQ files to analyze"
echo ""

echo "Running FastQC on raw FASTQ files..."
echo "This will take ~5-10 minutes..."
echo ""

for file in ${DATADIR}/*.fastq.gz; do
    echo "Processing $(basename $file)..."
    fastqc $file -o ${OUTDIR}
done

# Run MultiQC to aggregate results
echo "Running MultiQC to aggregate results..."
multiqc ${OUTDIR} \
    -o ${MULTIQC_OUT} \
    -n multiqc_report \
    --force

echo ""
echo "✓ MultiQC complete!"
echo ""
echo "=========================================="
echo "Quality Control Complete!"
echo "=========================================="
echo ""
echo "Results saved to:"
echo "  - Individual reports: ${OUTDIR}/"
echo "  - Aggregated report:  ${MULTIQC_OUT}/multiqc_report.html"
echo ""
echo "To view the report on your local computer:"
echo "  scp tran.nhi2@login.explorer.northeastern.edu:~/gbba-biohack-2026-cellsquad/challenge1-qc-detective/${MULTIQC_OUT}/multiqc_report.html ~/Downloads/"
echo ""
echo "Then open ~/Downloads/multiqc_report.html in your browser"
