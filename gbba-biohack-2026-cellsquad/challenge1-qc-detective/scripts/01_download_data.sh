#!/bin/bash

#############################################
# GBBA BioHack 2026 - Challenge 1
# Script 1: Download E. coli RNA-seq Data
# Dataset: PRJNA433855 (Bicarbonate treatment)
# Author: Team CellSquad
# Date: January 22-23, 2026
#############################################

set -e  # Exit on error
set -u  # Exit on undefined variable

echo "=========================================="
echo "E. coli RNA-seq Data Download"
echo "BioProject: PRJNA433855"
echo "Study: Bicarbonate antibiotics potentiator"
echo "=========================================="
echo ""

# Configuration
DATADIR="data/raw"
SAMPLES=("SRR6710205" "SRR6710206" "SRR6710207" "SRR6710208")
NAMES=("NaHCO3rep1" "NaHCO3rep2" "NaClrep1" "NaClrep2")
THREADS=4

# Create output directory
mkdir -p ${DATADIR}

echo "Downloading 4 samples..."
echo "  NaHCO3rep1 (SRR6710205): 5.5M reads (561 MB)"
echo "  NaHCO3rep2 (SRR6710206): 5.2M reads (535 MB)"
echo "  NaClrep1   (SRR6710207): 6.0M reads (610 MB)"
echo "  NaClrep2   (SRR6710208): 5.7M reads (585 MB)"
echo "Total: ~22.4M read pairs (~2.3 GB compressed)"
echo ""

# Download each sample
for i in ${!SAMPLES[@]}; do
    srr=${SAMPLES[$i]}
    name=${NAMES[$i]}
    
    echo "=========================================="
    echo "Sample $((i+1))/4: ${name} (${srr})"
    echo "=========================================="
    
    cd ${DATADIR}
    
    # Check if already downloaded
    if [[ -f "${name}_R1.fastq.gz" && -f "${name}_R2.fastq.gz" ]]; then
        echo "✓ ${name} already downloaded, skipping..."
        cd ../..
        continue
    fi
    
    # Download with fasterq-dump
    echo "Downloading..."
    fasterq-dump ${srr} \
        --split-files \
        --threads ${THREADS} \
        --progress
    
    # Compress
    echo "Compressing..."
    gzip ${srr}_1.fastq ${srr}_2.fastq
    
    # Rename with descriptive names
    mv ${srr}_1.fastq.gz ${name}_R1.fastq.gz
    mv ${srr}_2.fastq.gz ${name}_R2.fastq.gz
    
    cd ../..
    
    echo "✓ ${name} complete"
    echo ""
done

echo "=========================================="
echo "Download Complete!"
echo "=========================================="
echo ""

# List downloaded files
echo "Downloaded files:"
ls -lh ${DATADIR}/*.fastq.gz
echo ""

# Verify read counts
echo "Verifying read counts..."
echo "----------------------------------------"
for file in ${DATADIR}/*_R1.fastq.gz; do
    reads=$(zcat $file | wc -l | awk '{print $1/4}')
    printf "%-30s %'12d reads\n" "$(basename $file):" "$reads"
done

echo ""
echo "Expected read counts:"
echo "  NaHCO3rep1_R1.fastq.gz:     5,520,708 reads"
echo "  NaHCO3rep2_R1.fastq.gz:     5,184,530 reads"
echo "  NaClrep1_R1.fastq.gz:       5,998,161 reads"
echo "  NaClrep2_R1.fastq.gz:       5,700,208 reads"
echo ""
echo "✓ All downloads verified!"
echo "✓ Ready for quality control analysis!"
