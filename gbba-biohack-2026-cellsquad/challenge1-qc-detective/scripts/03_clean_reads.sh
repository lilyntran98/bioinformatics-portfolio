#!/bin/bash

#############################################
# GBBA BioHack 2026 - Challenge 1
# Script 3: Clean Reads
# Uses fastp for adapter trimming, quality filtering, and deduplication
#############################################

set -e

echo "=========================================="
echo "Cleaning RNA-seq Reads"
echo "=========================================="
echo ""

# Configuration
RAWDIR="data/raw"
CLEANDIR="data/cleaned"
REPORTDIR="results/fastp_reports"
THREADS=4

# Create output directories
mkdir -p ${CLEANDIR}
mkdir -p ${REPORTDIR}

# Sample names
SAMPLES=("NaHCO3rep1" "NaHCO3rep2" "NaClrep1" "NaClrep2")

echo "Cleaning strategy:"
echo "  1. Remove adapters"
echo "  2. Filter low quality bases (Q<20)"
echo "  3. Remove short reads (<50bp)"
echo "  4. Deduplicate (remove PCR duplicates)"
echo ""

# Clean each sample
for sample in ${SAMPLES[@]}; do
    echo "=========================================="
    echo "Cleaning: ${sample}"
    echo "=========================================="
    
    fastp \
        -i ${RAWDIR}/${sample}_R1.fastq.gz \
        -I ${RAWDIR}/${sample}_R2.fastq.gz \
        -o ${CLEANDIR}/${sample}_clean_R1.fastq.gz \
        -O ${CLEANDIR}/${sample}_clean_R2.fastq.gz \
        --detect_adapter_for_pe \
        --qualified_quality_phred 20 \
        --length_required 50 \
        --dedup \
        --dup_calc_accuracy 6 \
        --thread ${THREADS} \
        --html ${REPORTDIR}/${sample}_fastp.html \
        --json ${REPORTDIR}/${sample}_fastp.json
    
    echo "✓ ${sample} cleaned!"
    echo ""
done

echo "=========================================="
echo "Cleaning Complete!"
echo "=========================================="
echo ""

# Summary statistics
echo "Summary:"
echo "----------------------------------------"
for sample in ${SAMPLES[@]}; do
    raw_reads=$(zcat ${RAWDIR}/${sample}_R1.fastq.gz | wc -l | awk '{print $1/4}')
    clean_reads=$(zcat ${CLEANDIR}/${sample}_clean_R1.fastq.gz | wc -l | awk '{print $1/4}')
    retention=$(awk "BEGIN {printf \"%.1f\", ($clean_reads/$raw_reads)*100}")
    
    printf "%-15s Raw: %'10.0f  Clean: %'10.0f  Retained: %5.1f%%\n" \
        "${sample}:" "$raw_reads" "$clean_reads" "$retention"
done

echo ""
echo "Cleaned files: ${CLEANDIR}/"
echo "Reports: ${REPORTDIR}/"
echo ""
echo "Next: Run post-cleaning QC with scripts/04_post_qc.sh"
