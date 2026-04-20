#!/bin/bash

set -e

echo "=========================================="
echo "Downloading from ENA"
echo "Fast alternative to NCBI SRA"
echo "=========================================="

DATADIR="data/raw"
mkdir -p ${DATADIR}
cd ${DATADIR}

# Download NaHCO3rep1 (SRR6710205)
echo "Downloading NaHCO3rep1..."
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/005/SRR6710205/SRR6710205_1.fastq.gz \
    -O NaHCO3rep1_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/005/SRR6710205/SRR6710205_2.fastq.gz \
    -O NaHCO3rep1_R2.fastq.gz
echo "✓ NaHCO3rep1 complete"
echo ""

# Download NaHCO3rep2 (SRR6710206)
echo "Downloading NaHCO3rep2..."
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/006/SRR6710206/SRR6710206_1.fastq.gz \
    -O NaHCO3rep2_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/006/SRR6710206/SRR6710206_2.fastq.gz \
    -O NaHCO3rep2_R2.fastq.gz
echo "✓ NaHCO3rep2 complete"
echo ""

# Download NaClrep1 (SRR6710207)
echo "Downloading NaClrep1..."
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/007/SRR6710207/SRR6710207_1.fastq.gz \
    -O NaClrep1_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/007/SRR6710207/SRR6710207_2.fastq.gz \
    -O NaClrep1_R2.fastq.gz
echo "✓ NaClrep1 complete"
echo ""

# Download NaClrep2 (SRR6710208)
echo "Downloading NaClrep2..."
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/008/SRR6710208/SRR6710208_1.fastq.gz \
    -O NaClrep2_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/008/SRR6710208/SRR6710208_2.fastq.gz \
    -O NaClrep2_R2.fastq.gz
echo "✓ NaClrep2 complete"
echo ""

cd ../..

echo "=========================================="
echo "Download Complete!"
echo "=========================================="
echo ""
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
echo "Expected:"
echo "  NaHCO3rep1: 5,520,708 reads"
echo "  NaHCO3rep2: 5,184,530 reads"
echo "  NaClrep1:   5,998,161 reads"
echo "  NaClrep2:   5,700,208 reads"
echo ""
echo "✓ Ready for QC analysis!"
