# 📥 Data Download Instructions

## Dataset Information

**Study:** The antibiotics potentiator bicarbonate causes upregulation of tryptophanase and iron acquisition proteins in *Escherichia coli*

**BioProject:** PRJNA433855  
**SRA Study:** SRP132716  
**GEO Series:** GSE111324  
**Platform:** Illumina HiSeq 2000  
**Release Date:** November 2, 2018

---

## Selected Samples

We downloaded 4 RNA-seq samples representing two biological conditions with technical replicates:

| Sample ID | GEO Accession | SRA Experiment | SRA Run | Condition | Read Pairs | Size |
|-----------|---------------|----------------|---------|-----------|------------|------|
| NaHCO3rep1 | GSM2992380 | SRX3683647 | SRR6710205 | Bicarbonate treatment | 5,520,708 | 561 MB |
| NaHCO3rep2 | GSM2992381 | SRX3683648 | SRR6710206 | Bicarbonate treatment | 5,184,530 | 535 MB |
| NaClrep1 | GSM2992382 | SRX3683649 | SRR6710207 | NaCl control | 5,998,161 | 610 MB |
| NaClrep2 | GSM2992383 | SRX3683650 | SRR6710208 | NaCl control | 5,748,488 | 585 MB |

**Total Dataset:** 22,451,887 read pairs (~2.3 GB compressed)

---

## Dataset Selection Rationale

This dataset was selected for the Quality Detective challenge because:

✅ **Older Platform (HiSeq 2000, 2018)** - Known for quality issues, ideal for diagnostic practice  
✅ **Manageable Size** - 22.5M reads suitable for hackathon timeframe  
✅ **Biological Replicates** - Enables sample comparison and quality assessment  
✅ **Public Data** - Fully reproducible from NCBI/ENA repositories  
✅ **Expected Issues** - Anticipated duplication, adapters, quality challenges

---

## Download Method: ENA (European Nucleotide Archive)

We used ENA FTP download as it proved more reliable than NCBI SRA toolkit on the Discovery cluster.

### Required Software

None! Just `wget` which is pre-installed on Discovery.

---

## Download Commands

### Automated Script (Recommended)
```bash
# Run the download script
bash scripts/01_download_data_ena.sh
```

This will:
1. Download all 4 samples from ENA via FTP
2. Rename files with descriptive sample names
3. Verify read counts match expected values

---

### Manual Download

If you need to download samples individually:
```bash
# Navigate to data directory
cd data/raw

# Download NaHCO3rep1 (SRR6710205)
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/005/SRR6710205/SRR6710205_1.fastq.gz \
    -O NaHCO3rep1_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/005/SRR6710205/SRR6710205_2.fastq.gz \
    -O NaHCO3rep1_R2.fastq.gz

# Download NaHCO3rep2 (SRR6710206)
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/006/SRR6710206/SRR6710206_1.fastq.gz \
    -O NaHCO3rep2_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/006/SRR6710206/SRR6710206_2.fastq.gz \
    -O NaHCO3rep2_R2.fastq.gz

# Download NaClrep1 (SRR6710207)
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/007/SRR6710207/SRR6710207_1.fastq.gz \
    -O NaClrep1_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/007/SRR6710207/SRR6710207_2.fastq.gz \
    -O NaClrep1_R2.fastq.gz

# Download NaClrep2 (SRR6710208)
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/008/SRR6710208/SRR6710208_1.fastq.gz \
    -O NaClrep2_R1.fastq.gz
wget -q --show-progress \
    ftp://ftp.sra.ebi.ac.uk/vol1/fastq/SRR671/008/SRR6710208/SRR6710208_2.fastq.gz \
    -O NaClrep2_R2.fastq.gz
```

---

## Expected Files

After download, you should have:
```
data/raw/
├── NaHCO3rep1_R1.fastq.gz  (5,520,708 reads, ~280 MB)
├── NaHCO3rep1_R2.fastq.gz  (5,520,708 reads, ~281 MB)
├── NaHCO3rep2_R1.fastq.gz  (5,184,530 reads, ~267 MB)
├── NaHCO3rep2_R2.fastq.gz  (5,184,530 reads, ~268 MB)
├── NaClrep1_R1.fastq.gz    (5,998,161 reads, ~305 MB)
├── NaClrep1_R2.fastq.gz    (5,998,161 reads, ~305 MB)
├── NaClrep2_R1.fastq.gz    (5,748,488 reads, ~292 MB)
└── NaClrep2_R2.fastq.gz    (5,748,488 reads, ~293 MB)
```

**Total Size:** ~2.3 GB compressed

---

## Verification
```bash
# Check file sizes
ls -lh data/raw/*.fastq.gz

# Verify read counts
for file in data/raw/*_R1.fastq.gz; do
    reads=$(zcat $file | wc -l | awk '{print $1/4}')
    printf "%-30s %'d reads\n" "$(basename $file):" "$reads"
done
```

**Expected Output:**
```
NaHCO3rep1_R1.fastq.gz:     5,520,708 reads
NaHCO3rep2_R1.fastq.gz:     5,184,530 reads
NaClrep1_R1.fastq.gz:       5,998,161 reads
NaClrep2_R1.fastq.gz:       5,748,488 reads
```

---

## Troubleshooting

### Download Interrupted?
Simply re-run the wget command - it will resume where it left off.

### Need to Download Elsewhere?

**Alternative Source - NCBI SRA:**
```bash
# Install SRA toolkit
conda install -c bioconda sra-tools

# Download using SRA toolkit
fasterq-dump SRR6710205 --split-files --progress
fasterq-dump SRR6710206 --split-files --progress
fasterq-dump SRR6710207 --split-files --progress
fasterq-dump SRR6710208 --split-files --progress

# Compress and rename as needed
```

---

## Reference Genome (Optional)

If needed for downstream analysis:
```bash
cd data/reference

# Download E. coli K-12 MG1655 reference genome
wget https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/005/845/GCF_000005845.2_ASM584v2/GCF_000005845.2_ASM584v2_genomic.fna.gz

gunzip GCF_000005845.2_ASM584v2_genomic.fna.gz
mv GCF_000005845.2_ASM584v2_genomic.fna ecoli_K12_MG1655.fna
```

---

## Library Construction Protocol

From SRA metadata:  
*"At the end of the incubation period, cultures were transferred to ice and bacterial cells were collected at 4°C by centrifugation at 10,000 × g for a 10 min period. The collected cells were then subjected to RNA extraction."*

**Technical Details:**
- **Library Type:** Paired-end cDNA
- **Read Length:** 150 bp
- **Library Selection:** cDNA
- **Library Source:** Transcriptomic

---

## Quality Issues Found

After downloading, initial QC revealed:

⚠️ **Critical Issues:**
- Very high duplication: 68-71% (vs expected 20-40%)
- Root cause: Over-amplification during library prep

⚠️ **Moderate Issues:**
- Adapter contamination: ~4%

✅ **Positive Findings:**
- Excellent base quality: Q37 average
- Correct GC content: 48-50%

**See full diagnosis:** [`report/diagnosis_report.md`](../report/diagnosis_report.md)

---

## Data Availability Statement

All sequencing data used in this analysis are publicly available through:

- **NCBI SRA:** https://www.ncbi.nlm.nih.gov/bioproject/PRJNA433855
- **ENA:** https://www.ebi.ac.uk/ena/browser/view/PRJNA433855
- **GEO:** https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE111324

Individual sample accessions are listed in the table above.

---

## Time Estimate

- **Download time:** ~20-30 minutes (depends on network speed)
- **Storage required:** ~2.3 GB compressed, ~10 GB uncompressed

---

## Next Steps

After downloading:

1. **Run initial QC:** `bash scripts/02_run_fastqc.sh`
2. **View QC report:** `results/multiqc_raw/multiqc_report.html`
3. **Diagnose issues:** Analyze MultiQC report
4. **Clean data:** `bash scripts/03_clean_reads.sh`
5. **Verify cleaning:** `bash scripts/04_post_qc.sh`

---

**Data downloaded successfully on:** January 22-23, 2026  
**Platform:** NEU Discovery Cluster  
**Team:** CellSquad
