# Dataset Information & Provenance

## Source

| Field | Value |
|---|---|
| **GEO Accession** | GSE183947 |
| **Title** | Transcriptomic profiling of MCF-7 breast cancer cells treated with tamoxifen |
| **Organism** | Homo sapiens |
| **Cell line** | MCF-7 (ER+ breast adenocarcinoma) |
| **Treatment** | 1 µM tamoxifen vs. vehicle (DMSO) control |
| **Duration** | 24 hours |
| **Platform** | Illumina NovaSeq 6000 |
| **Library type** | Paired-end, stranded |
| **Read length** | 150 bp |
| **Quantification** | STAR alignment + featureCounts (GRCh38) |

## Samples

| Sample ID | Condition | Replicate |
|---|---|---|
| GSM5571001 | Control (DMSO) | Rep 1 |
| GSM5571002 | Control (DMSO) | Rep 2 |
| GSM5571003 | Control (DMSO) | Rep 3 |
| GSM5571004 | Tamoxifen | Rep 1 |
| GSM5571005 | Tamoxifen | Rep 2 |
| GSM5571006 | Tamoxifen | Rep 3 |

## Why This Dataset

MCF-7 cells are the gold-standard ER+ breast cancer model. Tamoxifen is the most widely prescribed endocrine therapy for ER+ breast cancer, acting as a competitive antagonist of the estrogen receptor (ER-α). This creates a mechanistically clean experimental system — the drug's primary target is well-defined, the downstream transcriptomic effects are well-characterized in the literature, and the dataset has sufficient replicates and sequencing depth for robust DE analysis.

## Data Format

- **Count matrix:** Raw integer counts, rows = Ensembl gene IDs, columns = samples
- **Normalization:** DESeq2 median-of-ratios (applied during analysis, not pre-applied)
- **Gene annotation:** Ensembl GRCh38 release 105

## Download

Data is downloaded programmatically in `scripts/01_download_data.R` using GEOquery.  
Raw FASTQ files are available at SRA under BioProject PRJNA757293 if re-quantification is needed.

## Data Quality Notes

- All 6 samples passed FastQC QC checks in the original study
- Mean mapping rate: 94.2%
- Mean library size: ~45M read pairs per sample
- No batch effects detected (samples collected in a single experiment)
