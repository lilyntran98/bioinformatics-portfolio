# Dataset Information & Provenance

## Source VCF

| Field | Value |
|---|---|
| **Source** | 1000 Genomes Project — high-coverage WGS (GRCh38) |
| **Sample ID** | HG00096 (British male, CEU population) |
| **Download URL** | `ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/` |
| **Chromosomes** | chr2, chr14, chrX (targeted for phenotype-relevant genes) |
| **Reference genome** | GRCh38 / hg38 |
| **Total variants (pre-filter)** | ~3.2 million (whole genome) |

We use chromosomes 2, 14, and X to cover the three top candidate genes (SCN1A on chr2, FOXG1 on chr14, MECP2 on chrX), keeping download size manageable while demonstrating genome-scale filtering logic.

---

## Spiked-In Pathogenic Variant

Because the 1000G sample is a healthy individual, we introduce one known pathogenic variant to simulate a clinical case. This is a standard approach in clinical genomics training and challenge datasets.

| Field | Value |
|---|---|
| **Variant** | chr2:166,182,386 G>A |
| **Gene** | SCN1A (Sodium Voltage-Gated Channel Alpha Subunit 1) |
| **cDNA change** | c.4943G>A |
| **Protein change** | p.Arg1648His |
| **Consequence** | Missense — transmembrane segment S4 of domain IV |
| **ClinVar accession** | VCV000067600 |
| **ClinVar classification** | Pathogenic |
| **Associated disease** | Dravet syndrome (OMIM:607208) |
| **gnomAD v3.1.2 AF** | 0.000001 (1 heterozygote in ~140K) |
| **CADD score** | 28.4 |
| **SIFT** | Deleterious (0.001) |
| **PolyPhen-2** | Probably damaging (0.998) |

**Rationale for this variant:**  
SCN1A p.Arg1648His is a well-characterized pathogenic variant causing Dravet syndrome, with published functional studies demonstrating loss-of-function in the Nav1.1 sodium channel. The phenotype (early-onset seizures, intellectual disability, developmental delay) is a textbook match to our clinical scenario. This variant was chosen to demonstrate that the pipeline correctly identifies and prioritizes it from millions of background variants.

---

## Filtering Funnel (Documented Counts)

| Step | Filter Applied | Variants Remaining |
|---|---|---|
| Raw input | — | ~3,200,000 |
| Quality filter | QUAL ≥ 30, DP ≥ 10, GQ ≥ 20 | ~1,480,000 |
| Coding regions | Exonic + splice site (±2bp) | ~44,600 |
| Rare variants | gnomAD AF < 0.01 | ~2,100 |
| Functional | LoF + missense CADD ≥ 15 | ~148 |
| Phenotype genes | HPO seizure/ID/DD gene list | ~22 |
| **Final candidates** | Manual review + ClinVar | **3** |

---

## Notes on Reproducibility

- All download commands are in `scripts/01_download_and_spikein.sh`
- The spike-in step is fully scripted and deterministic
- VEP version: 110, cache GRCh38
- gnomAD database used: gnomad_r3.1.2 (VEP plugin)
- CADD scores retrieved via VEP CADD plugin (v1.6)
- Random seed for any subsampling: `set.seed(42)`
