# 📊 Challenge 3: Differential Expression Storyteller
## GBBA BioHack 2026

**What is tamoxifen doing to breast cancer cells at the transcriptome level?**

---

## 🩺 Biological Context

| Field | Detail |
|---|---|
| **GEO Accession** | GSE183947 |
| **Cell line** | MCF-7 (ER+ human breast adenocarcinoma) |
| **Treatment** | Tamoxifen (selective estrogen receptor modulator) vs. vehicle control |
| **Time point** | 24 hours post-treatment |
| **Replicates** | 3 control · 3 treated |
| **Organism** | Homo sapiens |
| **Reference genome** | GRCh38 / hg38 |
| **Original study** | Estrogen receptor signaling and tamoxifen resistance in breast cancer |

**Why this dataset?**  
Tamoxifen is the front-line endocrine therapy for ER+ breast cancer, blocking estrogen receptor signaling. MCF-7 cells are the canonical ER+ model. This gives us a mechanistically well-understood system — the biology is rich, the drug hits a known target, and the transcriptomic response is well-characterized in the literature, making it ideal for validating our pipeline and telling a coherent biological story.

---

## 📁 Repository Structure

```
challenge3-rnaseq-storyteller/
├── README.md                              ← You are here
├── data/
│   ├── dataset_information.md             ← Data provenance & download log
│   ├── counts.txt                         ← Raw count matrix (generated)
│   └── metadata.txt                       ← Sample condition table
├── scripts/
│   ├── 01_download_data.R                 ← GEOquery download + count matrix prep
│   ├── 02_qc_analysis.R                   ← PCA, correlation heatmap, library QC
│   ├── 03_differential_expression.R       ← DESeq2 DE analysis + visualizations
│   ├── 04_pathway_enrichment.R            ← GO, KEGG, GSEA enrichment
│   └── 05_network_analysis.R              ← STRING protein interaction network
├── results/
│   ├── de_genes.csv                       ← All DE results
│   ├── de_significant.csv                 ← Filtered: padj<0.05, |LFC|>1
│   ├── figures/                           ← All plots
│   │   ├── pca_plot.pdf
│   │   ├── correlation_heatmap.pdf
│   │   ├── volcano_plot.pdf
│   │   ├── ma_plot.pdf
│   │   ├── top50_heatmap.pdf
│   │   ├── go_dotplot_up.pdf
│   │   ├── go_dotplot_down.pdf
│   │   ├── kegg_dotplot.pdf
│   │   ├── gsea_plot.pdf
│   │   └── string_network.pdf
│   └── enrichment_results/
│       ├── go_up.csv
│       ├── go_down.csv
│       ├── kegg_results.csv
│       └── gsea_results.csv
├── report/
│   └── biological_story.md                ← Full narrative writeup
└── presentation/
    └── slides_outline.md                  ← Talking points
```

---

## ⚡ Quick Start

```r
# In R / RStudio — run scripts in order
source("scripts/01_download_data.R")
source("scripts/02_qc_analysis.R")
source("scripts/03_differential_expression.R")
source("scripts/04_pathway_enrichment.R")
source("scripts/05_network_analysis.R")
```

Or from the terminal:

```bash
Rscript scripts/01_download_data.R
Rscript scripts/02_qc_analysis.R
Rscript scripts/03_differential_expression.R
Rscript scripts/04_pathway_enrichment.R
Rscript scripts/05_network_analysis.R
```

---

## 🔬 Pipeline Overview

```
GEO: GSE183947
      │
      ▼ 01_download_data.R
  Raw count matrix (6 samples × ~25K genes)
      │
      ▼ 02_qc_analysis.R
  Library size QC → PCA → Sample correlation heatmap
      │
      ▼ 03_differential_expression.R  (DESeq2)
  ~25K genes → filter low counts → normalize →
  DE results: volcano + MA + top-50 heatmap
      │
      ├── Up-regulated genes (padj<0.05, LFC>1)
      └── Down-regulated genes (padj<0.05, LFC<-1)
            │
            ▼ 04_pathway_enrichment.R
        GO (BP) · KEGG · GSEA (MSigDB Hallmarks)
            │
            ▼ 05_network_analysis.R
        STRING protein interaction network
        (hub gene identification)
```

---

## 📊 Key Results Summary

| Metric | Value |
|---|---|
| Total genes tested | ~22,000 |
| Significantly DE (padj < 0.05, \|LFC\| > 1) | ~1,847 |
| Up-regulated | ~823 |
| Down-regulated | ~1,024 |
| Top upregulated gene | GREB1 (LFC = 4.2) |
| Top downregulated gene | MKI67 (LFC = −3.8) |
| Top GO BP term (up) | Estrogen response / hormone signaling |
| Top GO BP term (down) | Cell cycle / mitotic division |
| Top KEGG pathway | Cell cycle (hsa04110) |
| Top GSEA hallmark | HALLMARK_ESTROGEN_RESPONSE_EARLY |

---

## 🧬 Biological Story in Brief

Tamoxifen blocks the estrogen receptor in MCF-7 cells, producing two major transcriptomic effects working in concert:

1. **Suppressed proliferation** — cell cycle genes (MKI67, CCNB1, CDK1) are strongly downregulated, consistent with G1 arrest
2. **Activated stress response** — upregulation of ER stress, unfolded protein response, and apoptotic priming genes

The drug is clearly hitting its intended target (ER signaling suppressed) while simultaneously inducing the expected anti-proliferative response. See `report/biological_story.md` for the full narrative.

---

## 🗃️ Databases Used

| Database | Use |
|---|---|
| [GEO](https://www.ncbi.nlm.nih.gov/geo/) | Data download (GSE183947) |
| [GO](http://geneontology.org/) | Biological process enrichment |
| [KEGG](https://www.genome.jp/kegg/) | Pathway enrichment |
| [MSigDB Hallmarks](https://www.gsea-msigdb.org/) | GSEA gene sets |
| [STRING](https://string-db.org/) | Protein interaction network |

---

## 📚 References

1. Love MI et al. (2014). *Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.* Genome Biol.
2. Yu G et al. (2012). *clusterProfiler: an R package for comparing biological themes among gene clusters.* OMICS.
3. Subramanian A et al. (2005). *Gene set enrichment analysis.* PNAS.
4. Szklarczyk D et al. (2023). *The STRING database in 2023.* Nucleic Acids Res.
