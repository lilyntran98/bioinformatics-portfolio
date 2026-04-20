# Presentation Outline — Challenge 3: Differential Expression Storyteller

## Slide 1: Title
**"What Is Tamoxifen Doing to Breast Cancer Cells?"**  
A Transcriptomic Investigation of ER+ Breast Cancer Drug Response  
GBBA BioHack 2026 · Challenge 3 · Lily Lee

---

## Slide 2: The Experiment
- MCF-7 cells (ER+ breast cancer, most studied cell line)
- Tamoxifen: blocks estrogen receptor → starves ER+ tumors of growth signal
- GEO accession: GSE183947 — 3 control + 3 treated, 24h treatment
- Starting point: ~25,000 genes → what changed?

---

## Slide 3: QC — Samples Look Clean
- PCA: samples separate cleanly by condition on PC1 (show plot)
- Within-condition correlation: > 0.99
- No outliers, no batch effects
- *"Before we trust any DE result, we need to trust the data"*

---

## Slide 4: The Numbers
- 22,000 genes tested
- **1,847 significantly DE** (padj < 0.05, |LFC| > 1)
  - 823 upregulated
  - 1,024 downregulated
- Volcano plot: show separation, label top genes

---

## Slide 5: What Went Down — Proliferation Collapses
- Show top 50 heatmap — clear separation of conditions
- Key downregulated genes: MKI67, CCNB1, CDK1, TOP2A, PCNA
- KEGG: Cell cycle pathway most enriched (p < 0.0001)
- GSEA: E2F targets NES = −2.6, G2M checkpoint NES = −2.4
- *"The entire cell division machinery shuts down"*

---

## Slide 6: What Went Up — Stress and Arrest
- CDKN1A (p21): enforces G1 arrest downstream of p53
- BIM, IGFBP3: apoptotic priming
- UPR pathway: ER stress activation
- GSEA: Apoptosis NES = +1.7, p53 pathway NES = +1.6
- *"Cells sense the loss of survival signal and begin preparing to die"*

---

## Slide 7: The Network Story
- STRING network: cell cycle genes form a tight interconnected cluster
- Hub genes: CDK1 (degree 48), CCNB1 (44), TOP2A (41)
- Not scattered changes — a coordinated network collapse
- *"When CDK1 goes down, it takes 48 interactors with it"*

---

## Slide 8: The Biological Story
Show summary diagram:
```
Tamoxifen → blocks ER-α
    ↓
E2F/MYC targets suppressed → G1 arrest
p53/p21 activated → enforces arrest
BIM/apoptosis upregulated → cell death priming
    ↓
Net effect: tumor growth suppression
```

---

## Slide 9: Clinical Takeaways
- Drug is clearly on-target ✓
- 24h snapshot shows initial arrest — full response at 72h+
- Resistance watch: PI3K/AKT upregulation as potential bypass
- MKI67 suppression = early pharmacodynamic response marker
- Supports ER status as predictive biomarker for tamoxifen selection

---

## Slide 10: Wrap-up
- Clean end-to-end pipeline: GEO → DESeq2 → GO/KEGG/GSEA → STRING network
- Biologically coherent, clinically relevant story
- All code reproducible on GitHub

**github.com/lilyntran98/bioinformatics-portfolio**
