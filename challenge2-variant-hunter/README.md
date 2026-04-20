# 🧬 GBBA BioHack 2026 — Challenge 2: Variant Hunter

**Team:** Lily Lee  
**Challenge:** Find the Pathogenic Mutation  
**Event:** GBBA BioHack 2026 · January 22–23, 2026  
**Points possible:** 100 + 10 bonus  

---

## 🩺 Clinical Scenario

| Field | Detail |
|---|---|
| Patient | 8-year-old child |
| Phenotypes | Intellectual disability (moderate), Seizures (onset age 2), Developmental delay, Facial dysmorphisms |
| Family history | Parents unaffected, no siblings |
| Prior testing | Chromosomal microarray — normal |
| Working hypothesis | De novo pathogenic variant in a neurodevelopmental gene |

**HPO terms used:** `HP:0001249` (Intellectual disability), `HP:0001250` (Seizures), `HP:0001263` (Developmental delay), `HP:0001畸` (Facial dysmorphism)

---

## 📁 Repository Structure

```
challenge2-variant-hunter/
├── README.md                          ← You are here
├── data/
│   └── dataset_information.md         ← Data provenance & spike-in log
├── scripts/
│   ├── 00_setup_environment.sh        ← Conda env + tool install
│   ├── 01_download_and_spikein.sh     ← Fetch 1000G VCF, spike pathogenic variant
│   ├── 02_annotate_variants.sh        ← VEP annotation pipeline
│   ├── 03_filter_variants.py          ← Python: QC → frequency → functional filtering
│   ├── 04_prioritize_candidates.R     ← R: scoring, ranking, ACMG rule application
│   └── 05_generate_report_table.R     ← R: candidate summary table for report
├── results/
│   ├── filtered_variants.vcf          ← Final filtered VCF (generated)
│   ├── top_candidates.xlsx            ← Ranked candidate table (generated)
│   └── igv_screenshots/               ← IGV images of candidate variants
├── shiny_app/
│   ├── app.R                          ← Interactive variant browser (Shiny)
│   └── README.md                      ← How to run the Shiny app
├── report/
│   └── clinical_report.md             ← Clinical interpretation writeup
└── presentation/
    └── slides_outline.md              ← Presentation talking points
```

---

## ⚡ Quick Start

```bash
# 1. Set up environment
bash scripts/00_setup_environment.sh

# 2. Download data and spike in pathogenic variant
bash scripts/01_download_and_spikein.sh

# 3. Annotate with VEP
bash scripts/02_annotate_variants.sh

# 4. Filter variants (Python)
python scripts/03_filter_variants.py \
  --input results/vep_annotated.vcf \
  --output results/filtered_variants.vcf \
  --max_af 0.01 \
  --min_cadd 15

# 5. Prioritize and score candidates (R)
Rscript scripts/04_prioritize_candidates.R

# 6. Launch interactive Shiny browser
cd shiny_app && Rscript -e "shiny::runApp('.')"
```

---

## 🔬 Pipeline Overview

```
Raw VCF (~3M variants)
        │
        ▼ 01_download_and_spikein.sh
  Patient VCF (spike-in added)
        │
        ▼ 02_annotate_variants.sh  (Ensembl VEP)
  Annotated VCF — gene, consequence, gnomAD AF, ClinVar, CADD, SIFT, PolyPhen
        │
        ▼ 03_filter_variants.py
  Step 1: Quality filter    QUAL≥30, DP≥10       → ~1.5M
  Step 2: Coding only       exonic + splice       → ~45K
  Step 3: Rare variants     gnomAD AF < 0.01      → ~2K
  Step 4: Functional        LoF + high CADD       → ~150
  Step 5: Phenotype genes   HPO seizure/ID genes  → ~20
        │
        ▼ 04_prioritize_candidates.R
  Scored & ranked candidates (ACMG rules applied)
        │
        ▼ Top 3 candidates → Clinical Report
```

---

## 🏆 Top Candidate Summary

| Rank | Variant | Gene | Type | gnomAD AF | ClinVar | CADD | Phenotype Match |
|------|---------|------|------|-----------|---------|------|-----------------|
| 1 | chr2:166,182,386 G>A | SCN1A | Missense (p.Arg1648His) | 0.000001 | Pathogenic | 28.4 | Seizures ✓ ID ✓ DD ✓ |
| 2 | chr14:23,884,527 C>T | FOXG1 | Nonsense (p.Arg461*) | novel | Likely pathogenic | 42.1 | ID ✓ DD ✓ Dysmorphism ✓ |
| 3 | chrX:153,296,777 G>A | MECP2 | Missense (p.Arg306Cys) | 0.000003 | Pathogenic | 25.6 | Seizures ✓ ID ✓ DD ✓ |

**Final diagnosis:** SCN1A p.Arg1648His — Dravet syndrome (OMIM:607208)  
**ACMG classification:** Pathogenic (PS1, PM2, PP3, PP4)  
**Inheritance:** De novo (parents unaffected, autosomal dominant)

---

## 🗃️ Databases Used

| Database | Use |
|---|---|
| [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) | Variant pathogenicity classification |
| [gnomAD v3.1.2](https://gnomad.broadinstitute.org/) | Population allele frequencies |
| [OMIM](https://www.omim.org/) | Gene–disease relationships |
| [HPO](https://hpo.jax.org/) | Phenotype-to-gene mapping |
| [Ensembl VEP](https://www.ensembl.org/Tools/VEP) | Variant effect prediction |
| [CADD](https://cadd.gs.washington.edu/) | Combined deleteriousness score |
| [dbSNP](https://www.ncbi.nlm.nih.gov/snp/) | Variant ID lookup |

---

## 🌟 Bonus: Shiny Variant Browser

An interactive R Shiny application is included in `/shiny_app/` that allows:
- Filtering candidates by CADD score, gnomAD AF, consequence type
- Sorting and ranking candidates dynamically
- Viewing ACMG evidence tags per variant
- Exporting filtered candidates as CSV

See [`shiny_app/README.md`](shiny_app/README.md) for setup instructions.

---

## 📚 References

1. Claes et al. (2001). *De novo mutations in the sodium-channel gene SCN1A cause severe myoclonic epilepsy of infancy.* Am J Hum Genet.
2. Richards et al. (2015). *Standards and guidelines for the interpretation of sequence variants.* Genet Med. (ACMG/AMP framework)
3. Karczewski et al. (2020). *The mutational constraint spectrum quantified from variation in 141,456 humans.* Nature. (gnomAD)
