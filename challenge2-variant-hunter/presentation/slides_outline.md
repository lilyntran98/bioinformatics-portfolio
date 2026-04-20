# Presentation Outline — Challenge 2: Variant Hunter

## Slide 1: Title
**"Finding the Pathogenic Needle in a 3-Million Variant Haystack"**  
GBBA BioHack 2026 · Challenge 2  
Team: Lily Lee

---

## Slide 2: The Clinical Problem
- 8-year-old with seizures, intellectual disability, developmental delay
- Chromosomal microarray: normal → need higher resolution
- WGS gives us ~3.2 million variants — 99.99% are benign noise
- **Our job:** systematically find the 1 that explains this child's disease

*Key message: This is exactly how real clinical genomics works.*

---

## Slide 3: Pipeline Overview (flowchart from README)
Walk through the 5-step funnel:
1. Quality → 2. Coding → 3. Rare → 4. Functional → 5. Phenotype  
3.2M → 1.48M → 44.6K → 2.1K → 148 → **22 → 3 candidates**

---

## Slide 4: Dataset & Spike-In Strategy
- Base: 1000 Genomes HG00096, chromosomes 2/14/X (GRCh38)
- Spike-in: SCN1A p.Arg1648His — ClinVar Pathogenic, Dravet syndrome
- Why spike-in: ethical access to real clinical WGS is restricted; this is standard in training pipelines
- Demonstrates the pipeline correctly recovers the known pathogenic variant

---

## Slide 5: Top 3 Candidates (table)
| Rank | Variant | Gene | CADD | ClinVar | ACMG |
|---|---|---|---|---|---|
| 1 | chr2:166,182,386 G>A | SCN1A p.Arg1648His | 28.4 | Pathogenic | **Pathogenic** |
| 2 | chr14:23,884,527 C>T | FOXG1 p.Arg461Ter | 44.0 | Likely Path | Likely Pathogenic |
| 3 | chrX:153,296,777 G>A | MECP2 p.Arg306Cys | 25.6 | Pathogenic | Pathogenic (partial match) |

---

## Slide 6: The Verdict — Dravet Syndrome
- **Gene:** SCN1A — Nav1.1 sodium channel, critical for inhibitory interneurons
- **Mechanism:** Haploinsufficiency → loss of GABAergic inhibition → refractory seizures
- **Phenotype match:** Seizures ✓ · ID ✓ · DD ✓ · Dysmorphisms ✓
- **ACMG criteria:** PS1 + PS2 + PM2 + PP2 + PP3 + PP4 = Pathogenic
- **Inheritance:** De novo autosomal dominant

---

## Slide 7: IGV Evidence
- Screenshot: SCN1A chr2:166,182,386 G>A in patient track
- Clear heterozygous call: ~50% alt reads
- Clean coverage, no artifacts

---

## Slide 8: Shiny Variant Browser (demo)
Live demo of `shiny_app/app.R`:
- Interactive filtering sliders (CADD, AF)
- ACMG heatmap
- Phenotype match visualization
- CSV export

*"We built this so future teams and clinicians can interactively explore the candidates without running command-line tools."*

---

## Slide 9: Clinical Implications
- **Do NOT use:** sodium channel blockers (carbamazepine, lamotrigine) — worsen Dravet
- **Consider:** valproate, clobazam, stiripentol (Dravet-specific regimen)
- **Recommend:** parental testing to confirm de novo, genetic counseling
- **Resources:** Dravet Syndrome Foundation, SUDEP risk counseling

---

## Slide 10: Takeaways
- Systematic filtering is essential — 3.2M → 3 candidates in 5 steps
- ACMG framework provides reproducible, evidence-based classification
- Open-source tools (VEP, bcftools, R/Bioconductor) make this accessible
- Phenotype-first thinking drives gene prioritization

**GitHub:** `github.com/[yourusername]/challenge2-variant-hunter`
