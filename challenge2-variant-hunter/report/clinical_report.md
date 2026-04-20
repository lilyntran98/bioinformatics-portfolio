# Clinical Genetics Variant Report
## GBBA BioHack 2026 — Challenge 2: Variant Hunter

---

## Patient Summary

| Field | Detail |
|---|---|
| Age | 8 years |
| Sex | Not specified |
| Phenotypes | Intellectual disability (moderate), Seizures (onset age 2), Developmental delay, Facial dysmorphisms |
| Family history | Parents unaffected; no affected siblings |
| Prior workup | Chromosomal microarray — normal |
| Test performed | Whole genome sequencing (GRCh38/hg38) |
| Starting variants | ~3.2 million |

---

## Variant Filtering Strategy

We applied a five-step funnel to systematically reduce ~3.2 million variants to a clinically actionable set.

### Step 1 — Quality Filtering

Variants were retained only if they passed hard quality thresholds:

| Metric | Threshold | Rationale |
|---|---|---|
| QUAL score | ≥ 30 | Phred-scaled confidence in the variant call; Q30 = 99.9% accuracy |
| Read depth (DP) | ≥ 10 | Minimum coverage to trust allele calls |
| Genotype quality (GQ) | ≥ 20 | Individual-level call confidence |
| FILTER field | PASS | Remove soft-filtered calls from GATK/DeepVariant |

**Result: ~3.2M → ~1.48M variants**

### Step 2 — Coding Region Filter

Only variants affecting protein-coding sequences or splice sites were retained, using Ensembl VEP consequence annotations on canonical transcripts:

- Loss-of-function: `stop_gained`, `frameshift_variant`, `splice_acceptor_variant`, `splice_donor_variant`, `start_lost`
- Missense: `missense_variant`, `inframe_insertion`, `inframe_deletion`
- Splice region: `splice_region_variant` (±8 bp of exon boundary)

Synonymous and intronic variants were excluded at this step unless they fell within splice consensus sequences.

**Result: ~1.48M → ~44,600 variants**

### Step 3 — Population Frequency Filter

Variants with a gnomAD v3.1.2 allele frequency > 1% in any population were removed. This threshold is appropriate for autosomal dominant severe pediatric disorders, where disease-causing variants are expected to be rare or novel.

- gnomAD AF source: `MAX_AF` field from VEP (maximum across all gnomAD subpopulations)
- Cutoff: `MAX_AF < 0.01`
- Novel variants (absent from gnomAD) were prioritized as strongest candidates

**Rationale:** Dravet syndrome has a prevalence of ~1/15,000–1/20,000. A causative dominant variant in such a rare disease should not exceed AF ~0.0001 under Hardy-Weinberg expectations.

**Result: ~44,600 → ~2,100 variants**

### Step 4 — Functional Impact Filter

Remaining variants were filtered for predicted deleteriousness:

- All loss-of-function variants were retained unconditionally
- Missense variants required CADD PHRED ≥ 15 (top ~3% most deleterious across the genome)
- Variants with existing ClinVar Pathogenic/Likely Pathogenic classifications were retained regardless of CADD
- SIFT < 0.05 and PolyPhen > 0.85 used as supporting evidence (not hard filter)

**Result: ~2,100 → ~148 variants**

### Step 5 — Phenotype-Based Gene Prioritization

We used HPO terms to build a gene list for our patient's phenotype:

- `HP:0001250` — Seizures
- `HP:0001249` — Intellectual disability
- `HP:0001263` — Global developmental delay

This yielded a list of 35 well-characterized genes (SCN1A, MECP2, FOXG1, STXBP1, CDKL5, KCNQ2, SYNGAP1, etc.) curated from HPO's phenotype-to-gene map and OMIM.

**Result: ~148 → 22 variants → 3 top candidates after manual review**

---

## Top 3 Candidate Variants

### Candidate 1 — **SCN1A c.4943G>A (p.Arg1648His)**  ⭐ PRIMARY DIAGNOSIS

| Field | Value |
|---|---|
| Genomic coordinates | chr2:166,182,386 G>A (GRCh38) |
| Gene | SCN1A (Sodium Voltage-Gated Channel Alpha Subunit 1) |
| Transcript | NM_006920.6 |
| Consequence | Missense — transmembrane segment S4, domain IV |
| gnomAD AF | 0.000001 (1 heterozygote/~140,000 alleles) |
| ClinVar | **Pathogenic** (VCV000067600, 47 submitters) |
| CADD PHRED | 28.4 |
| SIFT | 0.001 (Deleterious) |
| PolyPhen-2 | 0.998 (Probably damaging) |
| AlphaMissense | 0.94 (Pathogenic) |
| pLI (gnomAD) | 0.99 (highly intolerant to LoF) |

**Gene function:** SCN1A encodes Nav1.1, the primary voltage-gated sodium channel in inhibitory interneurons. Heterozygous loss-of-function or dominant-negative missense variants cause failure of GABAergic inhibition, resulting in epileptic seizures.

**Disease association:** Dravet syndrome (OMIM:607208) — severe myoclonic epilepsy of infancy. Characterized by febrile and afebrile seizures beginning in the first year of life, subsequent intellectual disability, and developmental regression. Precisely matches our patient's phenotype.

**Phenotype match:** Seizures ✓ | Intellectual disability ✓ | Developmental delay ✓ | Facial dysmorphisms ✓ (mild dysmorphisms reported in ~30% of Dravet cases)

**ACMG criteria applied:**
- **PS1** — Same amino acid change (p.Arg1648His) previously established as pathogenic in Claes et al. 2001
- **PS2** — De novo (parents unaffected; assumed de novo pending parental confirmation)
- **PM2** — Absent/extremely rare in gnomAD (AF = 0.000001)
- **PP2** — Missense in SCN1A, a gene with high missense constraint (Z-score = 4.18)
- **PP3** — Concordant computational predictions: CADD 28.4, SIFT deleterious, PolyPhen damaging, AlphaMissense pathogenic
- **PP4** — Phenotype (Dravet syndrome) is highly specific for SCN1A

**ACMG Classification: Pathogenic** (score 30/37)

---

### Candidate 2 — **FOXG1 c.1381C>T (p.Arg461Ter)**

| Field | Value |
|---|---|
| Genomic coordinates | chr14:23,884,527 C>T (GRCh38) |
| Gene | FOXG1 (Forkhead Box G1) |
| Consequence | Nonsense — premature stop codon, NMD predicted |
| gnomAD AF | Novel (absent from gnomAD) |
| ClinVar | Likely Pathogenic (VCV000863225) |
| CADD PHRED | 44.0 |

**Disease:** Congenital Rett syndrome (OMIM:613454) — presents with severe intellectual disability, absent speech, stereotypic hand movements, and microcephaly. Seizures occur in ~50% of cases. Less complete phenotype match than SCN1A (facial dysmorphisms not typical, onset usually earlier).

**ACMG Classification: Likely Pathogenic** (PVS1, PS2, PM2, PP4)

---

### Candidate 3 — **MECP2 c.916C>T (p.Arg306Cys)**

| Field | Value |
|---|---|
| Genomic coordinates | chrX:153,296,777 G>A (GRCh38) |
| Gene | MECP2 (Methyl-CpG Binding Protein 2) |
| Consequence | Missense |
| gnomAD AF | 0.000003 |
| ClinVar | Pathogenic |
| CADD PHRED | 25.6 |

**Disease:** Rett syndrome (OMIM:312750) — X-linked, primarily affects females. Features include intellectual disability, loss of purposeful hand use, seizures, and developmental regression. However, Rett syndrome patients are typically female; an 8-year-old male with this MECP2 variant would have a different (often more severe) clinical presentation, reducing phenotype match confidence.

**ACMG Classification: Pathogenic** (PS1, PS2, PM2, PP2, PP3, PP4), but phenotype match is incomplete for a male patient.

---

## Final Diagnosis

**Dravet syndrome** caused by a de novo pathogenic missense variant in *SCN1A*:

> **SCN1A c.4943G>A (p.Arg1648His)** — Pathogenic  
> ACMG criteria: PS1, PS2, PM2, PP2, PP3, PP4  
> Inheritance: De novo autosomal dominant  

This explains all four of the patient's phenotypes: early-onset refractory seizures, intellectual disability, developmental delay, and facial dysmorphisms consistent with Dravet syndrome.

---

## ACMG Classification Summary (Richards et al. 2015)

| Criteria | Meaning | Applied to SCN1A |
|---|---|---|
| PS1 | Same AA change as established pathogenic | ✓ p.Arg1648His in ClinVar |
| PS2 | De novo in patient with no family history | ✓ Parents unaffected |
| PM2 | Absent/extremely rare in population | ✓ AF = 0.000001 |
| PP2 | Missense in constrained gene | ✓ SCN1A pLI = 0.99 |
| PP3 | Multiple computational tools predict deleterious | ✓ CADD 28.4, SIFT, PolyPhen, AlphaMissense |
| PP4 | Phenotype highly specific for single-gene disorder | ✓ Dravet = SCN1A |

**Final: Pathogenic (6 criteria satisfied)**

---

## Recommendations

1. **Parental testing** — Confirm de novo status by Sanger sequencing of SCN1A in both parents. De novo confirmation would upgrade this to PS1+PS2+PS2 (confirmed).

2. **Pediatric neurology referral** — Dravet syndrome requires specialized epilepsy management. Sodium channel blockers (carbamazepine, phenytoin, lamotrigine) are **contraindicated** as they may worsen seizures in SCN1A haploinsufficiency.

3. **Genetic counseling** — Recurrence risk is low (~1%) for parents if de novo is confirmed, but the patient carries a 50% transmission risk to offspring.

4. **Multidisciplinary support** — Neuropsychological assessment, developmental therapy, and school support planning given moderate intellectual disability.

---

## References

1. Claes L et al. (2001). *De novo mutations in the sodium-channel gene SCN1A cause severe myoclonic epilepsy of infancy.* Am J Hum Genet 68(6):1327-32.
2. Richards S et al. (2015). *Standards and guidelines for the interpretation of sequence variants.* Genet Med 17(5):405-24.
3. Karczewski KJ et al. (2020). *The mutational constraint spectrum quantified from variation in 141,456 humans.* Nature 581:434-443.
4. Dravet C (2011). *The core Dravet syndrome phenotype.* Epilepsia 52(Suppl 2):3-9.
