# Challenge 1: Sequence Quality Detective - Problem Diagnosis

**Team:** CellSquad  
**Date:** January 22-23, 2026  
**Dataset:** E. coli bicarbonate treatment (PRJNA433855)

---

## 1. Dataset Information

**Study:** The antibiotics potentiator bicarbonate causes upregulation of tryptophanase and iron acquisition proteins in *Escherichia coli*

**BioProject:** PRJNA433855  
**Platform:** Illumina HiSeq 2000  
**Release Date:** November 2, 2018

**Samples Analyzed:**

| Sample | SRA Accession | Condition | Read Pairs | Size |
|--------|---------------|-----------|------------|------|
| NaHCO3rep1 | SRR6710205 | Bicarbonate treatment | 5,520,708 | 561 MB |
| NaHCO3rep2 | SRR6710206 | Bicarbonate treatment | 5,184,530 | 535 MB |
| NaClrep1 | SRR6710207 | NaCl control | 5,998,161 | 610 MB |
| NaClrep2 | SRR6710208 | NaCl control | 5,748,488 | 585 MB |

**Total Data:** 22.4M read pairs (~2.3 GB compressed)

---

## 2. Quality Issues Identified

### A. Critical Issues

#### Issue 1: Very High Sequence Duplication ⚠️

**Severity:** CRITICAL

**Affected Samples:** ALL samples (68-71%)

**Metrics:**
- NaHCO3rep1: 70.9% duplicated
- NaHCO3rep2: 69.7% duplicated
- NaClrep1: 68.6% duplicated
- NaClrep2: 69.3% duplicated

**Evidence:** FastQC General Statistics table shows 68-71% duplication across all samples, compared to expected 20-40% for RNA-seq datasets.

**Root Cause Analysis:**

The uniformly high duplication across all samples (68-71%) indicates a systematic library preparation issue rather than sample-specific problems. Most likely causes:

1. **Over-amplification during library prep:** The library was likely prepared with 15-20 PCR cycles instead of the optimal 8-12 cycles for RNA-seq
2. **Low RNA input:** Insufficient starting material required excessive amplification to generate enough library for sequencing
3. **Platform characteristics:** HiSeq 2000 (2018) represents older sequencing technology known for higher duplication rates compared to modern platforms

**Impact:**
- Severely reduces effective library diversity (only ~30% unique sequences)
- May bias gene expression quantification toward highly amplified fragments
- Inflates read counts for abundant transcripts
- Limits ability to accurately quantify lowly expressed genes
- Reduces statistical power for differential expression analysis

---

### B. Moderate Issues

#### Issue 2: Adapter Contamination

**Severity:** MODERATE

**Affected Samples:** ALL samples (~4% maximum)

**Metrics:**
- Maximum adapter content: ~4% at positions 30-100
- Multiple adapter types detected (Illumina Universal, TruSeq)
- Contamination starts at position ~10 and plateaus by position 40

**Evidence:** Adapter Content plot shows multiple adapter lines (primarily teal line at 4%, with others at 1.5-2%) rising throughout read length.

**Root Cause Analysis:**

Two possible explanations:

1. **Incomplete adapter trimming:** Library preparation protocol did not fully remove adapter sequences before sequencing
2. **Short insert sizes:** Some RNA fragments were shorter than the 150bp read length, causing reads to extend past the insert into the adapter sequence

Given that this was HiSeq 2000 in 2018, older adapter removal protocols were likely less efficient than current methods.

**Impact:**
- Approximately 4% of sequences contain adapter artifacts
- Adapter-containing reads will fail to map to the reference genome
- Reduces usable data by ~4%
- May cause spurious mappings if adapters align by chance

---

### C. Minor Issues

#### Issue 3: Overrepresented Sequences

**Severity:** MINOR (acceptable for RNA-seq)

**Affected Samples:** ALL samples (5-7.5% total)

**Metrics:**
- Top overrepresented sequence: <2.5% per sample
- Total overrepresented sequences: 5-7.5% per sample
- Consistent pattern across all samples

**Evidence:** Overrepresented Sequences plot shows uniform low-level overrepresentation across all samples.

**Root Cause Analysis:**

This level of sequence overrepresentation is **normal and expected** for RNA-seq data:

1. **Highly expressed genes:** Some transcripts (e.g., ribosomal proteins, housekeeping genes, stress response genes) are naturally very abundant
2. **Biological signal:** In E. coli under stress conditions, certain genes are upregulated 100-1000 fold
3. **Not contamination:** No single sequence dominates, ruling out technical contamination

**Impact:**
- Minimal - this is biological signal, not technical artifact
- Does not require correction or removal
- Expected behavior for transcriptomic sequencing

---

## 3. Positive Quality Aspects

Despite the critical duplication issue, several quality metrics are excellent:

✅ **High Per-Sequence Quality Scores**
- Mean quality: Q35-38 across all samples
- Per Sequence Quality Scores plot shows sharp peak in green zone (Q>30)
- Very few reads below Q20
- **This provides a strong foundation for analysis after cleaning**

✅ **Correct GC Content**
- Observed: 48-50% across samples
- Expected for *E. coli*: ~50-51%
- **Indicates no major contamination from other organisms**

✅ **No Per-Base Quality Degradation**
- Quality remains consistently high throughout read length
- No characteristic HiSeq quality drop-off at read ends
- **Unusual for HiSeq 2000 - actually better than expected**

✅ **Consistent Sample Quality**
- All samples show similar quality profiles
- No outlier samples requiring exclusion
- **Good experimental consistency**

✅ **Adequate Sequencing Depth**
- 5.2-6.0M read pairs per sample
- After cleaning, expect 1.5-2M usable reads per sample
- **Sufficient for differential expression analysis**

---

## 4. Sample-by-Sample Comparison

| Sample | Reads (M) | %Dups | %GC | Mean Q | Adapters | Overall Assessment |
|--------|-----------|-------|-----|--------|----------|-------------------|
| NaHCO3rep1 | 5.5 | 70.9% | 49% | Q37 | ~4% | Fair (highest dup) |
| NaHCO3rep2 | 5.2 | 69.7% | 50% | Q37 | ~4% | Fair |
| NaClrep1 | 6.0 | 68.6% | 49% | Q37 | ~4% | Fair (best overall) |
| NaClrep2 | 5.7 | 69.3% | 48% | Q37 | ~4% | Fair |

**Key Observations:**

**Most affected sample:** NaHCO3rep1 (70.9% duplication, lowest read count)

**Best quality sample:** NaClrep1 (68.6% duplication, highest read count)

**Treatment vs Control:** Bicarbonate treatment samples show slightly higher duplication (70.3% avg) compared to NaCl controls (68.9% avg). This 1.4% difference is minor and likely due to random variation rather than biological factors.

**Sample consistency:** All samples are remarkably similar in quality metrics, indicating:
- Uniform library preparation protocol applied to all samples
- Batch processing effect (all samples prepared together)
- No sample-specific quality failures

---

## 5. Root Cause Summary

### Primary Issue: Library Preparation Protocol

The critical 68-71% duplication across ALL samples indicates a **systematic protocol issue**:

**Evidence:**
- Uniform duplication rate across all samples (68.6-70.9%)
- No sample-specific variation that would indicate individual sample problems
- All samples prepared in same batch with same protocol

**Most Likely Cause:**
- Library preparation used **15-20 PCR cycles** instead of optimal 8-12 cycles
- This level of over-amplification would produce exactly this duplication pattern
- Possibly compensating for low RNA input or inefficient library prep

**Alternative Explanations (less likely):**
- Low RNA quality/quantity requiring extra amplification
- Inefficient cDNA synthesis requiring compensation
- Protocol optimized for older platform (HiSeq 2000) rather than current standards

### Secondary Issue: Adapter Trimming

The consistent 4% adapter contamination indicates:
- Standard for HiSeq 2000 era (2018)
- Older adapter removal protocols less efficient than modern methods
- Not unusual for this platform/timeframe

### Platform Factor: HiSeq 2000 (2018)

This older sequencing platform was **intentionally selected** for this challenge because:
- Known for higher duplication rates than modern platforms
- Provides realistic quality issues to diagnose
- Represents real-world data analysts encounter with archived datasets

---

## 6. Recommendations

### Should this data be used?

**Decision:** ⚠️ **PROCEED WITH CAUTION** (after aggressive cleaning and deduplication)

**Justification:**

**Reasons to proceed:**
- High base quality (Q37) provides excellent foundation
- Adapter contamination is cleanable
- Duplication can be addressed through computational deduplication
- Read depth is adequate even after expected 65-70% loss to deduplication
- All samples show similar quality (good for comparative analysis)

**Reasons for caution:**
- Will lose ~70% of reads to deduplication
- May reduce power to detect lowly expressed genes
- Duplication pattern may not be fully random (could bias results)

**Overall:** Data is salvageable but requires aggressive cleaning and awareness of limitations in final analysis.

---

### Cleaning Strategy

**Step 1: Adapter Removal**
- Tool: fastp
- Parameters: `--detect_adapter_for_pe`
- Expected outcome: Remove ~4% contaminated reads

**Step 2: Quality Filtering**
- Tool: fastp
- Parameters: `--qualified_quality_phred 20`, `--length_required 50`
- Expected outcome: Remove low-quality bases/reads

**Step 3: Deduplication**
- Tool: fastp
- Parameters: `--dedup` (removes PCR duplicates)
- Expected outcome: Reduce duplication from 70% to <10%

**Expected Final Metrics:**
- Reads retained: ~25-35% of original (5.6-7.8M total reads)
- Final duplication rate: <10%
- Adapter content: <0.1%
- Mean quality: Q>30 throughout

**This is acceptable for differential expression analysis with appropriate statistical methods.**

---

### Salvageable Samples

**All samples are salvageable:**

✅ **NaClrep1** - Best quality (68.6% dup, 6.0M reads)  
✅ **NaClrep2** - Good quality (69.3% dup, 5.7M reads)  
✅ **NaHCO3rep2** - Good quality (69.7% dup, 5.2M reads)  
✅ **NaHCO3rep1** - Acceptable (70.9% dup, 5.5M reads)

**Recommendation:** Include all four samples in analysis
- Quality differences are minimal (1.4% range)
- Statistical power benefits from including all replicates
- After cleaning, quality will be comparable

---

### If Resequencing Were Needed

**What to do differently:**

1. **Optimize PCR cycles**
   - Current: Likely 15-20 cycles
   - Recommended: 8-12 cycles
   - Test: Titrate cycle number for each RNA input level

2. **Increase RNA input**
   - Use more starting material to reduce amplification needs
   - Minimum: 100ng total RNA for standard protocols
   - Consider concentration methods if sample limited

3. **Improve adapter removal**
   - Use enzymatic adapter removal (e.g., USER enzyme)
   - Increase adapter removal incubation time
   - Size-select library to avoid adapter dimers

4. **Use modern sequencing platform**
   - NovaSeq or NextSeq for lower duplication rates
   - Better chemistry reduces PCR bias
   - Higher throughput for better coverage

5. **Quality checkpoints**
   - Bioanalyzer/TapeStation after library prep
   - qPCR quantification before sequencing
   - Pilot run (1M reads) to check duplication before full sequencing

---

### Future QC Checkpoints

**1. After RNA extraction:**
- Check RNA integrity (RIN score >7)
- Quantify accurately (Qubit, not NanoDrop)
- Assess 260/280 and 260/230 ratios

**2. After library preparation:**
- Bioanalyzer/TapeStation to check library size distribution
- Should see peak at ~300-400 bp for 150bp PE sequencing
- Verify no adapter dimers (~120 bp peak)

**3. After adapter removal:**
- Test subset of library with FastQC
- Verify <1% adapter contamination
- If >1%, repeat adapter removal step

**4. Before sequencing:**
- Pilot sequencing (1M reads) to check duplication
- If >30% duplication, adjust protocol before full run
- Document PCR cycle number for reproducibility

**5. Immediately after data receipt:**
- Run FastQC within 24 hours
- Check duplication, adapters, quality scores
- Flag issues before data ages

---

## 7. Data Cleaning Implementation

**Command to be executed:**
```bash
# Using fastp for all-in-one cleaning
fastp \
    -i raw_R1.fastq.gz \
    -I raw_R2.fastq.gz \
    -o clean_R1.fastq.gz \
    -O clean_R2.fastq.gz \
    --detect_adapter_for_pe \
    --qualified_quality_phred 20 \
    --length_required 50 \
    --dedup \
    --dup_calc_accuracy 6 \
    --thread 4 \
    --html sample_report.html
```

**Expected outcomes:**
- Adapter removal: 96% reads retained (~4% loss)
- Quality filtering: 98% of remaining reads retained (~2% loss)
- Deduplication: 30% of remaining reads retained (~70% loss)
- **Overall retention: ~28-30% of original reads**

---

## 8. Limitations and Caveats

**Analysis Limitations:**

1. **Reduced statistical power:** Lower read counts after deduplication may limit ability to detect lowly expressed genes

2. **Potential quantification bias:** If duplication is non-random, expression estimates may be biased toward certain transcripts

3. **Uncertainty in deduplication:** Cannot perfectly distinguish biological duplicates (naturally identical molecules) from PCR duplicates

**Recommended approach:**
- Use robust statistical methods (DESeq2, edgeR) that account for technical variation
- Focus on highly expressed genes and large fold-changes
- Validate key findings with qRT-PCR
- Report duplication issue in methods section of any publication

---

## 9. Conclusions

**Summary:** This E. coli RNA-seq dataset suffers from **critical sequence duplication (68-71%)** due to over-amplification during library preparation, along with **moderate adapter contamination (4%)**. However, the data maintains **excellent base quality (Q37)** and is **salvageable through aggressive cleaning and deduplication**.

**Key Findings:**
1. ⚠️ Critical: 70% duplication from excessive PCR
2. ⚠️ Moderate: 4% adapter contamination  
3. ✅ Excellent: Q37 mean quality scores
4. ✅ Good: Consistent quality across samples

**Recommendation:** **Proceed with data analysis** after implementing comprehensive cleaning pipeline (adapter removal, quality filtering, deduplication). Expected to retain 25-35% of reads with high quality suitable for differential expression analysis.

**Lesson learned:** This dataset exemplifies why QC should be performed **immediately after sequencing** rather than after library prep completion. Early detection would have allowed protocol optimization before committing full sequencing resources.

---

**Report prepared by:** Team CellSquad  
**Date:** January 23, 2026  
**Tools used:** FastQC v0.12.1, MultiQC v1.14  
**Analysis platform:** NEU Explorer Cluster
