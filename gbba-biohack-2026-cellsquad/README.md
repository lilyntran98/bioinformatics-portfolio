# 🧬 GBBA BioHack 2026 - Team Cells Squad

## 🏆 Awards

![Best Presentation] (https://img.shields.io/badge/GBBA%20BioHack%202026-Best%20Presentation-gold?style=for-the-badge)

![Technical Depth] (https://img.shields.io/badge/GBBA%20BioHack%202026-Technical%20Depth-blue?style=for-the-badge)


**Team Members:**
- Lily Lee
- Yash Kiran
- Anushka Shetty

**Challenge:** 
🔍 Challenge 1: Sequence Quality Detective
You're a bioinformatics consultant hired to investigate why a research lab's RNA-seq experiment failed. They sequenced bacterial samples but their analysis pipeline crashed. Your job: diagnose the problem, clean the data, and determine if anything can be salvaged.

---

## 📢 Event Information
This repository was created for **GBBA BioHack 2026**, organized by the **Graduate Biotechnology Business Association (GBBA)** at Northeastern University.

**Stay Connected:**
- 📸 Instagram: [@gbba_neu](https://https://www.instagram.com/gbba_neu/)
- 💼 LinkedIn: [GBBA Northeastern](https://www.linkedin.com/company/gbbaneu/)

---

## About This Repository
This repository will contain our team's submission for GBBA BioHack 2026.

## 🎤 Presentation

**Google Slides:** [View Presentation](https://docs.google.com/presentation/d/14taL3mDlIl9JrqeznElNdrTOwxLQ-0GJ1jj2jcQKDrg/edit?slide=id.g3ba48983ece_0_5#slide=id.g3ba48983ece_0_5)

*5-minute presentation covering dataset selection, problems identified, cleaning results, and recommendations*

---

## 🔬 Our Investigation

### Dataset Selected
**Study:** *E. coli* response to bicarbonate as antibiotics potentiator  
**BioProject:** [PRJNA433855](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA433855)  
**Platform:** Illumina HiSeq 2000 (2018)  
**Samples:** 4 (2 bicarbonate treatment, 2 NaCl control)  
**Total Data:** 22.5M read pairs (~2.3 GB)

**Why this dataset?**
- Older HiSeq 2000 technology (known for quality issues)
- Real-world data with authentic problems
- Appropriate size for hackathon timeframe
- Biological replicates for quality comparison

---

## 🔍 Problems Identified

### Critical Issue: Very High Duplication
- **Found:** 68-71% sequence duplication across all samples
- **Expected:** 20-40% for RNA-seq
- **Root Cause:** Over-amplification during library preparation (too many PCR cycles)
- **Impact:** Severely reduced library diversity

### Moderate Issue: Adapter Contamination
- **Found:** ~4% adapter content
- **Root Cause:** Incomplete adapter trimming or short insert sizes
- **Impact:** 4% of reads unusable without cleaning

### Positive Finding: Excellent Base Quality
- **Found:** Q35-38 average quality scores
- **GC Content:** 48-50% (correct for *E. coli*)
- **Consistency:** All samples showed similar quality profiles

---

## 🧹 Our Solution

### Cleaning Pipeline (fastp)
1. **Adapter removal** - Auto-detect and trim paired-end adapters
2. **Quality filtering** - Remove bases <Q20
3. **Length filtering** - Discard reads <50bp after trimming
4. **Aggressive deduplication** - Remove PCR duplicates

---

## 📊 Results - Exceeded All Expectations!

### Overall Performance

| Metric | Before Cleaning | After Cleaning | Improvement |
|--------|----------------|----------------|-------------|
| **Total Reads** | 22.5M pairs | **13.7M pairs** | **61.2% retained** ✓ |
| **Duplication** | 68-71% | **35-37%** | **-32%** ✓ |
| **Adapter Content** | ~4% | **<0.1%** | **-99%** ✓ |
| **Mean Quality** | Q37 | **Q37** | **Maintained** ✓ |
| **Q30 Bases** | N/A | **>95%** | **Excellent** ✓ |

### 🏆 Major Achievement
**Predicted retention:** 25-35% (worst case for 70% duplication)  
**Actual retention:** 61.2%  
**EXCEEDED prediction by >25 percentage points!**

### Per-Sample Results

| Sample | Condition | Input | Output | Retention | Dup: Before → After |
|--------|-----------|-------|--------|-----------|---------------------|
| NaHCO3rep1 | Bicarbonate | 5.5M | 3.3M | 59.8% | 70.9% → 37.8% |
| NaHCO3rep2 | Bicarbonate | 5.2M | 3.2M | 60.9% | 69.7% → ~37% |
| NaClrep1 | NaCl control | 6.0M | 3.7M | 61.7% | 68.6% → ~36% |
| NaClrep2 | NaCl control | 5.7M | 3.6M | 62.3% | 69.3% → 35.5% |

---

## ✅ Final Recommendation

### **PROCEED WITH DATA**

**Justification:**
- ✅ Adequate sequencing depth: 13.7M read pairs (avg 3.4M per sample)
- ✅ Excellent quality maintained: Q37 throughout
- ✅ Acceptable duplication: 37% (within 20-40% range)
- ✅ No adapter contamination: <0.1%
- ✅ All 4 samples consistently cleaned
- ✅ Suitable for differential expression analysis

**Data Quality:** Analysis-ready for publication-quality research

---

## 📁 Repository Structure
```
gbba-biohack-2026-cellsquad/
├── README.md                          # This file
└── challenge1-qc-detective/
    ├── README.md                      # Challenge overview
    ├── data/
    │   ├── download_instructions.md   # How to obtain data
    │   ├── raw/                       # Raw FASTQ files (22.5M pairs, not in repo)
    │   └── cleaned/                   # Cleaned FASTQ files (13.7M pairs, not in repo)
    ├── scripts/
    │   ├── 01_download_data_ena.sh    # Download from ENA
    │   ├── 02_run_fastqc.sh           # Initial QC
    │   ├── 03_clean_reads.sh          # Data cleaning (fastp)
    │   ├── 04_post_qc.sh              # Post-cleaning QC
    │   └── 05_compare_qc.sh           # Before/after comparison
    ├── results/
    │   ├── multiqc_raw/               # Initial QC report
    │   ├── multiqc_cleaned/           # Post-cleaning QC report
    │   └── fastp_reports/             # Cleaning reports (4 samples)
    ├── report/
    │   └── diagnosis_report.md        # Comprehensive analysis
    └── presentation/
        ├── PRESENTATION_LINK.md       # Google Slides link
        └── slides.md                  # Presentation content
```

---

## 🛠️ Tools & Technologies

### Analysis Platform
- **HPC Cluster:** Northeastern Discovery Cluster
- **Environment:** Conda + module system
- **Version Control:** Git + GitHub

### Bioinformatics Tools
- **SRA Toolkit** - Data download from NCBI/ENA
- **FastQC v0.12.1** - Per-sample quality control
- **MultiQC v1.14** - Aggregate quality reports
- **fastp v0.23.4** - Read cleaning and deduplication

---

## 📈 Key Achievements

✅ **Exceeded retention prediction by >25%** (61% vs 25-35% predicted)  
✅ **Reduced duplication by 33%** (69% → 37%)  
✅ **Eliminated 99% of adapter contamination** (4% → <0.1%)  
✅ **Maintained excellent quality** (Q37 throughout)  
✅ **All samples salvageable** (4/4 suitable for analysis)  
✅ **Complete documentation** (reproducible pipeline)

---

## 📚 Deliverables Completed

### Required Deliverables (100 points)
- [x] **Initial QC Report (25 pts)** - MultiQC aggregate of raw data
- [x] **Problem Diagnosis (30 pts)** - 8-section comprehensive report identifying root causes
- [x] **Data Cleaning & Results (30 pts)** - Cleaned FASTQ files + before/after comparison
- [x] **Recommendations (10 pts)** - Evidence-based assessment with future QC checkpoints
- [x] **Code & Documentation (5 pts)** - Well-commented scripts, reproducible workflow

### Bonus Deliverables
- [x] **Statistical comparison** - Documented 61% retention vs 25-35% prediction
- [x] **Multiple analysis reports** - FastQC, MultiQC, fastp reports
- [x] **Comprehensive documentation** - Complete README, download instructions, analysis reports
- [x] **Professional presentation** - 5-minute Google Slides with actual results
---

## 🔗 Quick Links

- **Challenge Directory:** [challenge1-qc-detective/](challenge1-qc-detective/)
- **Diagnosis Report:** [Full Analysis](challenge1-qc-detective/report/diagnosis_report.md)
- **Download Instructions:** [Data Acquisition Guide](challenge1-qc-detective/data/download_instructions.md)
- **Scripts:** [Reproducible Pipeline](challenge1-qc-detective/scripts/)
- **QC Reports:**
  - [Raw Data QC](challenge1-qc-detective/results/multiqc_raw/multiqc_report.html)
  - [Cleaned Data QC](challenge1-qc-detective/results/multiqc_cleaned/multiqc_report.html)

---

## 🎓 Learning Outcomes

Through this challenge, we demonstrated:
- **Diagnostic Skills:** Identified critical quality issues in real-world sequencing data
- **Root Cause Analysis:** Determined over-amplification as primary issue
- **Problem Solving:** Implemented effective cleaning strategy exceeding predictions
- **Technical Expertise:** Used industry-standard bioinformatics tools (FastQC, MultiQC, fastp)
- **Documentation:** Created comprehensive, reproducible analysis pipeline
- **Communication:** Presented findings clearly with evidence-based recommendations

  
<img width="500" height="500" alt="BioHack" src="https://github.com/user-attachments/assets/beae234f-f0b2-41c6-a2e8-357f3e92476d" />

