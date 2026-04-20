#!/usr/bin/env bash
# =============================================================================
# 01_download_and_spikein.sh
# Download 1000 Genomes WGS VCF (chr2, chr14, chrX) for sample HG00096,
# then spike in a known pathogenic SCN1A variant to simulate a clinical case.
#
# Output: data/patient_wgs.vcf.gz (+ .tbi index)
# =============================================================================

set -euo pipefail

mkdir -p data results

# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------
SAMPLE="HG00096"
BASE_URL="ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/data_collections/1000G_2504_high_coverage/working/20201028_3202_phased"

# Chromosomes covering SCN1A (chr2), FOXG1 (chr14), MECP2 (chrX)
CHROMS=("chr2" "chr14" "chrX")

# ClinVar-confirmed pathogenic variant to spike in
# SCN1A p.Arg1648His — Dravet syndrome (ClinVar: VCV000067600)
SPIKE_CHROM="chr2"
SPIKE_POS="166182386"
SPIKE_REF="G"
SPIKE_ALT="A"
SPIKE_ID="rs121917885"
SPIKE_QUAL="999"
SPIKE_INFO="AF=0.5;SPIKEIN=SCN1A_DravetSyndrome;ClinVar=Pathogenic;OMIM=607208"

# ---------------------------------------------------------------------------
# Step 1: Download and extract single-sample VCF per chromosome
# ---------------------------------------------------------------------------
echo ">>> Downloading 1000G VCFs for sample ${SAMPLE}..."

VCF_FILES=()
for CHR in "${CHROMS[@]}"; do
    REMOTE_FILE="${BASE_URL}/CCDG_14151_B01_GRM_WGS_2020-08-05_${CHR}.filtered.shapeit2-duohmm-phased.vcf.gz"
    LOCAL_FILE="data/${CHR}_1000G.vcf.gz"

    if [ ! -f "${LOCAL_FILE}" ]; then
        echo "  Downloading ${CHR}..."
        wget -q -O "${LOCAL_FILE}" "${REMOTE_FILE}"
        tabix -p vcf "${LOCAL_FILE}"
    else
        echo "  ${CHR} already downloaded, skipping."
    fi

    # Extract single sample
    SAMPLE_FILE="data/${CHR}_${SAMPLE}.vcf.gz"
    bcftools view -s "${SAMPLE}" "${LOCAL_FILE}" -O z -o "${SAMPLE_FILE}"
    tabix -p vcf "${SAMPLE_FILE}"
    VCF_FILES+=("${SAMPLE_FILE}")
done

# ---------------------------------------------------------------------------
# Step 2: Merge chromosomes into one VCF
# ---------------------------------------------------------------------------
echo ">>> Merging chromosomes..."
bcftools concat "${VCF_FILES[@]}" --allow-overlaps -O z -o data/base_sample.vcf.gz
tabix -p vcf data/base_sample.vcf.gz

TOTAL_VARIANTS=$(bcftools stats data/base_sample.vcf.gz | grep "^SN" | grep "number of SNPs" | awk '{print $4}')
echo "  Total variants in base VCF: ${TOTAL_VARIANTS}"

# ---------------------------------------------------------------------------
# Step 3: Spike in the pathogenic SCN1A variant
# ---------------------------------------------------------------------------
echo ">>> Spiking in pathogenic variant: SCN1A ${SPIKE_CHROM}:${SPIKE_POS} ${SPIKE_REF}>${SPIKE_ALT}"

# Create a single-line VCF for the spike-in variant
cat > /tmp/spikein.vcf <<EOF
##fileformat=VCFv4.2
##FILTER=<ID=PASS,Description="All filters passed">
##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">
##INFO=<ID=SPIKEIN,Number=1,Type=String,Description="Spiked-in variant annotation">
##INFO=<ID=ClinVar,Number=1,Type=String,Description="ClinVar classification">
##INFO=<ID=OMIM,Number=1,Type=String,Description="OMIM disease ID">
##FORMAT=<ID=GT,Number=1,Type=String,Description="Genotype">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	${SAMPLE}
${SPIKE_CHROM}	${SPIKE_POS}	${SPIKE_ID}	${SPIKE_REF}	${SPIKE_ALT}	${SPIKE_QUAL}	PASS	${SPIKE_INFO}	GT	0/1
EOF

bgzip /tmp/spikein.vcf
tabix -p vcf /tmp/spikein.vcf.gz

# Merge base VCF with spike-in
bcftools concat data/base_sample.vcf.gz /tmp/spikein.vcf.gz \
    --allow-overlaps \
    -O z -o data/patient_wgs_unsorted.vcf.gz

bcftools sort data/patient_wgs_unsorted.vcf.gz -O z -o data/patient_wgs.vcf.gz
tabix -p vcf data/patient_wgs.vcf.gz

# ---------------------------------------------------------------------------
# Step 4: Verify spike-in
# ---------------------------------------------------------------------------
echo ">>> Verifying spike-in..."
FOUND=$(bcftools view data/patient_wgs.vcf.gz "${SPIKE_CHROM}:${SPIKE_POS}-${SPIKE_POS}" | grep -v "^#" | wc -l)

if [ "${FOUND}" -ge 1 ]; then
    echo "  ✓ Spike-in variant confirmed present in patient_wgs.vcf.gz"
else
    echo "  ✗ ERROR: Spike-in variant not found — check coordinates"
    exit 1
fi

FINAL_COUNT=$(bcftools stats data/patient_wgs.vcf.gz | grep "number of records" | awk '{print $4}')
echo ""
echo ">>> Done. Final VCF: data/patient_wgs.vcf.gz"
echo "    Total records: ${FINAL_COUNT}"
echo "    Spike-in: SCN1A p.Arg1648His (Dravet syndrome, Pathogenic)"
