#!/usr/bin/env bash
# =============================================================================
# 02_annotate_variants.sh
# Annotate patient VCF using Ensembl VEP with gnomAD, ClinVar, CADD plugins
#
# Prerequisites:
#   - VEP installed with GRCh38 cache (~14 GB, download once)
#   - CADD plugin data files
#   - gnomAD VEP plugin
#
# Input:  data/patient_wgs.vcf.gz
# Output: results/vep_annotated.vcf.gz
# =============================================================================

set -euo pipefail

VCF_IN="data/patient_wgs.vcf.gz"
VCF_OUT="results/vep_annotated.vcf"
VEP_CACHE="${HOME}/.vep"
CADD_SNVS="${VEP_CACHE}/Plugins/CADD/whole_genome_SNVs.tsv.gz"
CADD_INDELS="${VEP_CACHE}/Plugins/CADD/gnomad.genomes.r3.0.indel.tsv.gz"

mkdir -p results

# ---------------------------------------------------------------------------
# (First-time only) Download VEP cache
# ---------------------------------------------------------------------------
if [ ! -d "${VEP_CACHE}/homo_sapiens/110_GRCh38" ]; then
    echo ">>> Downloading VEP GRCh38 cache (this takes ~20 min on first run)..."
    vep_install -a cf -s homo_sapiens -y GRCh38 -c "${VEP_CACHE}" --CONVERT
fi

# ---------------------------------------------------------------------------
# Run VEP
# ---------------------------------------------------------------------------
echo ">>> Running Ensembl VEP annotation..."
echo "    Input:  ${VCF_IN}"
echo "    Output: ${VCF_OUT}.gz"

vep \
  --input_file "${VCF_IN}" \
  --output_file "${VCF_OUT}" \
  --format vcf \
  --vcf \
  --everything \
  --assembly GRCh38 \
  --cache \
  --dir_cache "${VEP_CACHE}" \
  --fork 4 \
  --species homo_sapiens \
  --canonical \
  --symbol \
  --hgvs \
  --af_gnomadg \
  --af_gnomade \
  --max_af \
  --check_existing \
  --variant_class \
  --sift b \
  --polyphen b \
  --plugin CADD,"${CADD_SNVS}","${CADD_INDELS}" \
  --plugin AlphaMissense,file="${VEP_CACHE}/Plugins/AlphaMissense/AlphaMissense_hg38.tsv.gz" \
  --filter_common \
  --stats_file "results/vep_stats.html" \
  --warning_file "results/vep_warnings.txt" \
  --no_progress

bgzip "${VCF_OUT}"
tabix -p vcf "${VCF_OUT}.gz"

# ---------------------------------------------------------------------------
# Quick summary
# ---------------------------------------------------------------------------
echo ""
echo ">>> VEP annotation complete."
echo "    Output: results/vep_annotated.vcf.gz"
echo ""
echo "    Consequence summary:"
bcftools +split-vep "results/vep_annotated.vcf.gz" -f "%Consequence\n" -A tab 2>/dev/null \
  | sort | uniq -c | sort -rn | head -15 \
  || echo "    (run bcftools split-vep to inspect consequences)"

echo ""
echo ">>> VEP stats written to: results/vep_stats.html"
echo "    Open in a browser to review annotation coverage."
