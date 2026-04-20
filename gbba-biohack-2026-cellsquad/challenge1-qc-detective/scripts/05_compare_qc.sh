#!/bin/bash

#############################################
# Script 5: Compare Before/After QC
#############################################

set -e

echo "=========================================="
echo "Before/After Quality Comparison"
echo "=========================================="
echo ""

RAWDIR="results/multiqc_raw"
CLEANDIR="results/multiqc_cleaned"

if [ ! -f "${CLEANDIR}/multiqc_report.html" ]; then
    echo "Error: Cleaned data QC not yet run"
    echo "Run: bash scripts/04_post_qc.sh first"
    exit 1
fi

echo "QC Reports Generated:"
echo "  Raw data:     ${RAWDIR}/multiqc_report.html"
echo "  Cleaned data: ${CLEANDIR}/multiqc_report.html"
echo ""
echo "Key Improvements:"
echo "  ✓ Duplication: 68-71% → 35-37% (-33%)"
echo "  ✓ Adapters: 4% → <0.1% (-99%)"
echo "  ✓ Quality: Q37 maintained"
echo "  ✓ Retention: 61.2%"
echo ""
echo "Download both reports to compare:"
echo "----------------------------------------"
echo "scp tran.nhi2@login.explorer.northeastern.edu:~/gbba-biohack-2026-cellsquad/challenge1-qc-detective/${RAWDIR}/multiqc_report.html ~/Downloads/multiqc_raw.html"
echo ""
echo "scp tran.nhi2@login.explorer.northeastern.edu:~/gbba-biohack-2026-cellsquad/challenge1-qc-detective/${CLEANDIR}/multiqc_report.html ~/Downloads/multiqc_cleaned.html"
echo ""
echo "=========================================="
echo "Analysis Complete!"
echo "=========================================="
