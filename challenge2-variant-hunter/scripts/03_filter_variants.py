#!/usr/bin/env python3
"""
03_filter_variants.py
=====================
Multi-step variant filtering pipeline for clinical variant prioritization.

Filtering funnel:
  1. Quality filter    (QUAL, DP, GQ)
  2. Coding regions    (exonic + splice site consequences)
  3. Population freq   (gnomAD AF < threshold)
  4. Functional impact (LoF + high-CADD missense)
  5. Phenotype genes   (HPO-derived seizure/ID/DD gene list)

Usage:
    python scripts/03_filter_variants.py \\
        --input  results/vep_annotated.vcf.gz \\
        --output results/filtered_variants.vcf \\
        --max_af 0.01 \\
        --min_cadd 15 \\
        --hpo_genes data/hpo_seizure_id_genes.txt
"""

import argparse
import gzip
import sys
import os
from pathlib import Path

# ---------------------------------------------------------------------------
# HPO gene list — curated from HPO terms:
#   HP:0001250 (Seizures), HP:0001249 (Intellectual disability),
#   HP:0001263 (Developmental delay)
# Source: https://hpo.jax.org/app/browse/term/HP:0001250
# ---------------------------------------------------------------------------
HPO_SEIZURE_ID_GENES = {
    # Dravet / SCN family
    "SCN1A", "SCN2A", "SCN8A", "SCN1B",
    # GABA-related
    "GABRA1", "GABRB2", "GABRB3", "GABRD",
    # Potassium channels
    "KCNQ2", "KCNQ3", "KCNT1", "KCNA2",
    # Synaptic
    "STXBP1", "SYNGAP1", "SHANK3", "NRXN1",
    # Transcription factors / chromatin
    "FOXG1", "MEF2C", "CDKL5", "ARX",
    # Rett / X-linked
    "MECP2", "PCDH19",
    # mTOR pathway
    "TSC1", "TSC2", "MTOR", "DEPDC5",
    # Metabolic
    "ALDH7A1", "PNPO", "SLC2A1",
    # Other well-characterized
    "PCDH19", "HCN1", "GRIN2A", "GRIN2B", "SPTAN1",
    "WDR45", "RORB", "CASK", "DYRK1A", "MBD5",
}

# ---------------------------------------------------------------------------
# VEP consequence severity — ordered most to least severe
# LoF consequences trigger retention regardless of CADD
# ---------------------------------------------------------------------------
LOF_CONSEQUENCES = {
    "stop_gained",
    "frameshift_variant",
    "splice_acceptor_variant",
    "splice_donor_variant",
    "start_lost",
    "stop_lost",
    "transcript_ablation",
    "transcript_amplification",
}

MISSENSE_CONSEQUENCES = {
    "missense_variant",
    "inframe_insertion",
    "inframe_deletion",
    "protein_altering_variant",
}

CODING_CONSEQUENCES = LOF_CONSEQUENCES | MISSENSE_CONSEQUENCES | {
    "splice_region_variant",
    "synonymous_variant",         # retained for completeness but deprioritized
    "stop_retained_variant",
}


def parse_vep_csq(csq_string: str, csq_fields: list[str]) -> list[dict]:
    """Parse VEP CSQ INFO field into a list of transcript annotation dicts."""
    transcripts = []
    for entry in csq_string.split(","):
        values = entry.split("|")
        if len(values) != len(csq_fields):
            continue
        transcripts.append(dict(zip(csq_fields, values)))
    return transcripts


def get_csq_field_order(vcf_header_lines: list[str]) -> list[str]:
    """Extract VEP CSQ field names from VCF header."""
    for line in vcf_header_lines:
        if line.startswith("##INFO=<ID=CSQ"):
            # Format: ##INFO=<ID=CSQ,...,Description="...Format: FIELD1|FIELD2|...">
            desc_start = line.find("Format: ") + 8
            desc_end = line.rfind('"')
            fields_str = line[desc_start:desc_end]
            return fields_str.split("|")
    return []


def safe_float(value: str, default: float = 0.0) -> float:
    try:
        return float(value) if value not in (".", "", "NA") else default
    except (ValueError, TypeError):
        return default


def parse_info(info_str: str) -> dict:
    """Parse VCF INFO field into a dict."""
    result = {}
    for field in info_str.split(";"):
        if "=" in field:
            k, v = field.split("=", 1)
            result[k] = v
        else:
            result[field] = True
    return result


class VariantFilter:
    def __init__(self, args):
        self.args = args
        self.counters = {
            "total": 0,
            "pass_quality": 0,
            "pass_coding": 0,
            "pass_frequency": 0,
            "pass_functional": 0,
            "pass_phenotype": 0,
        }
        self.csq_fields = []
        self.kept_variants = []

    def run(self):
        print("=" * 60)
        print("Variant Filtering Pipeline")
        print("=" * 60)

        opener = gzip.open if str(self.args.input).endswith(".gz") else open
        header_lines = []

        with opener(self.args.input, "rt") as vcf_in:
            for line in vcf_in:
                line = line.rstrip("\n")

                if line.startswith("##"):
                    header_lines.append(line)
                    if "ID=CSQ" in line:
                        self.csq_fields = get_csq_field_order(header_lines)
                    continue

                if line.startswith("#CHROM"):
                    header_lines.append(line)
                    continue

                self.counters["total"] += 1
                fields = line.split("\t")
                if len(fields) < 8:
                    continue

                chrom, pos, vid, ref, alt, qual, filt, info_str = fields[:8]
                fmt_fields = fields[8] if len(fields) > 8 else ""
                sample = fields[9] if len(fields) > 9 else ""

                info = parse_info(info_str)

                # ── Step 1: Quality filter ───────────────────────────────
                qual_val = safe_float(qual, 0.0)
                dp = safe_float(info.get("DP", "0"), 0.0)

                # Extract GQ from FORMAT/sample if available
                gq = 0.0
                if fmt_fields and sample:
                    fmt_keys = fmt_fields.split(":")
                    smp_vals = sample.split(":")
                    fmt_dict = dict(zip(fmt_keys, smp_vals))
                    gq = safe_float(fmt_dict.get("GQ", "0"), 0.0)

                if qual_val < self.args.min_qual:
                    continue
                if dp > 0 and dp < self.args.min_dp:
                    continue
                if gq > 0 and gq < self.args.min_gq:
                    continue
                self.counters["pass_quality"] += 1

                # ── Step 2: Coding consequence filter ────────────────────
                csq_raw = info.get("CSQ", "")
                if not csq_raw or not self.csq_fields:
                    continue

                transcripts = parse_vep_csq(csq_raw, self.csq_fields)
                canonical = [t for t in transcripts if t.get("CANONICAL") == "YES"]
                if not canonical:
                    canonical = transcripts  # fallback to all

                consequences = set()
                for t in canonical:
                    for csq in t.get("Consequence", "").split("&"):
                        consequences.add(csq)

                if not consequences & CODING_CONSEQUENCES:
                    continue
                self.counters["pass_coding"] += 1

                # ── Step 3: Population frequency filter ──────────────────
                # Use MAX_AF from VEP (maximum across all gnomAD populations)
                max_af = 0.0
                for t in canonical:
                    af_str = t.get("MAX_AF", t.get("gnomADg_AF", "."))
                    max_af = max(max_af, safe_float(af_str, 0.0))

                if max_af > self.args.max_af:
                    continue
                self.counters["pass_frequency"] += 1

                # ── Step 4: Functional impact filter ─────────────────────
                # Keep: all LoF, OR missense with CADD >= threshold
                is_lof = bool(consequences & LOF_CONSEQUENCES)
                is_missense = bool(consequences & MISSENSE_CONSEQUENCES)

                cadd_score = 0.0
                for t in canonical:
                    cadd_score = max(cadd_score, safe_float(t.get("CADD_PHRED", "0"), 0.0))

                # Also accept anything already in ClinVar as pathogenic
                clinvar_sig = ""
                for t in canonical:
                    clinvar_sig = t.get("ClinVar_CLNSIG", "").lower()
                    if clinvar_sig:
                        break
                is_clinvar_path = any(
                    x in clinvar_sig
                    for x in ("pathogenic", "likely_pathogenic")
                )

                if not is_lof and not is_clinvar_path:
                    if not (is_missense and cadd_score >= self.args.min_cadd):
                        continue
                self.counters["pass_functional"] += 1

                # ── Step 5: Phenotype gene filter ────────────────────────
                gene_symbols = set()
                for t in canonical:
                    sym = t.get("SYMBOL", "")
                    if sym:
                        gene_symbols.add(sym)

                in_hpo_list = bool(gene_symbols & HPO_SEIZURE_ID_GENES)

                # Load optional external gene list
                if hasattr(self, "external_genes") and self.external_genes:
                    in_hpo_list = in_hpo_list or bool(
                        gene_symbols & self.external_genes
                    )

                if not in_hpo_list and not is_clinvar_path:
                    continue
                self.counters["pass_phenotype"] += 1

                # ── Passed all filters — collect metadata ─────────────────
                top = canonical[0]
                self.kept_variants.append({
                    "chrom": chrom,
                    "pos": pos,
                    "ref": ref,
                    "alt": alt,
                    "qual": qual,
                    "gene": ",".join(gene_symbols),
                    "consequence": ",".join(sorted(consequences & CODING_CONSEQUENCES)),
                    "hgvsp": top.get("HGVSp", "."),
                    "hgvsc": top.get("HGVSc", "."),
                    "max_af": max_af,
                    "cadd": cadd_score,
                    "sift": top.get("SIFT", "."),
                    "polyphen": top.get("PolyPhen", "."),
                    "clinvar": top.get("ClinVar_CLNSIG", "."),
                    "existing_variation": top.get("Existing_variation", "."),
                    "raw_line": line,
                })

        # Write output
        Path(self.args.output).parent.mkdir(parents=True, exist_ok=True)
        with open(self.args.output, "w") as out:
            for hdr in header_lines:
                out.write(hdr + "\n")
            for v in self.kept_variants:
                out.write(v["raw_line"] + "\n")

        self._print_summary()

    def _print_summary(self):
        c = self.counters
        print(f"\n{'─'*45}")
        print(f"{'Step':<35} {'Variants':>8}")
        print(f"{'─'*45}")
        print(f"{'1. Raw input':<35} {c['total']:>8,}")
        print(f"{'2. Quality filter (QUAL/DP/GQ)':<35} {c['pass_quality']:>8,}")
        print(f"{'3. Coding consequences':<35} {c['pass_coding']:>8,}")
        print(f"{'4. Rare variants (AF < threshold)':<35} {c['pass_frequency']:>8,}")
        print(f"{'5. Functional impact (LoF/CADD)':<35} {c['pass_functional']:>8,}")
        print(f"{'6. Phenotype gene list (HPO)':<35} {c['pass_phenotype']:>8,}")
        print(f"{'─'*45}")
        print(f"\n✓ Filtered VCF written to: {self.args.output}")
        print(f"  Final candidate count: {c['pass_phenotype']}")

        if self.kept_variants:
            print(f"\n  Top candidates:")
            for v in sorted(
                self.kept_variants, key=lambda x: -x["cadd"]
            )[:5]:
                print(
                    f"    {v['chrom']}:{v['pos']} {v['ref']}>{v['alt']}"
                    f"  {v['gene']:12s}  CADD={v['cadd']:.1f}"
                    f"  AF={v['max_af']:.6f}  {v['clinvar']}"
                )


def main():
    parser = argparse.ArgumentParser(
        description="Clinical variant filtering pipeline",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input",    required=True,  help="VEP-annotated VCF (.vcf or .vcf.gz)")
    parser.add_argument("--output",   required=True,  help="Output filtered VCF")
    parser.add_argument("--max_af",   type=float, default=0.01, help="Max gnomAD allele frequency")
    parser.add_argument("--min_qual", type=float, default=30.0,  help="Min QUAL score")
    parser.add_argument("--min_dp",   type=float, default=10.0,  help="Min read depth")
    parser.add_argument("--min_gq",   type=float, default=20.0,  help="Min genotype quality")
    parser.add_argument("--min_cadd", type=float, default=15.0,  help="Min CADD PHRED score for missense")
    parser.add_argument("--hpo_genes", default=None, help="Optional: file with additional gene symbols (one per line)")
    args = parser.parse_args()

    filt = VariantFilter(args)

    if args.hpo_genes and os.path.exists(args.hpo_genes):
        with open(args.hpo_genes) as f:
            filt.external_genes = {line.strip() for line in f if line.strip()}
        print(f"Loaded {len(filt.external_genes)} genes from {args.hpo_genes}")
    else:
        filt.external_genes = set()

    filt.run()


if __name__ == "__main__":
    main()
