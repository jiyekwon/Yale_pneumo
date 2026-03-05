#!/usr/bin/env python3
"""
Pairwise identity check for multi-serotype amplicons.

Compares amplicon sequences ONLY between a assay's intended target serotypes
(parsed from the assay name). Cross-reactive hits — serotypes that produce an
amplicon but are NOT in the assay name — are reported separately.

Output:
  amplicon_audit/pairwise_identity.tsv  — pairwise comparisons (intended targets only)
  amplicon_audit/resolution_summary.tsv — per-assay resolution verdict
  amplicon_audit/cross_reactive_hits.tsv — amplicons on non-intended serotypes
"""

import re
import subprocess
import csv
from pathlib import Path
from difflib import SequenceMatcher
from itertools import combinations
from collections import Counter

PROJECT = Path("/Users/jkwon/Library/CloudStorage/OneDrive-YaleUniversity/Research/Yale_pneumo")
OUT_PW    = PROJECT / "amplicon_audit" / "pairwise_identity.tsv"
OUT_SUM   = PROJECT / "amplicon_audit" / "resolution_summary.tsv"
OUT_XREACT = PROJECT / "amplicon_audit" / "cross_reactive_hits.tsv"

# ─── Parse intended targets from assay name ───────────────────────────────────
def parse_intended(assay_name):
    """
    Extract intended target serotype labels from an assay name and normalize
    to SeroBA zero-padded format (e.g. '9A' → '09A', '1' → '01').

    Rules:
      - Strip parenthetical optional targets e.g. (E)
      - Strip label prefixes like '18F: '
      - Split on '/' and whitespace
      - Tokens that are just letters (A, B, C...) inherit the preceding number
      - Tokens that start with a new number reset the current base
    """
    name = assay_name

    # Strip label prefix like "18F: " (colon-separated heading)
    name = re.sub(r'^[^/]+:\s*', '', name)

    # Remove parenthetical content (optional targets)
    name = re.sub(r'\([^)]*\)', '', name)

    # Split on slashes and whitespace
    tokens = [t.strip() for t in re.split(r'[/\s]+', name) if t.strip()]

    targets = []
    current_num = None

    for token in tokens:
        m = re.match(r'^(\d+)([A-Za-z]*)$', token)
        if m:
            # Token has a number — it's a full serotype label
            current_num = m.group(1)
            targets.append(token.upper())
        elif re.match(r'^[A-Za-z]+$', token) and current_num:
            # Letters only — inherit the preceding number
            targets.append(current_num + token.upper())
        else:
            targets.append(token.upper())

    # Normalize to SeroBA zero-padded format (single-digit numbers → 0X)
    normalized = []
    for t in targets:
        m = re.match(r'^(\d+)([A-Za-z]*)$', t)
        if m and len(m.group(1)) == 1:
            normalized.append('0' + t)
        else:
            normalized.append(t)

    return normalized


# ─── Run seqkit amplicon (BED output, col7 = sequence) ───────────────────────
cmd = [
    "seqkit", "amplicon",
    "-p", str(PROJECT / "primers" / "primer_table_for_seqkit.tsv"),
    "-m", "3", "--bed",
    str(PROJECT / "reference" / "globalpneumoseq_seroba_v2_reference.fasta")
]
result = subprocess.run(cmd, capture_output=True, text=True)
bed_lines = [l for l in result.stdout.splitlines() if l.strip()]
print(f"BED lines: {len(bed_lines)}")

# BED cols: serotype, start, end, primer_name, score, strand, sequence
assay_to_seqs = {}   # {assay_name: {serotype: sequence}}

for line in bed_lines:
    parts = line.split("\t")
    if len(parts) < 7:
        continue
    serotype  = parts[0]
    assay_raw = parts[3]
    seq       = parts[6]

    assay_clean = re.sub(r'_dpo\d+r?\d*$|_rdpo\d+$', '', assay_raw)

    if assay_clean not in assay_to_seqs:
        assay_to_seqs[assay_clean] = {}
    if serotype not in assay_to_seqs[assay_clean] or \
       len(seq) > len(assay_to_seqs[assay_clean][serotype]):
        assay_to_seqs[assay_clean][serotype] = seq

print(f"Assays with sequences: {len(assay_to_seqs)}")

# ─── Pairwise identity (intended targets only) ────────────────────────────────
def pct_identity(a, b):
    return round(SequenceMatcher(None, a, b).ratio() * 100, 1)

pw_rows    = []
sum_rows   = []
xreact_rows = []

for assay, sero_seqs in sorted(assay_to_seqs.items()):
    all_hit      = set(sero_seqs.keys())
    intended_raw = parse_intended(assay)
    # Only keep intended targets that actually produced an amplicon
    intended = [s for s in intended_raw if s in all_hit]
    # Cross-reactive = produced an amplicon but not in intended list
    xreactive = sorted(all_hit - set(intended))

    # Log cross-reactive hits
    for s in xreactive:
        xreact_rows.append({
            "Assay_Set":           assay,
            "Intended_Targets":    ", ".join(intended_raw),
            "Cross_Reactive_Sero": s,
            "Amplicon_bp":         len(sero_seqs[s]),
        })

    n = len(intended)

    if n == 0:
        # All hits are cross-reactive — no intended target found in reference
        sum_rows.append({
            "Assay_Set": assay, "N_Intended": 0,
            "N_Amplified": len(all_hit),
            "Min_Identity": "N/A", "Max_Identity": "N/A",
            "Verdict": "NO_INTENDED_HIT",
            "Identical_Pairs": "", "Distinct_Pairs": "",
            "Cross_Reactive": ", ".join(sorted(xreactive)),
        })
        continue

    if n == 1:
        sum_rows.append({
            "Assay_Set": assay, "N_Intended": 1,
            "N_Amplified": len(all_hit),
            "Min_Identity": "N/A", "Max_Identity": "N/A",
            "Verdict": "SINGLE_TARGET",
            "Identical_Pairs": "", "Distinct_Pairs": "",
            "Cross_Reactive": ", ".join(sorted(xreactive)),
        })
        continue

    pairs = list(combinations(sorted(intended), 2))
    identities = []
    identical_pairs = []
    distinct_pairs  = []

    for s1, s2 in pairs:
        pct = pct_identity(sero_seqs[s1], sero_seqs[s2])
        identities.append(pct)
        pw_rows.append({
            "Assay_Set":    assay,
            "Serotype_A":   s1,
            "Serotype_B":   s2,
            "Len_A":        len(sero_seqs[s1]),
            "Len_B":        len(sero_seqs[s2]),
            "Pct_Identity": pct,
            "Pair_Type":    "INTENDED",
        })
        if pct == 100.0:
            identical_pairs.append(f"{s1}/{s2}")
        else:
            distinct_pairs.append(f"{s1}/{s2}({pct}%)")

    min_id = min(identities)
    max_id = max(identities)

    if min_id == 100.0:
        verdict = "INDISTINGUISHABLE"
    elif min_id >= 99.0:
        verdict = "NEAR-IDENTICAL"
    elif min_id >= 95.0:
        verdict = "LIKELY_RESOLVABLE"
    else:
        verdict = "RESOLVABLE"

    sum_rows.append({
        "Assay_Set":      assay,
        "N_Intended":     n,
        "N_Amplified":    len(all_hit),
        "Min_Identity":   min_id,
        "Max_Identity":   max_id,
        "Verdict":        verdict,
        "Identical_Pairs":  "; ".join(identical_pairs),
        "Distinct_Pairs":   "; ".join(distinct_pairs),
        "Cross_Reactive":   ", ".join(sorted(xreactive)),
    })

# ─── Write outputs ────────────────────────────────────────────────────────────
fieldnames_pw = ["Assay_Set","Serotype_A","Serotype_B","Len_A","Len_B",
                 "Pct_Identity","Pair_Type"]
with open(OUT_PW, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames_pw, delimiter="\t")
    w.writeheader(); w.writerows(pw_rows)

fieldnames_sum = ["Assay_Set","N_Intended","N_Amplified","Min_Identity","Max_Identity",
                  "Verdict","Identical_Pairs","Distinct_Pairs","Cross_Reactive"]
with open(OUT_SUM, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=fieldnames_sum, delimiter="\t")
    w.writeheader(); w.writerows(sum_rows)

with open(OUT_XREACT, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=["Assay_Set","Intended_Targets",
                                       "Cross_Reactive_Sero","Amplicon_bp"],
                       delimiter="\t")
    w.writeheader(); w.writerows(xreact_rows)

# ─── Print summary ────────────────────────────────────────────────────────────
verdict_counts = Counter(r["Verdict"] for r in sum_rows)

print(f"\n=== RESOLUTION SUMMARY (intended targets only) ===")
for v, c in sorted(verdict_counts.items()):
    print(f"  {v:22s}: {c}")

print(f"\nINDISTINGUISHABLE:")
for r in sum_rows:
    if r["Verdict"] == "INDISTINGUISHABLE":
        print(f"  {r['Assay_Set']:30s}  {r['Identical_Pairs']}")

print(f"\nNEAR-IDENTICAL (>99%):")
for r in sum_rows:
    if r["Verdict"] == "NEAR-IDENTICAL":
        print(f"  {r['Assay_Set']:30s}  {r['Distinct_Pairs']}")

print(f"\nLIKELY_RESOLVABLE (95-99%):")
for r in sum_rows:
    if r["Verdict"] == "LIKELY_RESOLVABLE":
        print(f"  {r['Assay_Set']:30s}  min={r['Min_Identity']}%  n={r['N_Intended']}")

print(f"\nRESOLVABLE (<95% min identity):")
for r in sum_rows:
    if r["Verdict"] == "RESOLVABLE":
        print(f"  {r['Assay_Set']:30s}  min={r['Min_Identity']}%  n={r['N_Intended']}")

print(f"\nCross-reactive hits: {len(xreact_rows)} total")
xr_by_assay = Counter(r["Assay_Set"] for r in xreact_rows)
for assay, count in xr_by_assay.most_common(10):
    print(f"  {assay:30s}  {count} non-intended serotypes amplified")

print(f"\nOutputs: {OUT_SUM.name}, {OUT_PW.name}, {OUT_XREACT.name}")
