#!/usr/bin/env python3
"""
In silico PCR for pneumo serotyping primer panel.
Uses seqkit amplicon against SeroBA CPS reference.

For DPO primers with multiple space-separated forward sequences,
each sub-sequence is tested independently.

Output:
  - primers/primer_table_for_seqkit.tsv   — input for seqkit
  - amplicon_audit/amplicon_results.tsv   — parsed results table
"""

import csv
import subprocess
import os
import re
from pathlib import Path

# Paths
PROJECT = Path("/Users/jkwon/Library/CloudStorage/OneDrive-YaleUniversity/Research/Yale_pneumo")
PRIMERS_CSV = PROJECT / "primers" / "downs2023_pneumo_primers.csv"
REFERENCE = PROJECT / "reference" / "globalpneumoseq_seroba_v2_reference.fasta"
PRIMER_TSV = PROJECT / "primers" / "primer_table_for_seqkit.tsv"
OUT_DIR = PROJECT / "amplicon_audit"
AMPLICON_FASTA = OUT_DIR / "amplicon_raw.fasta"
RESULTS_TSV = OUT_DIR / "amplicon_results.tsv"

OUT_DIR.mkdir(exist_ok=True)

# ─── Step 1: Build primer table for seqkit ───────────────────────────────────
print("Building primer table...")

entries = []  # list of (entry_name, fwd, rev, assay_set, pool, oligo_type)

with open(PRIMERS_CSV, encoding="utf-8-sig") as f:
    for row in csv.DictReader(f):
        assay = row["Assay_Set"].strip()
        rev   = row["Reverse_Primer"].strip()
        pool  = row["Pool"].strip()
        oligo = row["Oligo_Type"].strip()

        fwd_seqs = row["Forward_Primer"].strip().split()  # split on whitespace
        rev_seqs = row["Reverse_Primer"].strip().split()  # split on whitespace (DPO rev)

        # Build all fwd × rev sub-sequence combinations
        for fi, fwd in enumerate(fwd_seqs, 1):
            for ri, rev_seq in enumerate(rev_seqs, 1):
                if len(fwd_seqs) == 1 and len(rev_seqs) == 1:
                    entry_name = assay
                elif len(fwd_seqs) > 1 and len(rev_seqs) == 1:
                    entry_name = f"{assay}_dpo{fi}"
                elif len(fwd_seqs) == 1 and len(rev_seqs) > 1:
                    entry_name = f"{assay}_rdpo{ri}"
                else:
                    entry_name = f"{assay}_dpo{fi}r{ri}"
                entries.append((entry_name, fwd, rev_seq, assay, pool, oligo))

# Write seqkit primer file (3 columns: name, fwd, rev — tab separated)
with open(PRIMER_TSV, "w") as f:
    for entry_name, fwd, rev, *_ in entries:
        f.write(f"{entry_name}\t{fwd}\t{rev}\n")

print(f"  Wrote {len(entries)} primer entries to {PRIMER_TSV}")

# ─── Step 2: Run seqkit amplicon ─────────────────────────────────────────────
print("\nRunning seqkit amplicon...")

cmd = [
    "seqkit", "amplicon",
    "-p", str(PRIMER_TSV),
    "-m", "3",           # allow up to 3 mismatches
    "--bed",             # output BED format (gives coordinates)
    str(REFERENCE)
]

result = subprocess.run(cmd, capture_output=True, text=True)

if result.returncode != 0:
    print("seqkit stderr:", result.stderr[:500])
    raise RuntimeError("seqkit amplicon failed")

bed_output = result.stdout
print(f"  seqkit done. Lines of output: {len(bed_output.splitlines())}")

# ─── Step 3: Also get amplicon sequences ─────────────────────────────────────
print("\nRunning seqkit amplicon (sequence output)...")

cmd_seq = [
    "seqkit", "amplicon",
    "-p", str(PRIMER_TSV),
    "-m", "3",
    str(REFERENCE)
]

result_seq = subprocess.run(cmd_seq, capture_output=True, text=True)
with open(AMPLICON_FASTA, "w") as f:
    f.write(result_seq.stdout)
print(f"  Amplicon FASTA saved to {AMPLICON_FASTA}")

# ─── Step 4: Parse BED output ────────────────────────────────────────────────
print("\nParsing results...")

# BED columns: chrom, start, end, name, score, strand
# chrom = serotype name from reference
# name  = primer entry name (assay or assay_dpoN)

# Build lookup: entry_name -> (assay_set, pool, oligo_type)
meta = {}
for entry_name, fwd, rev, assay, pool, oligo in entries:
    meta[entry_name] = {"assay": assay, "pool": pool, "oligo": oligo, "fwd": fwd, "rev": rev}

rows_out = []
for line in bed_output.splitlines():
    if not line.strip():
        continue
    parts = line.split("\t")
    if len(parts) < 6:
        continue
    serotype, start, end, primer_name = parts[0], int(parts[1]), int(parts[2]), parts[3]
    amplicon_bp = end - start
    m = meta.get(primer_name, {})

    # Flag size
    if amplicon_bp < 150:
        flag = "TOO_SHORT"
    elif amplicon_bp > 600:
        flag = "TOO_LONG"
    else:
        flag = "OK"

    rows_out.append({
        "Assay_Set":    m.get("assay", primer_name),
        "Primer_Entry": primer_name,
        "Pool":         m.get("pool", ""),
        "Oligo_Type":   m.get("oligo", ""),
        "Reference_Serotype": serotype,
        "Amplicon_Start": start,
        "Amplicon_End":   end,
        "Amplicon_bp":    amplicon_bp,
        "Size_Flag":      flag,
    })

# Write results
rows_out.sort(key=lambda r: (r["Assay_Set"], r["Reference_Serotype"]))

with open(RESULTS_TSV, "w", newline="", encoding="utf-8-sig") as f:
    writer = csv.DictWriter(f, fieldnames=[
        "Assay_Set","Primer_Entry","Pool","Oligo_Type",
        "Reference_Serotype","Amplicon_Start","Amplicon_End","Amplicon_bp","Size_Flag"
    ], delimiter="\t")
    writer.writeheader()
    writer.writerows(rows_out)

# ─── Step 5: Summary ─────────────────────────────────────────────────────────
total = len(rows_out)
ok    = sum(1 for r in rows_out if r["Size_Flag"] == "OK")
short = sum(1 for r in rows_out if r["Size_Flag"] == "TOO_SHORT")
long_ = sum(1 for r in rows_out if r["Size_Flag"] == "TOO_LONG")

# Assays with no hit at all
hit_assays = {r["Assay_Set"] for r in rows_out}
all_assays = {m["assay"] for m in meta.values()}
no_hit = all_assays - hit_assays

print(f"\n=== RESULTS ===")
print(f"  Total amplicons found : {total}")
print(f"  OK (150-600 bp)       : {ok}")
print(f"  Too short (<150 bp)   : {short}")
print(f"  Too long (>600 bp)    : {long_}")
print(f"  Assays with NO hit    : {len(no_hit)}")
if no_hit:
    for a in sorted(no_hit):
        print(f"    - {a}")
print(f"\nResults saved to: {RESULTS_TSV}")
print(f"Amplicon FASTA : {AMPLICON_FASTA}")
