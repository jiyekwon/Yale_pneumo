#!/usr/bin/env python3
"""
Combinatorial resolution analysis.

For each pair of serotypes that appear indistinguishable by their primary assay,
check ALL other assays that produce amplicons for both serotypes. If any other
assay gives different sequences for the two, the combination resolves them.

This mirrors the Downs panel design: no single assay needs to distinguish
everything — the pattern of results across the full panel does the typing.

Output:
  amplicon_audit/combinatorial_resolution.tsv — per-serotype-pair verdict
"""

import re
import subprocess
import csv
from pathlib import Path
from difflib import SequenceMatcher
from itertools import combinations
from collections import defaultdict

PROJECT  = Path("/Users/jkwon/Library/CloudStorage/OneDrive-YaleUniversity/Research/Yale_pneumo")
OUT_FILE = PROJECT / "amplicon_audit" / "combinatorial_resolution.tsv"

# ─── Re-run seqkit to get all amplicons (BED with sequence) ──────────────────
cmd = [
    "seqkit", "amplicon",
    "-p", str(PROJECT / "primers" / "primer_table_for_seqkit.tsv"),
    "-m", "3", "--bed",
    str(PROJECT / "reference" / "globalpneumoseq_seroba_v2_reference.fasta")
]
result = subprocess.run(cmd, capture_output=True, text=True)
bed_lines = [l for l in result.stdout.splitlines() if l.strip()]

def pct_identity(a, b):
    return round(SequenceMatcher(None, a, b).ratio() * 100, 1)

def strip_dpo(name):
    return re.sub(r'_dpo\d+r?\d*$|_rdpo\d+$', '', name)

# ─── Build full map: {serotype: {assay: sequence}} ───────────────────────────
sero_to_assay_seq = defaultdict(dict)   # {serotype: {assay: seq}}

for line in bed_lines:
    parts = line.split("\t")
    if len(parts) < 7:
        continue
    serotype = parts[0]
    assay    = strip_dpo(parts[3])
    seq      = parts[6]
    # Keep longest match per assay/serotype combination
    if assay not in sero_to_assay_seq[serotype] or \
       len(seq) > len(sero_to_assay_seq[serotype][assay]):
        sero_to_assay_seq[serotype][assay] = seq

# ─── For each serotype pair: find all assays that amplify both ────────────────
all_serotypes = sorted(sero_to_assay_seq.keys())
print(f"Serotypes with amplicons: {len(all_serotypes)}")

rows = []
pair_count = 0

for s1, s2 in combinations(all_serotypes, 2):
    assays_s1 = set(sero_to_assay_seq[s1].keys())
    assays_s2 = set(sero_to_assay_seq[s2].keys())
    shared_assays = sorted(assays_s1 & assays_s2)

    if not shared_assays:
        continue   # No shared assay — cannot compare at all

    pair_count += 1
    identical_assays  = []
    distinct_assays   = []
    best_diff = 100.0  # track the most-different assay pair

    for assay in shared_assays:
        pct = pct_identity(sero_to_assay_seq[s1][assay],
                           sero_to_assay_seq[s2][assay])
        if pct == 100.0:
            identical_assays.append(assay)
        else:
            distinct_assays.append(f"{assay}({pct}%)")
            if pct < best_diff:
                best_diff = pct

    if distinct_assays:
        verdict = "RESOLVABLE_BY_COMBINATION"
        min_id  = best_diff
    elif identical_assays:
        verdict = "TRULY_INDISTINGUISHABLE"
        min_id  = 100.0
    else:
        verdict = "NO_SHARED_ASSAY"
        min_id  = None

    rows.append({
        "Serotype_A":        s1,
        "Serotype_B":        s2,
        "N_Shared_Assays":   len(shared_assays),
        "Verdict":           verdict,
        "Min_Identity":      min_id,
        "Distinguishing_Assays": "; ".join(distinct_assays),
        "Identical_Assays":      "; ".join(identical_assays),
    })

# ─── Write output ─────────────────────────────────────────────────────────────
rows.sort(key=lambda r: (r["Verdict"], r["Serotype_A"], r["Serotype_B"]))

with open(OUT_FILE, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=[
        "Serotype_A","Serotype_B","N_Shared_Assays","Verdict",
        "Min_Identity","Distinguishing_Assays","Identical_Assays"
    ], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

# ─── Summary ──────────────────────────────────────────────────────────────────
from collections import Counter
counts = Counter(r["Verdict"] for r in rows)

print(f"\n=== COMBINATORIAL RESOLUTION (across all shared assays) ===")
print(f"  Total serotype pairs with ≥1 shared assay: {pair_count}")
for v, c in sorted(counts.items()):
    print(f"  {v:35s}: {c}")

truly_indist = [r for r in rows if r["Verdict"] == "TRULY_INDISTINGUISHABLE"]
print(f"\nTRULY INDISTINGUISHABLE pairs ({len(truly_indist)}):")
print(f"  (identical amplicons across ALL shared assays — cannot be resolved by this panel)")
for r in truly_indist:
    print(f"  {r['Serotype_A']:6s} vs {r['Serotype_B']:6s}  "
          f"({r['N_Shared_Assays']} shared assays: {r['Identical_Assays']})")

# Highlight previously-indistinguishable pairs that ARE rescued by combination
prev_indist_pairs = [
    ("10C","10F"),("11A","11D"),("11B","11C"),("15B","15C"),("19B","19F"),
    ("24B","24F"),("25A","25F"),("25A","38"),("25F","38"),("28A","28F"),
    ("32A","32F"),("33A","33F"),("47A","47F"),("06C","06D"),("07A","07F"),
    ("09A","09V"),("09L","09N"),
]
# normalize to match reference format
prev_set = {(a.zfill(2) if a[0].isdigit() and len(a)<=2 else a,
             b.zfill(2) if b[0].isdigit() and len(b)<=2 else b)
            for a, b in prev_indist_pairs}

rescued = []
for r in rows:
    pair = (r["Serotype_A"], r["Serotype_B"])
    pair_rev = (r["Serotype_B"], r["Serotype_A"])
    if (pair in prev_set or pair_rev in prev_set) and \
       r["Verdict"] == "RESOLVABLE_BY_COMBINATION":
        rescued.append(r)

print(f"\nPreviously INDISTINGUISHABLE pairs rescued by combination ({len(rescued)}):")
for r in rescued:
    print(f"  {r['Serotype_A']:6s} vs {r['Serotype_B']:6s}  "
          f"→ distinguished by: {r['Distinguishing_Assays']}")

print(f"\nOutput: {OUT_FILE}")
