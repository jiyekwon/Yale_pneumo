#!/usr/bin/env python3
"""
Generate final_resolution_summary.tsv.

Combines:
  1. amplicon_audit/resolution_summary.tsv      — per-assay pairwise verdict
  2. amplicon_audit/combinatorial_resolution.tsv — per-pair rescue verdict
  3. amplicon_audit/genbank_amplicon_results.tsv — GenBank supplement hits
  4. primers/primer_table_for_seqkit.tsv         — to identify no-hit assays

For INDISTINGUISHABLE assays: checks whether ALL identical pairs are rescued
by combination (RESOLVABLE_BY_COMBINATION) → promotes to RESCUED (combinatorial).

For SeroBA no-hit assays: checks GenBank supplement for confirmation.
  - Valid hit = amplicon within 50–500 bp range (typical for qPCR)
  - Giant amplicons (thousands of bp) = non-specific, ignored.

Output:
  amplicon_audit/final_resolution_summary.tsv
"""

import csv
import re
from pathlib import Path
from collections import defaultdict

PROJECT  = Path("/Users/jkwon/Library/CloudStorage/OneDrive-YaleUniversity/Research/Yale_pneumo")
IN_SUM   = PROJECT / "amplicon_audit" / "resolution_summary.tsv"
IN_COMB  = PROJECT / "amplicon_audit" / "combinatorial_resolution.tsv"
IN_GB    = PROJECT / "amplicon_audit" / "genbank_amplicon_results.tsv"
PRIMER_TSV = PROJECT / "primers" / "primer_table_for_seqkit.tsv"
OUT_FILE = PROJECT / "amplicon_audit" / "final_resolution_summary.tsv"

# Max amplicon bp to consider a GenBank hit "real" (qPCR context)
GB_MAX_BP = 500

# ─── Helpers (mirrors parse_intended from pairwise_identity.py) ──────────────
def parse_intended(assay_name):
    name = assay_name
    name = re.sub(r'^[^/]+:\s*', '', name)
    name = re.sub(r'\([^)]*\)', '', name)
    tokens = [t.strip() for t in re.split(r'[/\s]+', name) if t.strip()]
    targets = []
    current_num = None
    for token in tokens:
        m = re.match(r'^(\d+)([A-Za-z]*)$', token)
        if m:
            current_num = m.group(1)
            targets.append(token.upper())
        elif re.match(r'^[A-Za-z]+$', token) and current_num:
            targets.append(current_num + token.upper())
        else:
            targets.append(token.upper())
    normalized = []
    for t in targets:
        m = re.match(r'^(\d+)([A-Za-z]*)$', t)
        if m and len(m.group(1)) == 1:
            normalized.append('0' + t)
        else:
            normalized.append(t)
    return normalized

def strip_dpo(name):
    return re.sub(r'_dpo\d+r?\d*$|_rdpo\d+$', '', name)

# ─── Load combinatorial resolution ───────────────────────────────────────────
# Build set of pairs that ARE resolvable by combination
comb_resolvable = set()   # {(sA, sB)}
comb_rescuer    = {}      # {(sA, sB): distinguishing_assays_string}

with open(IN_COMB, encoding="utf-8-sig") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["Verdict"] == "RESOLVABLE_BY_COMBINATION":
            pair = (row["Serotype_A"], row["Serotype_B"])
            comb_resolvable.add(pair)
            comb_resolvable.add((pair[1], pair[0]))  # both orderings
            comb_rescuer[pair] = row["Distinguishing_Assays"]
            comb_rescuer[(pair[1], pair[0])] = row["Distinguishing_Assays"]

# ─── Load GenBank supplement — valid hits only ───────────────────────────────
gb_valid_assays = defaultdict(set)  # {assay_base: set of GB filenames with valid hits}

with open(IN_GB, encoding="utf-8-sig") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        if row["Note"] != "HIT":
            continue
        try:
            bp = int(row["Amplicon_bp"])
        except ValueError:
            continue
        if bp < 0 or bp > GB_MAX_BP:
            continue
        assay_base = strip_dpo(row["Assay"])
        gb_valid_assays[assay_base].add(row["GB_File"])

# ─── Get all assay names from seqkit primer table ────────────────────────────
all_assays_in_panel = set()
with open(PRIMER_TSV) as f:
    for line in f:
        parts = line.strip().split("\t")
        if parts:
            all_assays_in_panel.add(strip_dpo(parts[0]))

# ─── Load per-assay resolution summary ───────────────────────────────────────
sum_rows = []
assays_with_hits = set()

with open(IN_SUM, encoding="utf-8-sig") as f:
    for row in csv.DictReader(f, delimiter="\t"):
        assays_with_hits.add(row["Assay_Set"])
        sum_rows.append(row)

# Assays in panel but NOT in SeroBA results
no_hit_assays = all_assays_in_panel - assays_with_hits

# ─── Build final summary rows ─────────────────────────────────────────────────
out_rows = []

for row in sum_rows:
    assay   = row["Assay_Set"]
    verdict = row["Verdict"]
    intended = parse_intended(assay)
    n_intended = len(intended)

    intended_str = ", ".join(intended) if intended else assay

    # --- Determine Final_Status ---
    if verdict == "SINGLE_TARGET":
        final_status = "UNIQUE"
        note = ""

    elif verdict in ("RESOLVABLE", "LIKELY_RESOLVABLE", "NEAR-IDENTICAL"):
        final_status = verdict
        if verdict == "RESOLVABLE":
            note = "Clear sequence differences"
        elif verdict == "LIKELY_RESOLVABLE":
            note = "SNPs present; needs good read depth"
        else:
            note = ">99% identity; SNP-level calling required"

    elif verdict == "NO_INTENDED_HIT":
        final_status = "NO_INTENDED_HIT"
        xr = row.get("Cross_Reactive", "")
        note = f"Amplifies {xr} — assay name does not match new DB; needs remapping" if xr else ""

    elif verdict == "INDISTINGUISHABLE":
        # Check if all identical pairs are rescued combinatorially
        identical_pairs_str = row.get("Identical_Pairs", "")
        identical_pairs = [p.strip() for p in identical_pairs_str.split(";") if p.strip()]

        # Parse pairs: "09A/09V" → ("09A", "09V")
        pairs = []
        for p in identical_pairs:
            parts = p.split("/")
            if len(parts) == 2:
                pairs.append((parts[0].strip(), parts[1].strip()))

        if not pairs:
            final_status = "INDISTINGUISHABLE"
            note = "Serogroup-level ceiling; no rescue assay"
        else:
            rescued_pairs = [p for p in pairs if p in comb_resolvable]
            if len(rescued_pairs) == len(pairs):
                # All pairs rescued
                rescuers = set()
                for p in pairs:
                    da = comb_rescuer.get(p, "")
                    for entry in da.split(";"):
                        entry = entry.strip()
                        if entry:
                            # Extract assay name before "("
                            assay_name = re.sub(r'\(.*\)', '', entry).strip()
                            rescuers.add(assay_name)
                final_status = "RESCUED (combinatorial)"
                note = "Rescued by: " + "; ".join(sorted(rescuers))
            else:
                final_status = "INDISTINGUISHABLE"
                note = "Serogroup-level ceiling; no rescue assay"

    else:
        final_status = verdict
        note = ""

    out_rows.append({
        "Assay_Set":         assay,
        "Intended_Serotypes": intended_str,
        "N_Serotypes":        n_intended,
        "Assay_Verdict":      verdict,
        "Final_Status":       final_status,
        "Note":               note,
    })

# ─── Add SeroBA no-hit assays ─────────────────────────────────────────────────
for assay in sorted(no_hit_assays):
    intended = parse_intended(assay)
    intended_str = ", ".join(intended) if intended else assay

    if assay in gb_valid_assays:
        gb_files = ", ".join(sorted(gb_valid_assays[assay]))
        final_status = "CONFIRMED (GenBank)"
        note = f"No hit in SeroBA CPS locus; confirmed via GenBank: {gb_files}"
    else:
        final_status = "NO_HIT"
        note = "No amplicon in SeroBA or GenBank; target gene likely absent from reference"

    out_rows.append({
        "Assay_Set":          assay,
        "Intended_Serotypes": intended_str,
        "N_Serotypes":        len(intended),
        "Assay_Verdict":      "NO_HIT_IN_SEROBA",
        "Final_Status":       final_status,
        "Note":               note,
    })

# ─── Sort and write ───────────────────────────────────────────────────────────
out_rows.sort(key=lambda r: r["Assay_Set"])

with open(OUT_FILE, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=[
        "Assay_Set", "Intended_Serotypes", "N_Serotypes",
        "Assay_Verdict", "Final_Status", "Note"
    ], delimiter="\t")
    w.writeheader()
    w.writerows(out_rows)

# ─── Print summary ────────────────────────────────────────────────────────────
from collections import Counter
counts = Counter(r["Final_Status"] for r in out_rows)
print(f"\n=== FINAL RESOLUTION SUMMARY ({len(out_rows)} assays) ===")
for status, count in sorted(counts.items()):
    print(f"  {status:35s}: {count}")

print(f"\nOutput: {OUT_FILE}")
