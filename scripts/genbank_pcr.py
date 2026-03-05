#!/usr/bin/env python3
"""
In silico PCR against GenBank (.gb) files in test_gb/.

Used for assays that produce NO HIT against the SeroBA CPS reference —
typically because the primer target gene lies outside the CPS locus.

For each .gb file:
  1. Extracts the full DNA sequence as FASTA
  2. Runs seqkit amplicon against it using the full primer table
  3. Reports all hits (no size filter — not doing amplicon sequencing)

Output:
  amplicon_audit/genbank_amplicon_results.tsv
"""

import re
import csv
import subprocess
import tempfile
from pathlib import Path

PROJECT   = Path("/Users/jkwon/Library/CloudStorage/OneDrive-YaleUniversity/Research/Yale_pneumo")
GB_DIR    = PROJECT / "test_gb"
PRIMER_TSV = PROJECT / "primers" / "primer_table_for_seqkit.tsv"
OUT_FILE  = PROJECT / "amplicon_audit" / "genbank_amplicon_results.tsv"

# ─── GenBank → FASTA extraction ───────────────────────────────────────────────
def gb_to_fasta(gb_path: Path) -> tuple[str, str, str]:
    """
    Parse a GenBank file and return (accession, definition, sequence).
    """
    text = gb_path.read_text(encoding="utf-8")

    accession  = re.search(r'^ACCESSION\s+(\S+)', text, re.MULTILINE)
    definition = re.search(r'^DEFINITION\s+(.+?)(?=\n\S)', text, re.DOTALL)
    origin     = re.search(r'^ORIGIN(.*?)^//', text, re.DOTALL | re.MULTILINE)

    acc = accession.group(1) if accession else gb_path.stem
    dfn = definition.group(1).replace('\n', ' ').strip() if definition else ""
    seq = ""
    if origin:
        seq = re.sub(r'[\d\s]', '', origin.group(1)).upper()

    return acc, dfn, seq


# ─── Run seqkit against each .gb file ─────────────────────────────────────────
gb_files = sorted(GB_DIR.glob("*.gb"))
if not gb_files:
    print(f"No .gb files found in {GB_DIR}")
    raise SystemExit(0)

print(f"Found {len(gb_files)} GenBank file(s) in {GB_DIR.name}/")

rows = []

for gb_path in gb_files:
    acc, dfn, seq = gb_to_fasta(gb_path)
    if not seq:
        print(f"  [WARN] No sequence extracted from {gb_path.name}")
        continue

    # Write temp FASTA
    with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as tmp:
        tmp.write(f">{acc} {dfn}\n{seq}\n")
        tmp_path = tmp.name

    # Run seqkit amplicon
    cmd = [
        "seqkit", "amplicon",
        "-p", str(PRIMER_TSV),
        "-m", "3", "--bed",
        tmp_path
    ]
    result = subprocess.run(cmd, capture_output=True, text=True)
    Path(tmp_path).unlink()

    hits = [l for l in result.stdout.splitlines() if l.strip()]

    if not hits:
        print(f"  {gb_path.name:30s}  → no amplicons")
        rows.append({
            "GB_File":     gb_path.name,
            "Accession":   acc,
            "Definition":  dfn[:80],
            "Assay":       "—",
            "Amplicon_bp": "—",
            "Sequence":    "—",
            "Note":        "NO_HIT",
        })
        continue

    for line in hits:
        parts = line.split("\t")
        if len(parts) < 7:
            continue
        start, end, assay, seq_out = int(parts[1]), int(parts[2]), parts[3], parts[6]
        bp = end - start
        print(f"  {gb_path.name:30s}  assay={assay:20s}  {bp} bp")
        rows.append({
            "GB_File":     gb_path.name,
            "Accession":   acc,
            "Definition":  dfn[:80],
            "Assay":       assay,
            "Amplicon_bp": bp,
            "Sequence":    seq_out,
            "Note":        "HIT",
        })

# ─── Write output ──────────────────────────────────────────────────────────────
with open(OUT_FILE, "w", newline="", encoding="utf-8-sig") as f:
    w = csv.DictWriter(f, fieldnames=[
        "GB_File", "Accession", "Definition",
        "Assay", "Amplicon_bp", "Sequence", "Note"
    ], delimiter="\t")
    w.writeheader()
    w.writerows(rows)

print(f"\nOutput: {OUT_FILE}")
print(f"Total hits: {sum(1 for r in rows if r['Note'] == 'HIT')}")
print(f"No-hit files: {sum(1 for r in rows if r['Note'] == 'NO_HIT')}")
