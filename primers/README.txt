PRIMERS DIRECTORY — DATA PROVENANCE & CHANGE LOG
=================================================
Project: PCR-Seq Pneumococcal Serotyping Pipeline
Last updated: 2026-02-24


FILES IN THIS DIRECTORY
------------------------
1. downs2023_pneumo_primers.csv
   - 75 pneumococcal serotype assay sets
   - Pool assignments (A/B/C) from updated Lancet 2023 panel

2. downs2023_other_bacteria_primers.csv
   - 24 non-pneumococcal bacterial species assay sets from Downs 2023
   - Includes S. pneumoniae reference genes (LytA, Xisco, Ply, PiaB),
     H. influenzae subtypes, and other respiratory pathogens

3. lancet2023_additional_assays.csv
   - 14 viral and additional assay sets from the updated Lancet 2023 panel
   - These are Mallery's additional primers not in the original Downs 2023 paper


SOURCES
--------
Primary source (Downs 2023):
  Downs SL, Madhi SA, Nunes MC, Olwagen CP.
  "Optimization of a high-throughput nanofluidic real-time PCR to detect and
  quantify 15 bacterial species and 92 Streptococcus pneumoniae serotypes."
  Scientific Reports, 13, 4588 (2023).
  DOI: https://doi.org/10.1038/s41598-023-31820-4
  Raw data: LitReview/41598_2023_31820_MOESM1_ESM (1).docx
  Tables used: Supplementary Table 1 (oligo sequences), Table 5 (pool assignments)

Updated source (Lancet 2023):
  Additional assays from Lancet Microbe 2023.
  DOI: https://doi.org/10.1016/S2666-5247(23)00260-4
  Raw data: LitReview/addtl_serotype_supp.pdf
  Tables used: Supplementary Table 2 (full target list), Table 3 (new oligo sequences),
               Table 4 (updated pool assignments)

Personal communication:
  Provided clarifications on discrepancies in the Downs 2023 supplementary data
  and confirmed the Lancet 2023 additional primers.


HOW THE CSVs WERE BUILT
------------------------
Step 1 — Extracted Supplementary Table 1 from the Downs 2023 .docx file
         using python-docx. All 97 rows were parsed.

Step 2 — Separated into pneumococcal serotype assays (73) vs. other bacterial
         species assays (24) based on assay set name.

Step 3 — Cleaned encoding: the original symbols dagger and double-dagger were
         stripped from the Assay_Set column and moved to a new Oligo_Type column:
           Standard         = original dagger symbol (standard TaqMan assay)
           Modified (DPO/LNA) = original double-dagger symbol (modified oligonucleotide)
         Files saved with UTF-8 BOM so Excel opens them without garbling characters.

Step 4 — Added Pool column (A/B/C) from Lancet 2023 Supplementary Table 4,
         which supersedes the pool assignments in Downs 2023.

Step 5 — Added 19C and 19F Atypical from Lancet 2023 Supplementary Table 3
         (these were not in the original Downs 2023 paper).

Step 6 — Added Notes column capturing clarifications

Step 7 — Created lancet2023_additional_assays.csv for the 14 viral/additional
         assay sets from Lancet 2023 Supplementary Table 3.


CLARIFICATIONS FROM EMAIL
-------------------------------------
1. ASSAY COUNT
   The paper states 96 assay sets but Supplementary Table 1 lists 97.
   Reason: Pneumocystis jirovecii was included in the table as it was described
   in the manuscript, but was replaced in the actual panel by an additional
   pneumococcal reference gene. The journal will be notified of the correction.
   Action: Pneumocystis jirovecii assay is flagged DO NOT USE in the CSV.

2. SEROGROUP 18 NAMING DISCREPANCY
   18A/B/C/F (oligo table) = 18A/B/C (pool table) — these are the same assay.
   The oligo is named 18A/B/C/F but sometimes fails to detect 18F, so the
   pool table uses the more accurate name 18A/B/C.
   Action: Flagged in Notes column.

3. ASSAY 18C/F
   18C/F was added to Pool C in place of 23A (23A was removed from the panel).
   18C/F was designed after the g-blocks were ordered, so it has no g-block
   positive control — bacterial cultures must be used to validate this assay.
   Action: Flagged in Notes column. If ordering new g-blocks, add 18C/F sequence.

4. SEROTYPE 23A
   23A was removed from the active panel. It can be inferred algorithmically
   from other assay results. It is a useful addition if space allows.
   Action: 23A/B/F assay retained in CSV but flagged in Notes column.

5. UPDATED PANEL SIZE
   The current panel (Lancet 2023) tests 92 pneumococcal serotypes plus
   26 other bacterial and viral targets. This exceeds 96 assays and requires
   the Biomark HD platform for multiplexing.
   The Downs 2023 panel was designed for the Standard BioTools Fluidigm system.

6. G-BLOCKS
   g-Blocks are synthetic double-stranded DNA fragments (IDT) used as positive
   controls. The Downs 2023 paper pooled multiple target sequences into single
   g-block constructs (A1, A2, etc.) for efficiency. Assay 18C/F is not
   represented in any existing g-block and requires bacterial culture controls.


KNOWN ISSUES & DECISIONS
--------------------------
- Some assay names use inconsistent capitalization of gene targets in the source
  docx (e.g., "Wzx" vs "wzx" vs "WZX"). These have been left as-is from the
  source to avoid introducing errors; normalize before use in bioinformatics.

- Some Modified (DPO/LNA) assays contain multiple forward primer sequences
  separated by spaces in the source document (e.g., 6A/C, 18C/F). These need
  to be parsed carefully for in silico PCR — each space-separated sequence
  is a distinct primer in a DPO design.

- hRV (Human Rhinovirus) in lancet2023_additional_assays.csv uses 3 forward
  primers simultaneously with degenerate bases throughout. This requires
  special handling in the amplicon audit step.

- Gene targets labeled "cps" for 19C and 19F Atypical are placeholders —
  the specific CPS gene targeted should be confirmed from the Lancet 2023
  full paper.

- downs2023_other_bacteria_primers.csv is retained for completeness but is
  not the focus of the PCR-Seq serotyping pipeline.


NEXT STEPS (as of 2026-02-24)
-------------------------------
1. Amplicon size audit — run in silico PCR for all 75 pneumo assays against
   NCBI pneumococcal reference genomes. Flag amplicons outside 150-400 bp.
2. Cross-reactivity screen — check all primer pairs against each other.
3. Confirm gene targets for 19C and 19F Atypical from Lancet 2023 full text.
4. Normalize gene target capitalization across CSVs before bioinformatics use.
