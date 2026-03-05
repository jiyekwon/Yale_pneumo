# Yale Pneumo PCR-Seq Serotyping Pipeline

In silico validation of a multiplex PCR primer panel for sequencing-based pneumococcal
serotyping. Adapts the Downs et al. (2023) qPCR assay set — originally designed for
nanofluidic real-time PCR — to an amplicon sequencing workflow.

## Background

*Streptococcus pneumoniae* (pneumococcus) serotyping. The Downs 2023 panel covers 92 pneumococcal
serotypes using 73 multiplexed primer sets in three PCR pools (A, B, C). This project
assesses whether the same primer panel can drive an amplicon sequencing workflow —
replacing probe-based detection with sequence-based identification.

## Primer Panel

| Source | Assays | Coverage |
|--------|--------|----------|
| Downs et al. (2023), *Sci Rep* | 73 assay sets | 92 pneumococcal serotypes |
| Olwagen et al. (2023), *Lancet Microbe* | 2 additional | 19C, 19F Atypical |
| **Total pneumococcal** | **75 assay sets** | **94 serotypes** |

Pool assignments (A/B/C) follow Lancet 2023 Supplementary Table 4, which supersedes
Downs 2023. Probes are stripped; only forward/reverse primers are used.

## Repository Structure

| Folder/File | Contents |
|-------------|----------|
| `primers/` | Primer sequences and metadata (CSVs + provenance README) |
| `reference/` | SeroBA v2 CPS reference sequences and metadata |
| `amplicon_audit/` | All in silico PCR outputs and resolution analysis results |
| `scripts/` | Python pipeline scripts (run in order; see below) |
| `test_gb/` | GenBank `.gb` files for no-hit assay validation |
| `pipeline_summary.txt` | Detailed 4-phase pipeline plan |

## Results Summary

**75 assay sets covering 94 pneumococcal serotypes** (against SeroBA v2 CPS reference):
as of March 5, 2026. 

| Final Status | Assays |
|--------------|-------:|
| UNIQUE | 35 |
| RESCUED (combinatorial) | 13 |
| RESOLVABLE | 7 |
| LIKELY_RESOLVABLE | 10 |
| NEAR-IDENTICAL | 2 |
| INDISTINGUISHABLE | 3 |
| NO_INTENDED_HIT | 3 |
| CONFIRMED (GenBank) | 1 |
| NO_HIT | 1 |
| **Total** | **75** |

Full output: `amplicon_audit/final_resolution_summary.tsv`

**Truly indistinguishable pairs** (identical amplicons across all shared assays):
- 11A / 11D
- 32A / 32F
- 33A / 33F / 37

## Running the Pipeline

```bash
# Prerequisites: seqkit in PATH, Python ≥ 3.10

python3 scripts/run_insilico_pcr.py       # build primer table + run seqkit amplicon
python3 scripts/pairwise_identity.py      # pairwise identity + cross-reactivity
python3 scripts/combinatorial_resolution.py  # cross-assay resolution matrix
python3 scripts/genbank_pcr.py            # GenBank supplement (add .gb files to test_gb/ first)
python3 scripts/generate_final_summary.py # integrate all → final_resolution_summary.tsv
```

## References

1. **Downs SL, Madhi SA, Nunes MC, Olwagen CP.** (2023). Optimization of a
   high-throughput nanofluidic real-time PCR to detect and quantify 15 bacterial
   species and 92 *Streptococcus pneumoniae* serotypes.
   *Scientific Reports*, 13, 4588.
   https://doi.org/10.1038/s41598-023-31820-4

2. **Olwagen CP, et al.** (2023). Lancet Microbe.
   (Additional assays: 19C, 19F Atypical, and 14 viral/bacterial targets;
   updated pool assignments.)
   https://doi.org/10.1016/S2666-5247(23)00260-4

3. **Lorenz O, et al.** (2025). SeroBA(v2.0) and SeroBAnk: a robust genome-based
   serotyping scheme and comprehensive atlas of capsular diversity in
   *Streptococcus pneumoniae*. *Microbial Genomics*, 11(10), 001483.
   https://doi.org/10.1099/mgen.0.001483
   GitHub: https://github.com/GlobalPneumoSeq/seroba

4. **seqkit** — A cross-platform and ultrafast toolkit for FASTA/Q file manipulation.
   https://doi.org/10.1371/journal.pone.0163962
