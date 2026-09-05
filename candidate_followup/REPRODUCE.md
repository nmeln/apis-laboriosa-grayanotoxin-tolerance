# Reproduce the candidate follow-up

## Requirements

- Linux x86-64
- Python 3.11 or newer; the independent GitHub runner uses Python 3.12
- About 5 GB free disk space
- Four logical CPU threads are useful; two jobs run with two threads each
- Internet access for the first fetch

Pinned analysis software:

| Tool | Version | Use |
| --- | --- | --- |
| Biopython | 1.88 | FASTA/GFF-derived sequence processing, translation, pairwise alignment |
| openpyxl | 3.1.5 | Read the published proteomics workbook |
| pyfamsa | 0.5.3 | Eight-sequence protein alignment; one worker, refinement disabled |
| minimap2 | 2.31-r1302 | Full-mRNA genome alignment and short-read reference swaps |
| miniprot | 0.18-r281 | Spliced protein-to-genome reconstruction |

The two executables are downloaded from official releases and checked against
fixed SHA-256 digests before installation. Python packages are pinned in
`requirements.txt`. No stochastic statistical model is used.

## Full run

From a fresh clone of the repository:

```bash
make -C candidate_followup all
```

Stages:

1. Create the base virtual environment and install the pinned dependencies.
2. Restore the 25 original inputs, the transcriptomic addendum's ten inputs,
   and this follow-up's 28 added source files.
3. Verify input byte counts, SHA-256 hashes, and eight frozen code/table dependencies.
4. Rebuild the pooled 19,632-CDS host reference with the preceding addendum's script.
5. Audit the saved assembly and RNA inventories.
6. Align both full candidate proteins and mRNAs to four genome assemblies;
   reconstruct and translate the genomic exon sequences.
7. Map both 2,000,000-read R1 prefixes to each species' complete primary CDS reference.
8. Count unique orthogroups, compute descriptive host normalization, and inspect coverage.
9. Count exact shared 100-nt candidate markers without an aligner.
10. Compare eight bee sequences and the retrieved human/fly ABC panels;
    project human ABCC4 features and recover the expressed AL ORF.
11. Generate explicitly untested substitution-panel sequences.
12. Extract candidate and comparison records from the gut proteomics workbook.
13. Assert the scientific headline values and verify every generated result hash.

All source retrieval occurs in fetch/setup scripts. Analysis scripts use local
checked files. The result pipeline does not run OrthoFinder again; it verifies
the existing frozen orthogroup membership and primary-protein tables. The
comparative addendum provides the complete orthology reconstruction separately.

## Separate stages

```bash
make -C candidate_followup setup
make -C candidate_followup fetch
make -C candidate_followup verify-inputs
make -C candidate_followup results
make -C candidate_followup validate
make -C candidate_followup verify-results
```

For a check of already generated files:

```bash
make -C candidate_followup check
```

## Data preservation

The [versioned input archive](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/candidate-inputs-v1)
contains exactly the 28 entries in `input_sources.tsv`. Its digest is recorded
in `input_snapshot.sha256`. The archive includes complete additional genome
files and the complete original proteomics workbook. Previously archived raw
R1 prefixes are restored from the transcriptomic release.

The manifest distinguishes analysis/supporting sources from retained
exploratory sources. The latter preserve searches that did not become tests,
including broad initial protein queries and an older protein-atlas article.
They are not claimed as independent confirmations.

The default fetch uses the archive so that changing database-query responses
cannot silently change the analysis. Original URLs can be tried explicitly:

```bash
../.venv/bin/python scripts/manage_inputs.py fetch --official
```

Run that command from this directory. A source that now returns different
bytes will fail verification. Preserve the discrepancy and version a new
analysis if it is scientifically warranted; do not silently update hashes.
The saved RNA and assembly inventory is dated 2026-09-05, not a live inventory.

To rebuild the identical archive locally:

```bash
make -C candidate_followup snapshot
```

The tar file has fixed file order, permissions, owners and timestamps. It
contains unmodified source bytes. The command refuses an unexplained archive
checksum change. Do not commit the tar file or downloaded inputs to Git.

## Result ownership

| Output | Script |
| --- | --- |
| Assembly/RNA inventories and source distinctions | `audit_inventory.py` |
| Four-genome models, exon coordinates, coding DNA and translated proteins | `analyze_genomes.py` |
| Reference-swap counts, read assignments and coverage | `align_references.py`, `analyze_reference_swaps.py` |
| Exact markers and counts | `exact_shared_markers.py` |
| Bee alignment, coding differences, source-family comparisons and expressed ORF | `analyze_protein.py` |
| Untested single and combined substitution designs | `design_variants.py` |
| Gut-protein accession match and replicate spectral counts | `analyze_gut_proteomics.py` |

`results.sha256` covers all generated files. `validate_claims.py` checks
headline numbers and interpretation-relevant invariants. `source_retrieval_log.json`
records successful and failed source retrievals.

## Boundaries that must survive reproduction

- One original untreated pool per species remains biological n=1.
- The two-million-read prefixes are technical subsets, not independent replicates.
- The exact-marker ratio and normalized alignment ratios use different counting rules.
- Two assemblies per species do not establish population fixation.
- Membrane-region predictions and human feature projections do not measure bee function.
- The gut counterpart is from mellifera, not a laboriosa gut sample.
- Gut spectral counts and reported peptide totals come from the authors' workbook.
  A fresh raw-spectrum or peptide-sequence validation was not completed.
- Variant FASTAs are untested designs and omit terminal stop codons.
- No GTX transport, causal protection or controlled laboriosa tolerance is established.
