# Reproducing the transcriptomic addendum

## Complete run

Requirements: Linux x86-64, Python 3.11 or newer, `make`, network access for
initial retrieval, four CPU threads, and about 4 GB of free disk space.

```bash
git clone https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance.git
cd apis-laboriosa-grayanotoxin-tolerance
make -C transcriptomic_addendum all
```

The workflow uses the root `requirements.txt` (Biopython 1.88, openpyxl 3.1.5)
and minimap2 2.31-r1302. `scripts/setup_tools.py` downloads the official Linux
binary archive and checks its SHA-256 before installing the executable.
The upstream software and source code are available from
[minimap2 v2.31](https://github.com/lh3/minimap2/releases/tag/v2.31).

This workflow consumes two committed comparative result tables. Their exact
hashes are checked in `dependencies.sha256`. To regenerate those tables from
the underlying genomes, use the separate
[comparative workflow](../comparative_addendum/REPRODUCE.md). Its OrthoFinder
run is unnecessary for rechecking this addendum against the frozen groups.

## Stages

```bash
make -C transcriptomic_addendum setup
make -C transcriptomic_addendum fetch
make -C transcriptomic_addendum verify-inputs
make -C transcriptomic_addendum results
make -C transcriptomic_addendum validate
make -C transcriptomic_addendum verify-results
```

`all` performs these sequentially. `results` rebuilds alignments and tables even
when old intermediates exist. It makes no network calls. `check` compiles the
scripts, validates claims from tables, and checks the output snapshot without
rerunning alignments.

## Input preservation

The ten added files total **317,928,612 bytes**. Each expected size and SHA-256
is in `input_sources.tsv`. Retrieval first uses
`transcriptomic-addendum-inputs-v1.tar` from the
[versioned release](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/transcriptomic-inputs-v1).
Its archive hash is in `input_snapshot.sha256`. Original repositories are the
fallback; to explicitly test those sources:

```bash
.venv/bin/python transcriptomic_addendum/scripts/manage_inputs.py fetch --official
```

Existing files are verified and reused. This command does not overwrite a
mismatch. Testing retrieval of a present file requires moving it aside first.
Do not edit a manifest to accept changed upstream bytes without investigating
and versioning a new snapshot.

The original base input archive is also required and is restored by `fetch`.
The new archive contains two 2,000,000-record R1 prefixes, source metadata,
two viral GenBank records, two fly protein FASTAs, and the independent study's
XML and selected-gene spreadsheet. It does not contain complete sequencing runs.

Prefix retrieval streams the original ENA FASTQ and stops after the fixed
number of complete records. Gzip output uses timestamp zero and compression
level six. Both compressed and decompressed SHA-256 values are preserved in
`results/raw_subset_metadata.json`. Gzip implementation differences can change
compressed bytes; the published snapshot is the exact-byte reproduction path.
Full-run MD5 values come from ENA and were **not** validated by a partial
download. No random seed applies because selection is a fixed prefix.

Rebuild the deterministic source archive from validated inputs:

```bash
make -C transcriptomic_addendum snapshot
```

## Analysis rules

| Stage | Parameters and counting unit |
| --- | --- |
| Unigene to RNA | minimap2 `-x asm5 -c --secondary=yes -N 50 -p 0.5 -t 4`; each assembly against both RNA references |
| Unigene acceptance | At least 200 aligned nt, 95% identity, 80% query coverage; repeat at 50% coverage |
| CDS construction | One frozen primary protein per gene; exact full protein translation in annotated RNA; one unique CDS sequence or recorded exclusion |
| Host R1 alignment | minimap2 `-x sr -c --secondary=yes -N 20 -p 0.8 -t 4` against both species' pooled primary CDSs |
| Host R1 acceptance | At least 100 aligned nt, 90% identity, 80% read coverage; best other orthogroup score below 98% of the best score |
| Counts | One read per accepted orthogroup; combine alternative species references and within-group paralogues |
| Host normalization | Median log2 AL/AD ratio across one-to-one groups with at least ten reads in both pools |
| Virus alignment | Two pinned BQCV genomes; `-x sr -c --secondary=yes -N 5 -p 0.8 -t 2`; same length, identity and coverage filters |
| Viral counts | Union of read identifiers across references; separately check credible competing bee-CDS alignments |
| Viral coverage | Best alignment to OR496406.1 per read; alignment spans, not exact per-base pileups |
| Viral contig | minimap2 `-x asm5 -c --secondary=yes -N 5 -t 2` against both viral references |
| Chitin-synthase comparison | Local protein alignment; BLOSUM62, gap open -10, gap extension -0.5; all fly references >=1,300 aa, bee sequences >=1,000 aa |

The host reference excludes 200 proteins without a unique recoverable exact
CDS. This and the use of coding sequence only leave some reads uncounted.
The viral competition check tests the included bee CDSs, not every possible
bee genomic or noncoding sequence. Near-full viral contig identity supplies
independent support for the viral assignment.

## Verification and expected results

The complete entry point reproduced all 24 first-pass output files exactly;
the additional raw-prefix audit brings the final manifest to **25 files**.
All original 33 root result files also retained their hashes.

`validate_claims.py` checks tissue metadata, viral counts and identity, focal
gene counts, normalization, defined versus undefined ratios, enzyme comparison,
selected-supplement restrictions, and decompressed raw-prefix hashes.
`verify.py --results` requires an exact output file set, exact SHA-256 values,
and consistent TSV column counts. Neither command estimates biological
uncertainty from technical reads.

Alignment logs and regenerable PAF files are in ignored `work/`. The viral
contig PAF and chitin-synthase protein alignments are committed results. Input
metadata, all gene counts, panel sensitivity, exclusions, and annotations are
retained so individual claims can be traced without rerunning the entire job.
