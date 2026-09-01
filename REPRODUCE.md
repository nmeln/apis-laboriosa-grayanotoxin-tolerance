# Reproducing the *Apis laboriosa* grayanotoxin analysis

Updated: 2026-09-01

## Fast path

Requirements:

- Linux or macOS
- Python 3.11 or newer
- internet access to the public repositories listed below
- about 1 GB of free disk space

From a fresh clone:

```bash
make all
```

This command:

1. creates `.venv` and installs the two pinned Python packages;
2. downloads the versioned input snapshot, with official repositories as a fallback;
3. checks every input byte count and SHA-256 hash;
4. runs the 12 analysis programs in a fixed order;
5. checks the headline numbers directly from the generated tables; and
6. checks every generated file against `results.sha256`.

The analysis uses no randomness.

## Run one stage at a time

```bash
make setup
make fetch
make verify-inputs
make results
make validate
make verify-results
```

`make check` is the short verification command for a clone whose inputs and
results are already present. It compiles every Python script, validates the
headline claims, checks every committed result hash, and checks TSV column
counts.

## Inputs

`data_sources.tsv` defines the analyzed input snapshot. Each row records:

- repository-relative destination;
- accession or dataset identifier;
- official source page or download URL;
- expected byte count; and
- expected SHA-256 hash.

`scripts/fetch_inputs.py` first tries the versioned GitHub Release archive. The
archive contains the exact 25 files used for the committed result tables. If
the archive is unavailable, the script retrieves the files from:

- NCBI RefSeq for genomes, annotations, proteomes, transcript sequences, rat
  Nav1.4, and the bee Para panel;
- Genome Warehouse for the independent eastern-Yunnan *A. laboriosa* assembly;
- NCBI GEO for GSE130963 processed unigene and expression files;
- Oxford University Press for the Cao et al. (2023) supplementary workbook; and
- UniProt for rat Nav1.4 topology annotations.

Large third-party files are excluded from Git history. The release snapshot is
about 266 MiB because most members are already compressed.

The fetcher writes each download to a temporary file. It installs the file only
after its byte count and SHA-256 hash match `data_sources.tsv`. A changed
upstream file causes a hard failure and leaves the expected destination
unchanged.

To replace even a valid cached file:

```bash
.venv/bin/python scripts/fetch_inputs.py --force
```

To bypass the release mirror and test every official source:

```bash
.venv/bin/python scripts/fetch_inputs.py --official-only --force
```

To require the release mirror:

```bash
.venv/bin/python scripts/fetch_inputs.py --snapshot-only --force
```

The archive itself is deterministic and is checked against
`input_snapshot.sha256`. Each extracted member is then checked against
`data_sources.tsv`.

## Analysis order

`scripts/run_all.py` is the canonical pipeline. It runs:

1. known Para grayanotoxin-site mapping and whole-channel comparison;
2. Para isoform extraction and transmembrane-block comparison;
3. exact cross-species DIII S3-S4 haplotype comparison;
4. processed Para transcript-consensus comparison;
5. extraction of Cao et al. Table S6;
6. exact-probe mapping of Para into the eastern assembly;
7. exact-probe cross-mapping of genes on selected scaffolds 8 and 25;
8. DSC1/60E and Nach candidate analysis;
9. detoxification-family copy-number screening;
10. barrier and clearance-family copy-number screening;
11. chitin-synthase locus quality control; and
12. descriptive mapping into the submitted pooled-worker expression tables.

Run it directly after inputs are present:

```bash
.venv/bin/python scripts/run_all.py
```

This form verifies the inputs before running and verifies the outputs when it
finishes.

## Independent checks

Headline assertions:

```bash
.venv/bin/python scripts/validate_claims.py
```

The script checks the exact values stated in `README.md`, including site
invariance, sequence identity, isoform count, transcript mismatches, scaffold
mapping, chitin-locus quality control, and expression ratios.

Input and result hashes:

```bash
.venv/bin/python scripts/verify_project.py --inputs
.venv/bin/python scripts/verify_project.py --results
```

The result check also rejects extra untracked files in `results/` and malformed
TSV row widths.

## Conda alternative

```bash
conda env create -f environment.yml
conda activate apis-laboriosa-gtx
python scripts/fetch_inputs.py
python scripts/run_all.py
```

## Interpretation limits

- The expression screen has one untreated pooled whole-worker library per
  species. Its FPKM and CPM ratios are descriptive.
- Processed whole-worker transcripts can reveal a fixed abundant consensus
  change. They cannot exclude rare or tissue-specific RNA editing.
- The eastern assembly lacks the MAKER GFF used for the published model IDs.
  Candidate identities on scaffolds 8 and 25 are scaffold-and-GO inferences.
- The sequence and annotation checks do not replace an oral or injection
  dose-response, tissue toxin measurements, or expressed-channel assays.

See `REPORT.md` for the results, evidence bounds, and source list.
