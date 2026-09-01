# Reproduce the comparative toxin-tolerance addendum

Updated: 2026-09-01

This workflow extends the base repository. It reruns the six-species orthology analysis, exact-state screens, candidate-family controls, Para checks, NHE2/3-labelled exchanger taxon panels, transcript checks, and a tree quality-control analysis.

The output is a hypothesis screen. It does not establish controlled grayanotoxin tolerance in *Apis laboriosa*, a causal gene, or a phenotype association.

## Fast path

Requirements:

- Linux x86-64
- micromamba, mamba, or conda
- internet access for the first input download
- about 2 GB free disk space after cloning
- four or more CPU cores

From a fresh clone:

```bash
git clone https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance.git
cd apis-laboriosa-grayanotoxin-tolerance
make -C comparative_addendum all
```

To use a mamba-compatible executable with another name:

```bash
make -C comparative_addendum all MAMBA=mamba
```

The command creates two environments: the small base-project Python environment at `.venv` and the pinned bioinformatics environment at `comparative_addendum/.env`.

## Stages

```bash
make -C comparative_addendum setup
make -C comparative_addendum base-inputs
make -C comparative_addendum fetch
make -C comparative_addendum verify-inputs
make -C comparative_addendum results
make -C comparative_addendum validate
make -C comparative_addendum verify-results
```

The `results` stage performs the following work:

1. verifies the base and addendum inputs;
2. selects the longest RefSeq protein per GFF gene;
3. constructs six-species orthogroups with OrthoFinder 3.1.5;
4. aligns 8,273 complete single-copy orthogroups and applies the same strict residue-sharing rule to each *Apis* focal species;
5. tests prespecified barrier, detoxification, and toxicokinetic categories against length- and divergence-matched permutations;
6. screens strict copy-number sharing across all 10,003 orthogroups;
7. checks all 16 predefined Para/Nav grayanotoxin-contact positions and all callable Para residues;
8. tests the NHE2/3-labelled exchanger lead against the older *B. terrestris* reference, seven more *Bombus* proteins, and five non-*Bombus* bees;
9. checks the four focal exchanger positions against the public *A. laboriosa* and *A. dorsata* worker transcript assemblies;
10. reports existing whole-worker abundance values as a descriptive comparison;
11. infers a maximum-likelihood exchanger tree for alignment and sequence-identity quality control;
12. validates the stated results and all recorded hashes.

Set the OrthoFinder process count when needed:

```bash
ADDENDUM_THREADS=8 make -C comparative_addendum results
```

The alignment and permutation stages use recorded seeds and fixed thread settings where output order or floating-point reduction could otherwise change files.

If OrthoFinder exits without a complete orthogroup table, the runner preserves the incomplete directory with a `.failed` suffix and makes one fresh retry with two processes. This handles truncated intermediate search files without modifying inputs or accepting partial output.

## Inputs

The analysis uses 18 files from the base repository snapshot and six added files. [`combined_input_manifest.tsv`](combined_input_manifest.tsv) records the complete set. [`input_sources.tsv`](input_sources.tsv) records the added URLs, accessions, byte counts, and SHA-256 hashes.

`scripts/fetch_inputs.py` tries the versioned GitHub Release archive first. It then tries official NCBI and Europe PMC locations if the archive is unavailable. A download is installed only after both byte count and SHA-256 match.

Europe PMC regenerates metadata in its outer archive. The fetcher extracts the stable publisher member `evab227_supplementary_data.zip` and verifies that file. This avoids tying reproduction to changing ZIP timestamps.

Force official locations:

```bash
comparative_addendum/.env/bin/python \
  comparative_addendum/scripts/fetch_inputs.py --official-only --force
```

Require the release snapshot:

```bash
comparative_addendum/.env/bin/python \
  comparative_addendum/scripts/fetch_inputs.py --snapshot-only --force
```

Verify without downloading:

```bash
comparative_addendum/.env/bin/python \
  comparative_addendum/scripts/verify_project.py --inputs
```

The added archive is `comparative-addendum-inputs-v1.tar`. Its expected digest is in [`input_snapshot.sha256`](input_snapshot.sha256). The uncompressed deterministic archive is about 19 MB.

## Pinned software

[`environment.yml`](environment.yml) records every requested package version. The principal versions are:

| Program | Version |
| --- | --- |
| Python | 3.11.16 |
| OrthoFinder | 3.1.5 |
| DIAMOND | 2.2.5 |
| FAMSA | 2.4.1 |
| IQ-TREE | 3.1.3 |
| BLAST | 2.17.0 |
| Biopython | 1.88 |
| SciPy | 1.17.1 |
| pyfamsa | 0.5.3 |

The tested platform is Linux x86-64. Exact solver build strings remain available in the environment installation transaction, while `environment.yml` pins the portable package versions required for a fresh solve.

## Independent checks

Compile scripts and check committed result files:

```bash
make -C comparative_addendum check
```

Check the scientific assertions directly:

```bash
comparative_addendum/.env/bin/python \
  comparative_addendum/scripts/validate_claims.py
```

Check all result, script, and generated-work manifests after a full run:

```bash
comparative_addendum/.env/bin/python \
  comparative_addendum/scripts/verify_project.py --results --scripts --work
```

The validator checks the exact values used in the addendum summary, including orthogroup counts, callable sites, focal and control hit counts, enrichment statistics, strict copy-number hits, Para results, exchanger-panel states, transcript residue matches, and expression ratios.

## Expected headline output

[`results/analysis_summary.json`](results/analysis_summary.json) should contain:

- 60,216 representative proteins;
- 59,381 proteins assigned to orthogroups;
- 10,003 orthogroups;
- 8,273 complete single-copy orthogroups;
- 4,902,144 callable amino-acid sites;
- 2,024 strict focal sites for *A. laboriosa*;
- 2,036 corresponding sites for the *A. dorsata* internal control.

[`results/para_strict_sharing_summary.json`](results/para_strict_sharing_summary.json) should report zero exact focal Para sites across 2,050 callable residues.

## Maintainer snapshot workflow

Downloaded inputs stay outside Git history. To rebuild the deterministic added-input archive after an intentional input revision:

```bash
make -C comparative_addendum snapshot
sha256sum comparative_addendum/comparative-addendum-inputs-v1.tar
```

Inspect the input manifest and all scientific result differences before updating `input_snapshot.sha256` or replacing the release asset.

To refresh committed result hashes after an intentional and inspected result change:

```bash
comparative_addendum/.env/bin/python \
  comparative_addendum/scripts/build_manifests.py \
  --refresh-results --refresh-scripts --refresh-work
```

Do not refresh a checksum solely to make a failed analysis pass.

## Interpretation limits

- *A. laboriosa* has no located controlled oral or injection grayanotoxin challenge.
- The same-toxin phenotype matrix has one directly measured tolerant lineage and one directly measured sensitive lineage.
- Exact residue sharing does not prove independent substitution, phenotype association, or function.
- The transcript comparison has one untreated pooled whole-worker library per species and no biological replication.
- Current Chinese reference samples may not represent Nepalese mad-honey populations.
- Ecological ingestion does not establish the dose reaching hemolymph, brain, or sodium channels.

See [`REPORT.md`](REPORT.md) for the full result and experiment design.
