# Agent instructions

## Objective

Maintain and extend a reproducible public-data analysis of proposed grayanotoxin-tolerance mechanisms in *Apis laboriosa*.

The repository is a hypothesis screen. It does not establish that *A. laboriosa* has systemic grayanotoxin tolerance, and it does not establish a causal mechanism.

## Read before changing anything

1. `README.md`
2. `REPORT.md`
3. `REPRODUCE.md`
4. `data_sources.tsv`
5. The script that generated any result being discussed

Do not infer a result from a filename or README sentence alone. Trace it to the generated TSV and the script that produced it.

For cross-bee orthology, convergence, or wider-taxon work, also read `comparative_addendum/AGENTS.md`, `comparative_addendum/REPORT.md`, and `comparative_addendum/phenotype_matrix.tsv` completely.

For any expression or physiological interpretation, read `transcriptomic_addendum/AGENTS.md` and `transcriptomic_addendum/REPORT.md`. That addendum supersedes the original gut-barrier interpretation and corrects the tissue metadata.

For the specific ABCC candidate, independent dorsata assemblies, tissue evidence or proposed sequence variants, also read `candidate_followup/AGENTS.md` and `candidate_followup/REPORT.md`. Its gut-protein evidence concerns mellifera; its variant FASTAs are untested designs. It does not establish GTX transport or a replicated laboriosa expression effect.

## Scientific guardrails

- Keep observation, inference, and speculation separate.
- Describe the pooled-worker transcriptome comparison as descriptive. It has one untreated pooled library per species and no biological replication. GEO specifies “whole body without belly”; do not call it whole-bee or defined midgut RNA. Preserve the BQCV imbalance and other sample confounding described in the transcriptomic addendum.
- Do not call the expression ratios differential expression.
- Do not claim that peritrophin blocks grayanotoxin. GTX-I and GTX-III are much smaller than molecules known to cross the honeybee peritrophic matrix.
- Do not describe scaffold 8 or scaffold 25 candidate identities as direct MAKER-model assignments. The relevant eastern-assembly GFF is unavailable. They are scaffold-and-GO candidates.
- Para exclusion from the selected scaffold set is stronger than candidate identification: 29 independent Para probes map to scaffold 105.
- Do not claim that the processed transcriptome excludes rare or tissue-specific RNA editing.
- Do not claim that absence of a broad gene-family expansion excludes altered regulation, localization, transport direction, or substrate specificity.
- Preserve the population caveat. The available Chinese samples may not represent Nepalese mad-honey populations.
- Preserve the premise caveat. Ecological exposure and toxic honey production do not measure hemolymph or brain exposure.
- Negative results apply to the tested reference sequences, annotations, isoforms, mapped sites, and processed transcript consensus.

## Reproduction workflow

For a fresh clone:

```bash
make all
```

Equivalent explicit stages:

```bash
make setup
make fetch
make verify-inputs
make results
make validate
make verify-results
```

Run the comparative addendum independently with:

```bash
make -C comparative_addendum all
```

Do not label *A. dorsata*, *A. cerana*, or *A. florea* grayanotoxin-susceptible. Their phenotypes are unknown in the addendum evidence matrix.

Expected behavior:

- `make fetch` downloads files only from the public sources encoded in `scripts/fetch_inputs.py` and `data_sources.tsv`.
- `make verify-inputs` must report all 25 analyzed input files as exact SHA-256 matches.
- `make results` must run every analysis in the order declared by `scripts/run_all.py`.
- `make validate` must derive the headline values from generated TSV files and fail on disagreement.
- `make verify-results` must match every committed file listed in `results.sha256`.

If an upstream source changes bytes for a pinned accession, stop. Record the mismatch. Do not silently replace the checksum or committed results.

`scripts/build_input_manifest.py` is a maintainer audit helper. By default it
writes ignored `local_input_manifest.tsv`. It refuses to overwrite
`data_sources.tsv`.

## Extending the analysis

When adding a dataset:

1. Use a stable accession, DOI, or official repository URL.
2. Add its expected path, source, byte size, and SHA-256 hash to `data_sources.tsv`.
3. Add deterministic retrieval logic to `scripts/fetch_inputs.py`.
4. Add the analysis as a focused script under `scripts/`.
5. Write outputs under `results/` in TSV, JSON, or plain text.
6. Add direct validation to `scripts/validate_claims.py` when the output changes a headline claim.
7. Rerun the full pipeline.
8. Regenerate `results.sha256` only after inspecting the diff and confirming the scientific reason for every changed file.
9. Update `README.md`, `REPORT.md`, and `REPRODUCE.md` together.

When testing a hypothesis, state all four items:

1. prediction
2. public data used
3. computational test
4. result and limitation

## Code rules

- Python 3.11 or newer.
- Keep dependencies pinned in `requirements.txt` and `environment.yml`.
- Prefer the standard library for orchestration, downloading, and verification.
- Keep scripts deterministic. If randomness becomes necessary, expose and record a seed.
- Use repository-relative paths derived from `Path(__file__)`.
- Never commit downloaded genomes, proteomes, transcriptomes, raw reads, temporary archives, `.venv`, or `.tools`.
- Keep third-party source terms and citations in `THIRD_PARTY_DATA.md`. Do not
  apply the repository software license to downloaded inputs.
- Write tabular outputs with stable column order and consistent line endings.
- Fail loudly on missing files, missing candidates, unexpected row counts, checksum mismatches, or ambiguous mappings.
- Do not weaken an assertion solely to make the pipeline pass.
- Avoid network access inside analysis scripts. Retrieval belongs in `scripts/fetch_inputs.py`.

The exact inputs are published as a GitHub Release asset. After an intentional
input update, rebuild the archive with `make snapshot`, inspect the manifest
and result diffs, update `input_snapshot.sha256`, and replace the versioned
release asset.

## Verification after edits

Run:

```bash
make check
```

For changes that affect analysis outputs, run:

```bash
make verify-inputs
make results
make validate
make verify-results
```

Inspect the result diff before committing:

```bash
git diff -- results/ results.sha256 README.md REPORT.md REPRODUCE.md
```

## Result ownership

Each headline claim has a primary machine-readable source:

| Claim | Primary result file |
| --- | --- |
| Six-bee Para site invariance | `results/gtx_target_residues.tsv` |
| 99.9025% laboriosa/dorsata Para identity | `results/laboriosa_pairwise_identity.tsv` |
| N95 and V465 laboriosa-specific Para residues | `results/laboriosa_unique_genomewide.tsv` |
| Seventeen laboriosa full-length Para isoforms | `results/para_isoform_species_summary.tsv` |
| Shared DIII S3-S4 haplotypes | `results/para_diii_s3_s4_haplotype_overlap.tsv` |
| Zero nonsynonymous laboriosa transcript mismatches | `results/para_transcriptome_alignment_summary.tsv` |
| Para on eastern scaffold 105 | `results/para_eastern_scaffold_probe_hits.tsv` |
| Scaffold 8 and 25 candidate cross-map | `results/population_sodium_candidate_crossmap.tsv` |
| DSC1 mapped-site invariance | `results/selected_60e_gtx_site_residues.tsv` |
| Copy-number screens | `results/detox_family_counts.tsv` and `results/barrier_clearance_gene_counts.tsv` |
| Chitin-locus correction | `results/chitin_synthase_locus_qc.tsv` |
| Original pooled-worker expression ratios | `results/constitutive_expression_contrasts.tsv`; interpretation corrected by `transcriptomic_addendum/REPORT.md` |

## Reporting style

Use direct technical language. Give the result before interpretation. Put the limitation next to the result. Avoid promotional language, certainty unsupported by the test, and mechanistic storytelling beyond the data.
