# Agent instructions for the comparative addendum

## Objective

Maintain a reproducible cross-bee screen of coding patterns that could bear on grayanotoxin tolerance in *Apis laboriosa*.

This directory is an addendum to the base analysis. Its phenotype matrix is underpowered for a formal convergence-association model. Treat the exact-state screen as candidate generation.

## Read first

Read these files completely before changing an analysis or claim:

1. `README.md`
2. `REPORT.md`
3. `REPRODUCE.md`
4. `phenotype_matrix.tsv`
5. `hypothesis_matrix.tsv`
6. `input_sources.tsv`
7. the script that owns the result under discussion

Also read the repository-root `AGENTS.md`, `REPORT.md`, and `data_sources.tsv`. The root scientific guardrails apply here.

## Scientific boundary

- Controlled oral tolerance is directly measured for *Bombus terrestris audax*.
- Controlled oral sensitivity is directly measured for *Apis mellifera mellifera*.
- *A. laboriosa* is a focal ecological hypothesis taxon. No controlled grayanotoxin challenge or tissue internal-dose study was located.
- *A. dorsata*, *A. cerana*, and *A. florea* have unknown grayanotoxin phenotypes in this analysis. Do not label them susceptible.
- One directly measured tolerant lineage and one directly measured sensitive lineage are insufficient for PCOC, RERconverge, CSUBST, or a related phenotype-convergence claim.
- Exact shared residues may result from parallel change, reversal, ancestral retention, lineage history, alignment error, annotation error, or chance.
- The NHE2/3 label is uncertain across references. No direct grayanotoxin link has been shown for this exchanger.
- Whole-worker abundance values come from one untreated pooled library per species. Use the word descriptive. Do not call the ratios differential expression.
- Chinese genome and transcriptome samples may not represent Nepalese mad-honey populations.
- Ecological ingestion does not measure the grayanotoxin concentration reaching hemolymph, brain, or sodium channels.
- Absence of broad family expansion does not exclude regulation, localization, transport direction, substrate specificity, or a single specialized protein.

Keep observations, biological inferences, and proposed experiments distinct.

## Canonical run

From the repository root:

```bash
make -C comparative_addendum all
```

Expected checks:

- the base fetch verifies 25 files;
- the addendum fetch verifies six files;
- the primary-proteome table reports 60,216 proteins;
- OrthoFinder reports 10,003 orthogroups and 8,273 complete single-copy groups;
- `validate_claims.py` exits successfully;
- `verify_project.py --results --scripts --work` verifies all manifests.

Run the shorter committed-artifact check after documentation or validation-only edits:

```bash
make -C comparative_addendum check
```

If an input hash changes, stop and report the exact path, expected digest, observed digest, and source response. Do not silently revise a checksum.

## Extending the work

For every new hypothesis, record:

1. biological prediction;
2. phenotype labels and their evidence;
3. public data and stable accessions;
4. computational rule decided before inspecting hits;
5. internal or negative controls;
6. multiple-testing treatment where applicable;
7. machine-readable output;
8. result and evidence limit;
9. experiment that could reject the interpretation.

For every new input:

1. add its path, accession, URL, byte count, and SHA-256 to `input_sources.tsv`;
2. add deterministic retrieval to `scripts/fetch_inputs.py`;
3. ensure both official-source and release-snapshot modes verify it;
4. update `combined_input_manifest.tsv` through `scripts/build_manifests.py`;
5. record source terms in `THIRD_PARTY_DATA.md`.

For every changed result:

1. rerun the full pipeline;
2. inspect the result diff;
3. add or update a direct assertion in `scripts/validate_claims.py`;
4. refresh `results.sha256` only after the scientific diff is understood;
5. update the addendum README, report, and reproduction guide together;
6. update root documents if the overall conclusion changes.

## Result ownership

| Claim | Primary file | Generating script |
| --- | --- | --- |
| Orthogroup and callable-site totals | `results/analysis_summary.json` | `scripts/screen_strict_convergence.py` |
| Strict focal and internal-control sites | `results/strict_site_aggregate.tsv` | `scripts/screen_strict_convergence.py` |
| Candidate-family enrichment | `results/candidate_site_enrichment.tsv` | `scripts/screen_strict_convergence.py` |
| Category-wide exploratory tests | `results/candidate_category_fisher.tsv` | `scripts/screen_category_enrichment.py` |
| Strict copy-number sharing | `results/strict_copy_number_hits.tsv` | `scripts/screen_strict_convergence.py` |
| Current *Bombus* Para contact sites | `results/current_bombus_para_gtx_sites.tsv` | `scripts/check_current_bombus_para.py` |
| Whole-Para strict screen | `results/para_strict_sharing_summary.json` | `scripts/screen_para_strict_sharing.py` |
| Wider *Bombus* exchanger panel | `results/nhe3_bombus_panel_sites.tsv` | `scripts/screen_nhe3_bombus_panel.py` |
| External-bee exchanger panel | `results/nhe3_external_bee_sites.tsv` | `scripts/screen_nhe3_external_bees.py` |
| Transcript residue check | `results/nhe3_transcript_residue_validation.tsv` | `scripts/validate_nhe3_transcripts.py` |
| Descriptive pooled-worker abundance; GEO says whole body without belly | `results/nhe3_constitutive_expression.tsv`; interpretation corrected in `../transcriptomic_addendum/REPORT.md` | `scripts/screen_nhe3_expression.py` |

## Code and data rules

- Use Python 3.11 and the versions in `environment.yml`.
- Derive paths from `Path(__file__)`; do not depend on a caller's current directory.
- Keep random seeds explicit and stable.
- Use one pyfamsa worker with refinement disabled for every alignment. Realign complete single-copy orthogroups from the pinned primary proteins; do not use OrthoFinder's saved multiple-sequence alignments for the strict-state screen.
- Write result tables with fixed columns and stable ordering.
- Fail on missing candidates, malformed tables, ambiguous mappings, and checksum disagreements.
- Keep network access inside `scripts/fetch_inputs.py`.
- Do not commit downloaded genomes, proteomes, transcriptomes, raw reads, environment directories, OrthoFinder work directories, or snapshot tar files.
- Commit the compact generated evidence needed to audit each claim.
- Keep third-party data terms separate from the repository software license.

## Reporting style

State the measured result first. Put the main limitation in the same section. Use direct technical language. Avoid promotional claims, mechanistic certainty, and significance language for uncorrected exploratory results.
