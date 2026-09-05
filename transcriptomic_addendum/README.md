# Transcriptomic addendum: sample composition and raw-read checks

Updated 2026-09-05. Computational analysis with Codex, using public data from
GEO, ENA, RefSeq, UniProt, and a published grayanotoxin experiment.

## Result

The earlier interpretation of a stronger gut barrier needs to be withdrawn.
Both source samples are described as **“whole body without belly.”** A large
viral RNA imbalance further confounds the comparison. The elevated
chitin-synthase protein resembles the fly's cuticle and trachea enzyme more
closely than its gut-matrix enzyme. These data cannot identify an enhanced
midgut barrier or a grayanotoxin-tolerance mechanism.

Several large RNA differences survive direct checks against the original
sequencing reads. They remain useful candidates for explaining the difference
between these particular pools. Species, location, tissue composition, and
physiological condition cannot be separated with this design.

## New tests

Each raw-read result uses the first **2,000,000 R1 reads per species** from
SRR9034695 and SRR9034696. These are technical subsets of one pooled library
per species. Zero counts refer only to the analyzed subset.

| Test | A. laboriosa | A. dorsata | Interpretation |
| --- | ---: | ---: | --- |
| Black queen cell virus matching reads | 188 (0.0094%) | 584,834 (29.2417%) | Major sample-composition difference; illness and causation unmeasured |
| Original peritrophin-like group OG0001109 | 2,950 | 0 | Strong RNA contrast survives direct read mapping; tissue and function unresolved |
| Additional inconsistently named chitin-binding group OG0001110 | 42,268 | 14 | Much larger candidate was missed by relying on gene names |
| Elevated chitin-synthase group OG0000293 | 465 | 2 | Protein is Kkv-like in the tested comparison; weak evidence for a gut program |
| ABCC/MRP group OG0000499 | 1,161 | 54 | Specific transporter candidate; 17.52-fold after host-gene normalization |
| OATP group OG0006494 | 197 | 6 | Additional transporter candidate; low comparator count limits precision |

The dominant dorsata assembly contig matches an 8,440-base BQCV reference over
8,437 aligned bases at **98.66% identity**. Raw reads cover the whole reference.
Viral abundance alone cannot explain the very large chitin-binding contrasts:
removing the measured viral fraction changes the dorsata denominator by about
1.41-fold. The biological consequences of that viral RNA are unknown.

Most adequately covered ABCC genes have little relative shift: the median
normalized gene ratio is **1.03** across seven one-to-one groups, while the
family's summed ratio is **2.23**. The candidate list should focus on individual
genes. No transporter here has been shown to move grayanotoxin.

The independent GTX-I study in *Spodoptera litura* supplies a selected
282-gene summary rather than the replicate count matrix. Of those rows,
43 have the supplied adjusted p-value below 0.05. The file cannot support a
fresh differential-expression or unbiased cross-species enrichment analysis.

## Reproduce

From the repository root on Linux x86-64, with Python 3.11+ and about 4 GB free:

```bash
make -C transcriptomic_addendum all
```

This verifies the original inputs and frozen orthology tables, installs pinned
minimap2, restores ten added source files, reruns the alignments and analyses,
checks headline claims, and compares every output with its committed SHA-256.
The existing OrthoFinder results are checked inputs; this command does not
rerun OrthoFinder. Its full independent workflow remains in the comparative
addendum.

The exact added inputs, including both read subsets, are preserved in the
[versioned release](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/transcriptomic-inputs-v1).
The complete original sequencing runs are not included.

## Files

| File | Purpose |
| --- | --- |
| [REPORT.md](REPORT.md) | Evidence, revised interpretation, and bounded next experiments |
| [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md) | Predictions and when follow-up tests were added |
| [REPRODUCE.md](REPRODUCE.md) | Commands, assignment rules, software, and verification |
| [AGENTS.md](AGENTS.md) | Instructions for extending the analysis |
| [input_sources.tsv](input_sources.tsv) | Ten source URLs, retrieval methods, sizes, and hashes |
| [dependencies.sha256](dependencies.sha256) | Exact committed orthology tables used |
| [panels.tsv](panels.tsv) | Functional grouping rules and their limits |
| [results/](results/) | All counts, alignments, source audits, and sensitivity analyses |
| [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md) | Source credit, terms, and archive contents |

The original 33 result files reproduce byte for byte. Their numerical values
remain available; this addendum corrects their tissue description and biological
interpretation. The negative Para sequence findings still stand within their
original scope. No controlled *A. laboriosa* grayanotoxin phenotype has been
established by this project.
