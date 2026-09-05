# Instructions for the transcriptomic addendum

Read the root AGENTS.md, this directory's README.md, REPORT.md,
ANALYSIS_PLAN.md, REPRODUCE.md, input_sources.tsv, and the relevant generating
scripts before extending the work.

## Scientific constraints

- The original gut-barrier interpretation is superseded by this source audit.
- Both GEO samples say `whole body without belly`. Preserve that exact wording.
  Do not describe these libraries as whole-bee or measured midgut RNA.
- One pool of three individuals per species is biological n=1 per species.
  Raw reads, assembly contigs, and technical subsets do not add replication.
- High BQCV RNA in dorsata is supported by sequence identity and raw reads.
  Do not infer active infection, symptoms, or a causal host-expression effect.
- Zero reads means zero in the specified prefix and mapping procedure. Preserve
  undefined ratios; do not invent pseudocount-based fold claims.
- Chitin-binding labels do not establish tissue or substrate. Kkv-like sequence
  similarity does not establish bee localization or a full gene phylogeny.
- Both tested laboriosa synthases favor Kkv in this comparison. Do not force
  the second into a gut-specific class by elimination.
- The strongest ABCC contrast is concentrated in OG0000499. Preserve the
  family sum and typical-gene distinction. GTX transport remains unmeasured.
- The Spodoptera supplement is a selected summary. Do not perform DESeq2 or
  unbiased pathway enrichment without replicate counts and the tested universe.
- No controlled laboriosa GTX phenotype is established here. Dorsata's GTX
  phenotype remains unknown. Viral read burden is not a GTX phenotype label.

## Reproduction

Run `make -C transcriptomic_addendum all` from the repository root. For a
results-only integrity check, use `make -C transcriptomic_addendum check`.
The two orthology tables are frozen checked dependencies. Read the comparative
addendum's instructions before changing their construction.

Keep downloads in `inputs/` and temporary alignments in `work/`, both ignored.
Add a source URL, retrieval method, byte count and SHA-256 for every new input.
Archive analyzed data as a versioned release asset. Record complete-run versus
subset boundaries explicitly. Never claim full-run verification from a prefix.

State each new hypothesis's prediction, data, test, result, and limitation.
Record tests added after inspecting outcomes. Keep exploratory exclusions and
negative results available. Avoid mechanistic claims based solely on a product
name, family total, or sequence similarity.

After a numerical change, inspect the result diff, rerun the entry point, run
claim validation, and only then intentionally version `results.sha256`.
Do not weaken a check to accept an unexplained discrepancy. Keep all downloads
outside analysis scripts. Do not overwrite the original numerical snapshot
to conceal its earlier interpretation; correct its accompanying documentation.
