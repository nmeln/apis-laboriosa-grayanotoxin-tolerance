# Instructions for the candidate follow-up

Read the root instructions, this directory's README, REPORT, ANALYSIS_PLAN,
REPRODUCE, input manifest and relevant generating scripts before extending it.
For RNA interpretation, read the transcriptomic addendum's instructions and
report. For orthology changes, read the comparative addendum's instructions.

## Evidence constraints

- `GCA_000469605.1` and `GCF_000469605.1` are paired versions of one assembly.
  The Thai `GCA_009792835.1` is an independent source. Both are scaffold-level
  whole-genome assemblies, not gapless finished chromosomes.
- The BQCV observation concerns the original Chinese RNA pool. It does not
  establish a problem with the independent Malaysian or Thai DNA assemblies.
  Do not call their source animals virus-free.
- The original tissue is exactly `whole body without belly`. The anatomical
  boundary is not specified. A gene's DNA can be recovered from head/thorax;
  tissue activity requires RNA/protein measurements from that tissue.
- The paper's stated contamination-avoidance reason concerns DNA extraction.
  The RNA-dissection reason is not separately given. Label explanations about
  microbes and food as inference.
- The archive searches are a dated inventory with an explicit scope, not proof
  that no other data exist anywhere.
- Reference swaps and exact markers use the same two unreplicated pools.
  Do not call their ratios differential expression or adaptation evidence.
- Different marker and alignment counting rules must keep their own denominators.
- L254, L549 and F1134 are supported sequence candidates in two AL assemblies;
  this is not evidence of population fixation, selection or GTX protection.
- The Q447 gap is in a variable alignment region. Do not assert a structural
  consequence from its precise gap placement.
- Human ABCC4 feature coordinates are homology projections. The fly panel has
  closer hits than the name-matching lethal(2)03659 protein. Similarity ranking
  is not a resolved gene tree or a basis for specific function transfer.
- Mellifera `XP_393750.5` is detected in the authors' gut-membrane tables.
  This cannot locate AL's belly-excluded RNA signal or demonstrate GTX transport.
- The 31 unique-peptide value is a stage-level table entry, not 31 independently
  unique peptides in every replicate. Preserve the failed peptide-download
  record and the S4/S6 abundance inconsistency. Do not claim a new spectral search.
- The variant panel contains untested sequence designs. Do not present it as
  synthesized constructs, expressed proteins, measured transport, or a validated
  expression protocol. CDS FASTAs omit terminal stop codons.

## Reproduction and changes

Run `make -C candidate_followup all` for full reproduction. Run
`make -C candidate_followup check` for committed-result checks. Keep all
downloaded inputs under ignored `inputs/` and temporary files under `work/`.
Archive unmodified analyzed inputs in a versioned release.

Keep Python and executable versions pinned. Use pyfamsa with one worker and
refinement disabled. Sequence extraction must use genomic bases and exon
phases; do not fill unsupported bases with query residues.

For every added hypothesis, state the prediction, public data, test, result
and limit. Keep negative results and controls. No significance test of the
two original pools is justified by additional reads or normalization.

Add source URL, byte count and SHA-256 to the manifest. Existing input or
result mismatches are errors to investigate, never a reason to weaken an
assertion. After changing an analysis, inspect the result diff, validate the
claims, rerun the entry point, then deliberately version new checksums.
Update README, REPORT and REPRODUCE together when the conclusions change.

Use factual prose and exact identifiers. Do not turn a general transporter
annotation, a docking score or an indirect stress response into a GTX mechanism.
