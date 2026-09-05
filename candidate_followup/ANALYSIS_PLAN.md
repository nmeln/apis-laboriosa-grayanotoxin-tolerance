# Independent data and ABCC candidate follow-up

Started 2026-09-05 from e65a11d on main. This plan follows the previous raw-read
screen that nominated OG0000499 (laboriosa XP_043798878.1, dorsata
XP_006621542.1). It is exploratory and does not constitute preregistration.

## Questions

| Hypothesis or data requirement | Prediction | Test and limitations |
| --- | --- | --- |
| The candidate depends on a deficient dorsata reference | Independent genome assemblies disagree on its structure or sequence | Inventory public assemblies, distinguish resubmissions from independent specimens, align the candidate to independent genomes and compare gene structure and coding support |
| The RNA contrast depends on asymmetric mapping or a small transcript segment | Reference swapping, shared callable sequence, or coverage checks erase the contrast | Examine per-base/read coverage, uniquely distinguishable regions, both-species references, and original assembled transcripts; retain low-expression and multimapping limits |
| The candidate's protein differs in a potentially relevant way | Laboriosa has changes in conserved or functionally interpretable regions | Compare bee homologues, independent genomes, and curated transporter references; report family assignment, domain architecture and all substitutions without inferring GTX transport from similarity |
| The candidate RNA difference recurs independently | Independent comparable laboriosa/dorsata samples reproduce the direction | Inventory RNA datasets and sample-level tissue metadata; analyze suitable original reads or matrices; maintain biological replication and batch boundaries |
| A gut or abdominal role is independently supported | Tissue-resolved data show expression in relevant organs | Look for laboriosa/dorsata first, then explicitly labelled other-bee evidence; avoid treating a related species as independent proof of laboriosa regulation |

Genome sequencing and RNA expression address different questions. High BQCV RNA
in SRR9034696 does not establish contamination of the independent reference
genome. No retrieved sample will be labelled virus-free without a corresponding
check; no negative sequence screen establishes absence of all pathogens.

The reason for removing the abdomen will be sought in the primary methods.
If not specified there, a plausible methodological explanation will be labelled
as an inference. Gut microbial or metagenomic DNA cannot substitute for host
RNA expression. Repeated runs of one library do not supply biological replicates.

All new data, exclusions, controls and negative findings will be retained with
source URLs and checksums. Additional tests motivated by outcomes will be
recorded as follow-ups. Retrieval and analysis remain separate.

## Follow-ups added after initial inspection

- After reference swaps retained the RNA contrast, add an aligner-free exact
  shared-marker count. Use the central 100 nt of each 150-nt R1 read; retain
  only markers shared by both candidate CDSs and absent from every competing
  primary CDS and both tested BQCV genomes. This is a technical sensitivity test.
- After finding the exact mellifera accession in Arora et al. 2025, extract
  its reported peptide totals and all replicate spectral counts with two
  other ABCC proteins as controls. Attempt peptide-export retrieval and preserve
  failures; do not imply a fresh raw-spectrum search from the summary workbook.
- After reconstructing both independent genomes, retain every coding difference
  and distinguish differences supported in both assemblies per species from
  within-species differences. Project human ABCC4 regions as inferred features.
- After three substitutions distinguish the two AL assemblies from the other
  five reference bees, define single AL-to-AD changes and their combination as
  untested sequence designs. This makes a future coding test explicit; it does
  not constitute a transport experiment or selection test.
