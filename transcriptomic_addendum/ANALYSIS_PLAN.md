# Transcript robustness investigation

Started 2026-09-05 from main commit 2f69bea6bb83beacba69c5f46d8e1315b3254d11.

This plan was written after reading the existing analysis and inspecting the
submitted abundance-table totals and top 15 transcripts. It is a prospective
analysis plan for the tests below, not an independent preregistration.

## Questions and predictions

| Hypothesis | Prediction | Test | Controls and limit |
| --- | --- | --- | --- |
| H1: The large peritrophin contrast depends on incomplete or asymmetric matching | Aligning all unigenes to all annotated transcripts changes the family contrast or finds overlooked dorsata homologues | Use minimap2 2.31 to map each assembly to both species' complete RefSeq RNA sets, collapse isoforms to genes and genes to the published orthogroups, and compare with the original exact-probe assignments | Own and swapped references; ambiguous gene/group hits excluded; examine aligned span and sequence identity; no new biological replication |
| H2: The gut candidates are part of a broad tissue or physiological difference between samples | Digestive, cuticle, gland, muscle, or metabolic marker groups shift together with the focal genes | Compare ortholog-matched expression ratios for predefined functional marker panels and genome-wide controls | Report effects and ranks, no species differential-expression p-values; functional annotations are not exclusive tissue labels |
| H3: Transcript length and library composition exaggerate the positive candidates | FPKM, count-per-length TPM, and normalization to typical one-to-one orthologues give materially different rankings or fold ratios | Audit submitted counts/FPKM against assembled lengths; compute alternative abundance measures; report sensitivity across mapping and normalization choices | RSEM effective lengths and gene-level isoform aggregation can differ from longest-unigene lengths; alternative estimates are sensitivity analyses |

## Fixed initial alignment and assignment rules

- Align nucleotide unigenes to annotated RNA, both orientations, with minimap2
  `-x asm5 -c --secondary=yes -N 50 -p 0.5 -t 4`.
- A high-coverage assignment requires at least 200 aligned nucleotides, at least
  95% aligned identity, and at least 80% query coverage in one alignment.
- Collapse multiple transcript isoforms belonging to the same gene. Collapse
  genes using the already published OrthoFinder group membership.
- Require the best alternative gene/group alignment score to be below 98% of
  the best score to call the assignment unique. Record unresolved assignments.
- Sum each unigene's abundance once within a group; do not count isoforms twice.
- Primary species comparison uses groups with one annotated gene in each of
  laboriosa and dorsata. Also retain family totals and original focal genes.
- Report relaxed query coverage (50%) as sensitivity analysis, and use direct
  pairwise alignments to inspect important mapping disagreements.
- Use only positive-abundance orthologues for log ratios. For normalization,
  take the median ratio among mapped one-to-one orthologues with at least 10
  submitted counts in each sample. Do not turn this into a biological p-value.

## Additional evidence

Seek an independent grayanotoxin exposure transcriptome and replicated bee
tissue/physiology data. Record accession and design before deciding whether a
comparison is valid. An exposure response in a poisoned insect cannot establish
the mechanism of tolerance in a different species.

All negative and conflicting results will be retained. New results go in this
directory and will not overwrite the original result snapshot.

## Follow-up rules recorded after the first alignment screen

- GEO describes BOTH samples as `whole body without belly`. Treat this exact
  metadata as controlling tissue interpretation. The precise dissection is
  unspecified. Do not infer a midgut expression measurement.
- Check the first 2,000,000 R1 FASTQ records of each original run (SRR9034695,
  SRR9034696). These are deterministic technical subsets. They are not random
  samples or biological replicates.
- Build one exact translated CDS per primary gene per species from the pinned
  proteomes, GFF and RNAs. Pool both species' CDSs as a symmetric reference.
  Map raw R1 reads with minimap2 2.31 `-x sr -c --secondary=yes -N 20 -p 0.8`.
  Require at least 100 aligned bases, 90% identity, 80% read coverage, and a
  greater than 2% score lead over any other orthogroup. Count each read once.
- Report raw count ratios and normalization to the median of one-to-one
  orthologues with at least 10 reads per sample. Treat missing counts explicitly.
- Trace the dominant unclassified dorsata transcript by a public BLAST search,
  then retrieve the matching reference for a local reproducible comparison.
  No organism identity will be inferred solely from absence of a bee match.
- The independent Spodoptera supplement is a preselected 282-gene summary,
  despite its title saying read counts. It does not contain a replicate count
  matrix or all tested genes. Audit the supplied adjusted p-values and record
  these restrictions before considering any cross-species enrichment.

## H5: Chitin-synthase annotation does not identify a midgut program

After inspecting orthogroup annotations and raw counts, the elevated group
OG0000293 included both `chs-2` and `kkv` names. The other group OG0000875 also
used `chs-2`. Test the prediction that the elevated gene is closer to the
experimentally characterized cuticle/tracheal enzyme Kkv (CHS-A) than to the
midgut peritrophic-matrix enzyme Chs2 (CHS-B). Retrieve Drosophila melanogaster
reference proteins for both genes. Compare both full-length laboriosa proteins
and homologues from other Apis species to both references with the same local
protein alignment parameters (BLOSUM62, gap open -10, gap extension -0.5).
Report both scores, identity, coverage, and exact accessions. Sequence-based
class assignment does not establish tissue localization in A. laboriosa.

Functional comparator source: Bertran-Mas et al. 2025,
https://doi.org/10.1371/journal.pgen.1011847.
