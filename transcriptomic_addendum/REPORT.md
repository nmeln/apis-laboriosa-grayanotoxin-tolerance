# Transcriptomic robustness and sample-composition analysis

2026-09-05. Starting repository commit:
`2f69bea6bb83beacba69c5f46d8e1315b3254d11`.

## Current conclusion

The earlier gut-barrier interpretation exceeded the evidence. Direct source
metadata excludes an uncomplicated whole-bee interpretation; the reads reveal
a large black queen cell virus (BQCV) difference; and the elevated
chitin-synthase sequence favors a cuticle/trachea comparator. A specific
transporter mechanism remains possible. These samples cannot distinguish it
from species, collection-site, dissection, age/task, or health differences.

This work establishes several reproducible observations about the public
dataset. It does not establish grayanotoxin tolerance in *A. laboriosa*, a
causal role for any transcript, or a viral explanation for all host differences.
The existing negative Para sequence results reproduce and retain their stated
limits.

## Design and source audit

The original pipeline was rerun before extending it. All 25 original input
files matched their pinned hashes and all 33 original result files reproduced
byte for byte. The new analysis uses the original RefSeq annotations, RNA and
protein sequences, submitted Trinity assemblies/counts, and frozen comparative
orthogroup membership. Ten additional source files are in `input_sources.tsv`.

| Item | Observation | Consequence |
| --- | --- | --- |
| Biological replication | Three individuals contributed equal RNA to one pool per species | Biological n=1 per species; no species-expression p-values |
| Tissue description | Both GEO samples say “whole body without belly” | Exact dissection is unspecified; midgut expression was not measured as a defined tissue |
| Collection altitude | 2,150 m for laboriosa; 570 m for dorsata | Species and collection conditions are confounded |
| Exposure | No controlled grayanotoxin treatment in this comparison | Neither induction nor toxin protection is tested |
| Raw-read subset | First 2,000,000 R1 records of each original run | Deterministic technical check; no random sampling, R2 analysis, or biological replication |

Sources: [GSM3757258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3757258),
[GSM3757259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM3757259),
[PRJNA542114](https://www.ebi.ac.uk/ena/browser/view/PRJNA542114).
The complete retrieved GEO record and ENA run metadata are preserved in the
input snapshot. “Belly” suggests exclusion of the abdomen, but the record does
not define the anatomical boundary. We do not infer which small attached
structures may have remained.

## Hypotheses tested

| Hypothesis | Prediction | Test | Result |
| --- | --- | --- | --- |
| H1: Original mapping exaggerates the focal contrast | Symmetric mapping or raw reads erase the signal | All unigenes against each species' RNA; pooled two-species CDS reference for raw reads | Some assignments depend on reference and filters; original peritrophin contrast survives raw-read mapping |
| H2: A tissue/physiology difference accompanies the candidates | Other functional groups also shift | Fixed functional panels and source tissue audit | Muscle, gland-associated, immune, and chitin-related totals differ; no unique tissue or physiological cause identified |
| H3: Library composition explains the apparent elevation | Alternative normalization removes large contrasts | Submitted FPKM, rescaled TPM, count/length TPM, counts, and host-orthologue normalization | Magnitudes vary; large chitin-related and specific transport candidates persist |
| Follow-up: Dominant dorsata contig is viral | Near-full viral sequence identity plus many independent matching reads | Two pinned BQCV genomes, local contig alignment, raw-read coverage and host competition | 98.66% contig identity; 29.24% of dorsata R1 prefix matches BQCV |
| H5: Elevated chitin synthase is being assigned a gut function from its name | Sequence favors the cuticle/trachea comparator | Both full-length laboriosa proteins and other Apis homologues against both fly enzymes and all long reference isoforms | Elevated XP_043801530.1 strongly favors Kkv; bee tissue localization remains untested |
| Independent GTX response supports cross-species validation | Public replicate counts and tested background permit reanalysis | Audit the Spodoptera article and supplied spreadsheet | Selected summary only; no valid independent validation from that file |

`ANALYSIS_PLAN.md` records which tests followed inspection of earlier results.
This was an exploratory sequential investigation. The functional panels were
fixed for this screen; later viral and enzyme-class tests were motivated by
specific observations.

## 1. Viral RNA is a major unmodeled difference

The largest submitted dorsata contig, `AD|c34801_g1`, is 8,446 bases long.
It accounts for 14,177,610 of 61,355,778 submitted counts, or **23.1072%**.
A public sequence search identified BQCV as a candidate match. The reproducible
pipeline confirms the identity locally, without depending on that search.

| Reference | Identical/aligned bases | Identity | Contig coverage |
| --- | ---: | ---: | ---: |
| OR496406.1 | 8,324 / 8,437 | 98.6607% | 99.8934% |
| KY741959.1 | 8,311 / 8,437 | 98.5066% | 99.8934% |

References: [OR496406.1](https://www.ncbi.nlm.nih.gov/nuccore/OR496406.1),
[KY741959.1](https://www.ncbi.nlm.nih.gov/nuccore/KY741959.1).
The two matching-reference counts are combined by read identifier, never added
as though they were different reads.

| Raw-read check | Laboriosa | Dorsata |
| --- | ---: | ---: |
| Input R1 reads | 2,000,000 | 2,000,000 |
| Reads matching either BQCV reference | 188 | 584,834 |
| Fraction of input | 0.0094% | 29.2417% |
| OR496406.1 bases spanned by alignments | 7,953 / 8,440 | 8,440 / 8,440 |
| Credible competing bee-CDS reads | 0 | 0 |

This is strong evidence for abundant BQCV RNA in the sequenced dorsata pool.
It does not measure individual infection prevalence, active replication,
symptoms, or effects on host expression. Very low laboriosa counts could also
reflect trace contamination. BQCV detection in bees is established biology;
the result here is its measured contribution to these specific public data.

The submitted-contig fraction and raw-read fraction have different denominators
and assignment methods. They are separate estimates. Removing 29.2417% viral
reads changes the dorsata denominator by about 1.413-fold, so simple viral
dilution cannot explain several-hundred-fold host-gene contrasts. Broader
physiological effects remain inseparable from species in this comparison.

## 2. Symmetric raw-read mapping retains large, specific candidates

One primary protein per gene was taken from the frozen comparative analysis.
Its complete amino-acid sequence had to match an annotated RNA translation
exactly. This produced **19,632 CDS records** across the two species, with
**200 exclusions** documented individually. Both species' CDSs were pooled.
Each R1 read could contribute once to one orthogroup, after alignment-quality
and competing-group filters. Gene families sum grouped reads once; functional
panels can overlap.

The median laboriosa/dorsata count ratio across **6,157** adequately covered
one-to-one groups was **1.22727**. The following normalized ratios divide by
that value. This controls a typical host-gene shift; it assumes most eligible
groups provide a usable reference and does not prove equal RNA per bee.

| Group | Laboriosa protein | Label in pinned annotations | AL reads | AD reads | Normalized ratio |
| --- | --- | --- | ---: | ---: | ---: |
| OG0001109 | XP_043794209.1 | Peritrophin-1-like | 2,950 | 0 | Undefined |
| OG0001110 | XP_043794120.1 | PERK9 / peritrophin / early nodulin labels differ among species | 42,268 | 14 | 2,460.04 |
| OG0000293 | XP_043801530.1 | Chitin synthase, Kkv-like comparator | 465 | 2 | 189.44 |
| OG0000499 | XP_043798878.1 | Probable MRP lethal(2)03659 | 1,161 | 54 | 17.52 |
| OG0006494 | XP_043791287.1 | Organic anion transporter 74D | 197 | 6 | 26.75 |
| OG0005056 | XP_043791252.1 | Organic anion transporter 1A5 | 163 | 0 | Undefined |

The original peritrophin gene is short, 101 amino acids, so a small number of
distinct transcript regions can carry many reads. Its assignment table is
retained. Read counts measure library abundance, not independent molecules or
bees. The new 204-amino-acid OG0001110 protein is supported by the bee genome
and exact RNA translation. Inconsistent product labels must not be treated as
proof of plant origin, kinase activity, or gut localization.

| Functional panel | AL reads | AD reads | Normalized sum ratio | Median normalized ratio of eligible one-to-one groups |
| --- | ---: | ---: | ---: | ---: |
| ABCC | 2,316 | 847 | 2.228 | 1.028 (7 groups) |
| OATP | 640 | 168 | 3.104 | 1.197 (4 groups) |
| Chitin-binding/peritrophin labels | 45,270 | 86 | 428.915 | 34.416 (2 groups) |
| Muscle contractile | 17,227 | 28,127 | 0.499 | 0.572 (18 groups) |
| Royal-jelly labels | 24,475 | 4,706 | 4.238 | 2.939 (1 group) |

The ABCC excess is concentrated. After removing OG0000499, remaining ABCC
counts are 1,155 versus 793, a normalized ratio of about **1.19**. The group
therefore remains a focused candidate; broad ABCC activation is unsupported.
Transport direction, tissue localization, and GTX transport are unknown.
OATP ratios with six or zero comparator reads are particularly imprecise.

Panel patterns cannot diagnose age, task, health, or tissue proportions.
Chitin-binding annotations include proteins associated with structures outside
the gut. Other shifts, including a large CYP4G-related signal, also fit several
forms of physiological difference. No toxin-metabolism function is established.

### Assembly and normalization sensitivity

Mapping all unigenes to both species' complete RNA references exposed
reference-dependent assignments. At the 80% query-coverage threshold, only
about 35% of laboriosa and 26% of dorsata submitted counts mapped uniquely to
orthogroups using their own references; relaxing coverage to 50% raised these
to about 72% and 47%. These are incomplete subsets, especially when long
assembled UTRs do not align. Their family totals cannot replace direct raw-read
checks without qualification.

The alternate-reference screen found a low-abundance dorsata unigene for the
original peritrophin group (`AD|c31025_g1`, submitted FPKM 6.02). Failure to pass
one reference/filter combination is therefore not evidence that the species
lacks the transcript. All reference rotations and normalization variants are
retained in `panel_sensitivity.tsv`, `mapping_summary.tsv`, and
`orthogroup_expression.tsv`. The raw prefix's zero remains a subset observation.

## 3. Chitin-synthase name was an unreliable functional cue

In *Drosophila*, direct functional experiments distinguish Kkv (CHS-A,
cuticle/trachea) from Chs2 (CHS-B, peritrophic matrix). We compared both
full-length laboriosa synthases and homologues from five *Apis* species against
all four retrieved fly reference isoforms longer than 1,300 amino acids, using
the same local BLOSUM62 alignment parameters.

| Laboriosa protein | Best Kkv score | Best Chs2 score | Kkv aligned identity | Chs2 aligned identity |
| --- | ---: | ---: | ---: | ---: |
| XP_043801530.1, elevated group | 5,454 | 2,955 | 66.99% | 50.12% |
| XP_043794468.1, other group | 2,975 | 2,372 | 42.74% | 41.71% |

All tested long isoforms give the same preference. Scores, aligned spans,
coverage, and gapped sequences are retained. This supports describing the
elevated sequence as Kkv-like in this comparison. It does not establish a
gene tree, ancestral duplication history, or actual bee tissue expression.
The other bee protein also favors Kkv, more weakly; it cannot be assigned a
midgut role by elimination.

Functional reference:
[Bertran-Mas et al. 2025](https://doi.org/10.1371/journal.pgen.1011847).
This comparator test and the tissue metadata independently weaken the specific
interpretation of enhanced midgut chitin production.

## 4. Independent GTX-I transcriptome cannot provide the intended validation

[Zhou et al. 2023](https://doi.org/10.7717/peerj.16238) report an RNA-seq
experiment in *Spodoptera litura*. Supplement 2 is labelled as read counts,
but contains 282 selected genes with group means, fold ratios, p-values and
adjusted p-values. It contains neither individual replicate columns nor the
full tested background. The article reports 285 selected genes.

All 282 supplied rows have raw p<0.05; **43** have supplied adjusted p<0.05
and **60** have adjusted p<0.10. We did not recompute multiple-testing correction
within this selected subset. The article also gives incompatible concentration
labels: methods/results say 1.25%, while Figure 4 says 1.25 mg/L.

These restrictions prevent an independent expression reanalysis or unbiased
pathway-enrichment comparison from the supplied file. No public raw-read
accession was located in the retrieved article's data-availability material.
An exposure response in this insect would in any case require separate evidence
before being interpreted as a bee tolerance mechanism.

## What remains worth testing

1. Establish the *A. laboriosa* phenotype under controlled oral exposure, with
   dose, age/task, and viral burden recorded. Pairing oral and injected exposure
   with measured internal GTX would distinguish uptake from systemic handling.
2. Test whether the candidate RNA contrasts recur in independent colonies,
   documented tissues, and matched conditions. Quantify BQCV alongside host
   transcripts. A species effect that disappears after matching these factors
   would contradict a constitutive species-wide explanation.
3. For the specific ABCC and OATP genes, establish tissue localization and
   transport activity before claiming clearance. A perturbation that changes
   internal GTX exposure and the dose response would supply causal evidence.
4. For the chitin-binding candidates, establish protein identity and tissue
   localization. Midgut enrichment would support reconsideration; head/thorax
   structural localization would favor another explanation for these samples.

Additional unreplicated genome or expression scans cannot remove the current
confounding. The most useful computational next dataset would contain
independent colonies, named tissues, treatment/control reads, and sample-level
metadata. The present results justify a short candidate list and correct the
earlier interpretation. A mechanism discovery or literature-first novelty
claim would require evidence beyond this analysis.
