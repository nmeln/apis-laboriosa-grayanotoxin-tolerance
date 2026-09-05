# Population test of the ABCC candidate

The three candidate protein differences survive testing in public population
reads. Every accepted genotype carries its species' expected state. This makes
an isolated reference-genome error or an unusual reference individual a less
likely explanation. The effect on grayanotoxin transport remains unmeasured.

| Laboriosa protein position | Laboriosa result, 29 workers deposited | Dorsata result, 28 workers deposited |
| --- | --- | --- |
| 254 | 28/28 callable workers: leucine/leucine; 1 missing | 27/27: phenylalanine/phenylalanine; 1 missing |
| 549 | 28/28: leucine/leucine; 1 missing | 26/26: phenylalanine/phenylalanine; 2 missing |
| 1134 | 28/28: phenylalanine/phenylalanine; 1 missing | 28/28: leucine/leucine; 0 missing |

Calls must pass the stated quality rules against both complete species genomes.
The prespecified lower-quality sensitivity analysis gives the same primary-site
counts. Fifty-three workers have accepted calls at all three positions.

The two comparison positions vary within species. At position 452, laboriosa
has 24 T/T, two R/T and one R/R accepted genotypes. At position 822, dorsata has
25 I/I and three I/V. All three I/V workers come from Hainan. These observations
show that the method can recover heterozygotes and within-species variation.
They do not establish a geographical or toxin-related association.

The four workers with missing primary calls are listed in the
[report](REPORT.md#missing-calls). One excluded dorsata site contains a single
fragment with the laboriosa codon among 18 qualifying fragments. It remains
unresolved under the original rules. Missing calls are never filled from the
reference sequence.

## Evidence and checks

- All 114 complete FASTQs passed their original byte-count and MD5 checks:
  269,858,255,055 compressed bytes and 577,920,276,600 sequenced bases.
- The preserved evidence contains 86,296 original paired records. All tested
  codons were masked when constructing the capture baits.
- No accepted genotype changes between the two complete-genome mappings.
  BCFtools agrees at all 2,982 accepted codon-base comparisons, including 990
  at the three primary sites. It uses the same mapped evidence.
- Synthetic controls cover all possible target codons, known diploid genotypes,
  unmatched mates, overlapping mates and duplicate handling. Independent local
  collection of three libraries reproduces the selected reads byte for byte.

## Reproduce

```sh
git clone https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance.git
cd apis-laboriosa-grayanotoxin-tolerance
make -C population_addendum all
```

The [input release](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/population-inputs-v1)
preserves 296 analyzed files, including all four complete assemblies, selected
original reads, sample metadata and collection records. The archive is about
337 MB. Full-source recollection remains separately runnable.

| File | Contents |
| --- | --- |
| [REPORT.md](REPORT.md) | Findings, missing evidence, source disagreements and scientific limits |
| [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md) | Predictions and rules recorded before cohort genotype inspection, with later checks identified |
| [REPRODUCE.md](REPRODUCE.md) | Installation, complete-source collection, analysis and verification commands |
| [AGENTS.md](AGENTS.md) | Instructions for an agent rerunning or extending the work |
| [species_site_summary.tsv](results/species_site_summary.tsv) | Counts at all nine sites and both quality thresholds |
| [location_site_summary.tsv](results/location_site_summary.tsv) | Counts separated by collection locality |
| [consensus_genotypes.tsv](results/consensus_genotypes.tsv) | Every individual call and exclusion |
| [fragment_evidence.tsv](results/fragment_evidence.tsv) | Read-pair observations behind the calls |
| [input_sources.tsv](input_sources.tsv) | Source links, sizes and SHA-256 hashes |
| [THIRD_PARTY_DATA.md](THIRD_PARTY_DATA.md) | Original investigators, citations and data terms |

The candidate is OG0000499, laboriosa LOC122718161 / XP_043798878.1, described
in the [preceding follow-up](../candidate_followup/REPORT.md). Population reads
come from [Cao et al. (2023)](https://doi.org/10.1093/gbe/evad025),
[PRJNA931733](https://www.ebi.ac.uk/ena/browser/view/PRJNA931733).

This analysis was developed and run with Codex. It reuses published sequencing
data and performs no laboratory experiment. These Chinese workers have no
matched grayanotoxin phenotype. Cohort consistency does not establish
species-wide fixation, positive selection, Nepalese allele states or protection.
