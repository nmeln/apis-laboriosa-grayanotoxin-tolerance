# ABCC candidate genotypes in 57 giant-honeybee workers

Completed 2026-09-05. Analysis developed and run with Codex using published data.

## Result

The three candidate coding differences recur across the sampled populations.
Every accepted laboriosa genotype has L254, L549 and F1134. Every accepted
dorsata genotype has the contrasting F, F and L states. The sequence contrast
is supported beyond the four reference assemblies used in the preceding work.

| Site in XP_043798878.1 | Laboriosa diploid codon | Workers supporting it | Dorsata diploid codon | Workers supporting it |
| --- | --- | --- | --- | --- |
| 254 | CTT/CTT, L/L | 28 callable / 29 deposited | TTT/TTT, F/F | 27 callable / 28 deposited |
| 549 | TTG/TTG, L/L | 28 / 29 | TTT/TTT, F/F | 26 / 28 |
| 1134 | TTT/TTT, F/F | 28 / 29 | CTT/CTT, L/L | 28 / 28 |

L is leucine and F is phenylalanine. Dorsata protein positions corresponding
to these rows are 254, 546 and 1131. All numbering and reverse-strand genomic
coordinates are in [target_coordinates.tsv](results/target_coordinates.tsv).

Both reference mappings agree on every best codon genotype at all nine sites,
including calls excluded by quality rules. Accepted primary-site calls cover
165 of 171 worker/site combinations. Fifty-three workers, 28 laboriosa and 25
dorsata, are callable at all three sites. All five laboriosa and all six dorsata
archive localities have callable evidence for each candidate position. The
three-site profiles are unphased; distant alleles were not reconstructed onto
individual chromosomes.

The prespecified sensitivity analysis lowers base quality from 30 to 20 and
minimum depth from ten to eight fragments. It leaves primary-site genotypes
and callable denominators unchanged. At the six other tested positions it
adds two accepted laboriosa calls, at 452 and 727.

Sources: [species counts](results/species_site_summary.tsv),
[locality counts](results/location_site_summary.tsv),
[individual profiles](results/individual_primary_profiles.tsv),
[all calls](results/consensus_genotypes.tsv).

## Missing calls

These exclusions follow the rules recorded before cohort inspection. Depth
counts independent, nonduplicate fragments spanning the complete codon with
the required base and mapping qualities.

| Worker | Location | Site | Qualifying evidence | Reason excluded |
| --- | --- | --- | --- | --- |
| SRR23343455, HD17, laboriosa | Diqing | 254, 549, 1134 | Depths 1, 7 and 2 | Insufficient depth and genotype likelihood separation |
| SRR23343440, D12, dorsata | Puer | 549 | 15 TTT fragments, all in one orientation | One-strand-only evidence |
| SRR23343460, D16, dorsata | Xishuangbanna | 549 | 17 TTT and one TTG fragment | Genotype-quality score 19.392, below 30 |
| SRR23343476, D31, dorsata | Chongzuo | 254 | 15 TTT and one TTC fragment | Genotype-quality score 13.375, below 30 |

The single TTG observation in D16 encodes the laboriosa leucine state. Its
1/18 support does not meet the diploid heterozygote rules. The sensitivity
analysis has 20 TTT and the same single TTG fragment, and remains excluded.
This evidence cannot distinguish a sequencing or library artifact from a real
minor sequence contribution. It is not counted as an accepted cross-species
allele or silently discarded. D31's minority TTC codon is synonymous with TTT;
both encode phenylalanine. That codon call also remains missing.

HD17 has reduced captured coverage across the gene and surrounding sequence.
Against the laboriosa reference, mean fragment depths are 3.71 in coding
sequence, 4.33 in introns and 6.83/9.57 in the two 1 kb flanks. Its complete
source files passed the same checks as other workers. Coverage at this single
captured locus cannot identify the cause or establish a gene deletion. Its
low-depth apparent R/T mixture at comparison position 452 is excluded too.

Sources: [codon evidence and filters](results/sample_codon_calls.tsv),
[fragment observations](results/fragment_evidence.tsv),
[coverage summaries](results/locus_coverage_summary.tsv),
[100-base coverage bins](results/locus_coverage_100bp.tsv).

## Variation at comparison positions

The comparison sites were included because the earlier assemblies disagreed
within a species. They received the same capture, mapping and calling rules
as the primary candidates.

| Site | Species | Accepted amino-acid genotypes | Missing workers |
| --- | --- | --- | --- |
| 452 | Laboriosa | 24 T/T, two R/T, one R/R | 2 of 29 |
| 452 | Dorsata | 28 T/T | 0 of 28 |
| 822 | Laboriosa | 28 I/I | 1 of 29 |
| 822 | Dorsata | 25 I/I, three I/V | 0 of 28 |

T is threonine, R arginine, I isoleucine and V valine. The laboriosa R-bearing
workers are HD16 in Diqing and HD2/HD3 in Baoshan. The dorsata I/V workers are
D19 in Baoting and D22/D6 in Qiongzhong, all on Hainan. The sampled dorsata
V state appears in three heterozygotes and no accepted V/V individual. It was
the state in the Malaysian reference; the Thai reference has I.

These calls demonstrate recovery of within-species variation and heterozygous
sequence evidence. Geographic association, selective value and toxin effects
were not tested. The four secondary positions also retain their reference
contrasts wherever callable: laboriosa Q474/A717/A727/L766 versus dorsata
K/T/V/F. Their full denominators are in the species table.

## Predictions tested

| Prediction | Empirical result | Consequence |
| --- | --- | --- |
| The three candidate states recur across laboriosa workers | 28 accepted homozygotes for the expected state at each position | Supports cohort consistency beyond reference assemblies |
| Dorsata workers carry the contrasting states | All 27, 26 and 28 accepted site-specific genotypes have the contrasting state | Supports separation in the callable sampled workers |
| Reference-dependent mapping explains the contrast | No best-codon disagreement after swapping complete genomes | This explanation loses support in the tested evidence |
| The method can recover within-species variation | Both comparison sites vary, with accepted heterozygotes | A uniformly reference-copying result is inconsistent with these controls |
| These substitutions alter GTX transport | No transport measurement available | Unresolved |

Neutral divergence, selection on another substrate and linkage to another
causal change remain possible. All three primary sites were selected from an
earlier sequence comparison; they are not an unbiased genome-wide sample. No
selection statistic or toxin-associated significance test is inferred from
their species separation.

## Data and method

[Cao et al. (2023)](https://doi.org/10.1093/gbe/evad025) describe one randomly
selected diploid worker per colony, with abdomens removed for DNA extraction.
The deposited cohort contains 29 laboriosa and 28 dorsata libraries in
[PRJNA931733](https://www.ebi.ac.uk/ena/browser/view/PRJNA931733). Its 114 complete
FASTQs contain 577,920,276,600 sequenced bases in 269,858,255,055 compressed
bytes. Each complete file passed its ENA size and MD5 check before collection.

Four validated candidate loci plus 1 kb flanks supplied the sequence baits.
All nine tested codons were masked in all four assemblies. The 8,602 canonical
21-base baits exclude low-complexity sequences. BBDuk 39.91 selected a pair
when either mate matched. A second pass copied exact four-line records from
the verified originals because BBDuk changed quality scores at N bases in a
synthetic test. The 86,296 retained pairs occupy 58,049,708 FASTQ bytes.

Minimap2 2.31 maps the same retained pairs to each complete species genome,
preserving alternative whole-genome mapping competitors. SAMtools 1.21,
through pysam 0.23.3, marks duplicates within each worker's read group.
The codon caller excludes duplicates, nonprimary alignments, failed reads and
MAPQ below 30. Both eligible mates count as one fragment; conflicting mates
are excluded. Primary calls require all three base qualities at least 30,
at least ten fragments, support in both read orientations and genotype-quality
score at least 30. Heterozygotes additionally require three fragments per
allele and minor support of at least 20%.

The caller considers every unordered pair of the 64 codons. Its quality score
is the best-versus-second-best log-likelihood separation, capped at 99, using
equal allele contribution for diploid heterozygotes and no population prior.
It is not a calibrated posterior probability. Exact rules and their timing
are in [ANALYSIS_PLAN.md](ANALYSIS_PLAN.md).

BCFtools 1.21 independently calls each worker at every tested codon base,
including reference positions. All 2,982 comparisons for accepted codon calls
agree, including 990 at primary candidates. The other 96 comparisons belong
to excluded codon calls. This cross-check shares the mapped evidence and does
not independently validate sample identity, extraction or protein function.
See [the comparison table](results/bcftools_base_crosscheck.tsv),
[BCFtools documentation](https://samtools.github.io/bcftools/bcftools.html) and
[SAMtools duplicate-marking documentation](https://www.htslib.org/doc/samtools-markdup.html).

## Controls and preservation

- All 340,992 synthetic error-free 150-base read starts spanning the 36
  assembly/site combinations and 64 possible codons are captured.
- The actual BBDuk test retains 6,912 positive pairs and excludes 100 negative
  pairs. It checks both orientations, unmatched mates and exact record copying.
- Eleven known-truth genotype tests cover homozygotes, heterozygotes, unexpected
  codons, missing evidence, reverse-strand coordinates, overlapping mates,
  conflicting mates and duplicate marking within versus between workers.
- Independent local full-source collections of SRR23343443, SRR23343450 and
  SRR23343494 reproduce the archived selected FASTQs byte for byte.
- The two remote batches independently collected SRR23343446 and produced
  identical selected records and collection manifests. The merged cohort has
  all 57 expected accessions and no missing libraries.

The exhaustive capture check covers the tested codon changes on reference
backgrounds. Unrepresented surrounding variation, structural changes and
sequencing errors could still reduce capture. Missing calls are retained, and
the full-source collection remains separately reproducible.

The [input release](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/candidate-inputs-v1)
preserves 296 files in a 337,141,760-byte archive. It contains complete genomes,
the selected original paired reads, metadata, article XML, supplementary
workbook and collection records. It excludes the approximately 270 GB of
original population FASTQs; their complete URLs, sizes, MD5s and collection
scripts are preserved. [Source credit and terms](THIRD_PARTY_DATA.md) accompany
the archive. The original and accelerated collection runs are
[33993933619](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/actions/runs/33993933619)
and [33994894621](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/actions/runs/33994894621).

## Source disagreements

| Item | Published material | Deposited records | Treatment here |
| --- | --- | --- | --- |
| Tibetan locality totals | Table S2: Rikaze 5, Linzhi 5 | BioSamples: Rikaze 6, Linzhi 4 | Preserve both; use archive locality labels in generated tables |
| Mapping-table membership | Table S3 includes 60 aliases | 57 matching deposited runs | Exclude unavailable D2, D30 and HD192 |
| Sequencer | Article: HiSeq 2000 | ENA: NovaSeq 6000 | Record disagreement |
| Raw sequence amount | Article: approximately 538 Gb | ENA: 577.9202766 Gb | Use verified deposited file totals |
| Tissue | Methods: abdomen removed | BioSample: whole body | Preserve conflict; no tissue-expression claim |

This is the defined 57-worker cohort, not a census of all subsequent public
bee sequencing. Chinese population structure, finite sampling, unmeasured
relatedness and the absence of Nepalese workers constrain generalization.

## What this changes

The three substitutions remain sensible sequence changes to test in the
specific ABCC protein. Their consistency in callable population data weakens
the explanation that they arose solely from an unusual reference individual
or an isolated assembly error. The present results cannot rank their effects
or establish that this protein carries grayanotoxin. Matched-abundance direct
transport measurements would test that biochemical claim. The species'
grayanotoxin phenotype and internal toxin exposure also remain unresolved.
