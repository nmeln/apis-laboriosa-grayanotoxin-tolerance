# Independent genomes and a focused transporter test

2026-09-05. Starting commit: `e65a11dbef6f032303a46f8b768da80aa5a4f676`.

## Result

The ABCC transporter in `OG0000499` remains a credible candidate for a direct
functional test. We recovered an intact gene from an independent dorsata
genome, reproduced the RNA contrast using each species' reference separately,
and found independent experimental evidence for the corresponding mellifera
protein in midgut membrane preparations. Three laboriosa coding substitutions
also survive independent genome checks and distinguish the focal reference
from the other five bee references.

These results support two separate hypotheses: the amount of transporter may
differ, and the transporter proteins may function differently. Neither has
been connected experimentally to grayanotoxin (GTX). Higher RNA in the original
two pools is not a replicated species effect. The gut-protein observation is
from another species and cannot establish the location of the original
laboriosa RNA signal.

## Independent dorsata genomes

The NCBI assembly query returned three records representing two independent
assembly sources. The GCA/GCF versions of the Malaysian reference are a paired
GenBank/RefSeq submission of the same assembly.

| Assembly | Source material | Sequencing | Candidate recovered here |
| --- | --- | --- | --- |
| `GCF_000469605.1`, paired with `GCA_000469605.1` | Malaysian Sabah adult male material, collected 2007; `SAMN02954476`, `PRJNA174631` | Illumina | 14 coding exons; 1,333 aa; exact match to the annotated reference protein |
| `GCA_009792835.1`, `RUTG_Adors_1.0` | Pooled drones from Chiang Mai, Thailand, 2017; head and thorax; `SAMN13065485`, `PRJNA577936` | Illumina plus Oxford Nanopore | 14 coding exons; 1,333 aa; one amino-acid difference from the Malaysian reference |

These are whole-genome assemblies assembled into scaffolds, not finished
chromosome-by-chromosome sequences. The Thai assembly has a contig N50 of
30,868 bp and scaffold N50 of 34,588 bp. It is more fragmented at scaffold
level than the Malaysian reference, but the entire candidate is on a single
66,650-bp scaffold, `WJNQ01000516.1`. Its coding model spans positions
14,198 through 20,292, including the terminal stop.

Both reference sources are independent of the Chinese dorsata RNA pool with
abundant BQCV RNA. That RNA observation does not demonstrate a problem with
either DNA assembly. We did not classify either source animal as virus-free.

Sources: [Malaysian assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000469605.1/),
[Thai assembly](https://www.ncbi.nlm.nih.gov/datasets/genome/GCA_009792835.1/),
[Oppenheim et al., genome study](https://doi.org/10.1093/gbe/evz277).
Saved assembly metadata and the full additional genome files are in the input
archive. All assembly records, source samples and contiguity measures are in
[`genome_assembly_inventory.tsv`](results/genome_assembly_inventory.tsv).

## Belly and gut samples

Head and thorax cells contain the DNA for gut genes. They can therefore supply
a whole-genome sequence that includes a gut transporter. RNA and protein
measurements depend on the tissues actually collected.

The [Cao et al. methods](https://doi.org/10.1093/gbe/evad025) explicitly say
abdomens were removed to avoid contamination before DNA extraction. Reducing
microbes and food-derived material in the digestive tract is a plausible
reason for this choice. The paper does not separately explain the RNA
dissection. [GEO](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130963)
records both RNA samples as “whole body without belly”; it does not define the
precise anatomical boundary.

| Data search | What is available | What we can use it for |
| --- | --- | --- |
| Dorsata host RNA, ENA taxonomic subtree plus NCBI SRA | One run, `SRR9034696`, the original pool | No independent focal gut-expression confirmation |
| Laboriosa host RNA, same search | One run, `SRR9034695`, the original pool | Same limitation |
| Giant-bee population sequencing, `PRJNA931733` | 28 dorsata and 29 laboriosa WGS runs | Future population sequence analysis; these are DNA libraries |
| Mellifera midgut membrane proteomics, Arora et al. 2025 | Three adult and three larval preparations, including the exact counterpart of our candidate | Independent gut-associated protein evidence, analyzed below |
| Mellifera and cerana japonica, Yokoi et al. 2025 | Whole-individual RNA-seq across development, including the abdomen | Public follow-up resource; candidate expression was not reanalyzed here |
| Mellifera intestinal single-nucleus atlas, Liu et al. 2025 | Midgut, ileum and rectum study; publisher describes the atlas | No cell-level candidate result here; attempted supplementary downloads returned HTTP 403 |

This inventory describes the saved searches on 2026-09-05. It does not exclude
unindexed, differently labelled, private or future datasets. A gut microbiome
amplicon library measures microbial markers and cannot replace host RNA.

Sources: [focal RNA project](https://www.ebi.ac.uk/ena/browser/view/PRJNA542114),
[population project](https://www.ebi.ac.uk/ena/browser/view/PRJNA931733),
[Yokoi et al.](https://doi.org/10.1038/s41597-025-05279-z),
[Liu et al.](https://doi.org/10.1111/1744-7917.70157).
See [`focal_RNA_inventory.tsv`](results/focal_RNA_inventory.tsv) and
[`related_bee_abdominal_datasets.tsv`](results/related_bee_abdominal_datasets.tsv).

## Tests performed

| Hypothesis | Prediction | Data and computational test | Result and limit |
| --- | --- | --- | --- |
| A deficient dorsata reference creates the candidate | A second assembly lacks the gene or disagrees substantially | Both species' full proteins independently aligned with miniprot 0.18 to four assemblies; translate genomic exon bases | Same 14-exon structure from either query in every assembly; intact proteins. This is reference validation, not population sampling. |
| The RNA contrast is an asymmetric mapping result | It disappears when the reference changes | Each original 2-million-read R1 prefix mapped against all AL primary CDSs and all AD primary CDSs separately | Normalized ratios 17.43 and 17.60; no new biological replication |
| A short repeated segment drives the signal | Reads concentrate in a small portion of the CDS | Unique-group assignment and per-base CIGAR coverage | AL reads cover over 99.8% of either CDS; AD reads cover over 81.8% |
| Identical candidate sequence does not differ in the libraries | A reference-free shared-sequence count removes the contrast | Exact, identical 100-nt markers specific to the gene in both host references and absent from tested BQCV genomes | 258 AL versus 15 AD reads; a technical subset with different detection efficiency from the alignment count |
| The focal protein has reproducible coding differences | Differences survive independent genomes and the expressed transcript | Eight-sequence bee alignment; reconstructed genomic CDSs; six-frame search of the AL assembled transcript | Seven substitutions and three gap columns differ consistently between the two assemblies per species; three substitutions distinguish AL from the other five bee references |
| An abdominal role has independent support | The exact counterpart is measured in a relevant tissue | Extract the matching protein accession and controls from a separate midgut proteomics workbook | Present in all six reported preparations; this concerns mellifera, and no GTX substrate was measured |

The starting plan is in [`ANALYSIS_PLAN.md`](ANALYSIS_PLAN.md). Additional
controls and the untested variant panel are recorded there as exploratory
follow-ups. No phenotype association p-values were computed.

## RNA robustness

Each mapping reference retained every available primary host CDS for that
species, so other genes could compete for reads. Minimap2 used the same
short-read settings as the preceding addendum. Accepted alignments have at
least 100 aligned bases, at least 90% identity, and at least 80% read coverage.
A read is excluded when another orthogroup scores at least 98% of the best.

Normalization centers the AL/AD ratios of one-to-one orthogroups with at
least ten reads in each pool. This is a descriptive host-RNA scaling factor.

| Reference | AL candidate reads | AD candidate reads | Host normalization groups | Normalized AL/AD ratio |
| --- | ---: | ---: | ---: | ---: |
| Laboriosa CDSs | 1,161 | 54 | 6,141 | 17.4324 |
| Dorsata CDSs | 1,158 | 54 | 6,040 | 17.5954 |

With the laboriosa reference, AL reads cover 4,001 of 4,008 coding bases;
AD reads cover 3,280. There are 833 distinct strand/start/end combinations
among the AL candidate reads and 50 among AD. These are not deduplicated
molecule estimates, since these libraries have no unique molecular identifiers.

The aligner-free control uses 1,064 canonical 100-nt sequences shared exactly
by the two candidate CDSs. None occurs in any other pooled primary CDS or
either tested BQCV genome. The central 100 bases of each 150-nt R1 read are
matched in either orientation. Counts are 258 AL and 15 AD, raw ratio 17.2.
This ratio is not the same estimand as the normalized whole-CDS count.

The original pool, dissection, collection and BQCV caveats still apply. These
tests make mapping and assembly artifacts less likely explanations for this
particular RNA contrast. They cannot establish a normal species difference,
toxin induction, or protection.

## Protein differences and a precise coding hypothesis

The eight-bee alignment includes the six original orthogroup representatives
plus the independent Thai dorsata and Pingbian laboriosa sequences. It uses
pyfamsa 0.5.3, one thread and refinement disabled. The laboriosa reference is
99.40% identical to Malaysian dorsata across 1,333 residue-paired columns,
with three additional gap columns. The two laboriosa assemblies differ at
T452/R452. The two dorsata assemblies differ at V819/I819, corresponding to
laboriosa position 822. Thus position 822 does not distinguish these species
consistently even in this small sample.

| AL reference position | AL state in both assemblies | Dorsata state in both assemblies | Other reference bees | Region projected from human ABCC4 |
| --- | --- | --- | --- | --- |
| 254 | L | F | F in mellifera, cerana, florea and terrestris | Third membrane-spanning helix |
| 549 | L | F | F in mellifera and cerana; M in florea and terrestris | First nucleotide-binding domain |
| 1134 | F | L | L in all four | Second nucleotide-binding domain |

These are candidates for altered transport activity or substrate handling.
They are not demonstrated tolerance substitutions, and two assemblies per
species do not establish population fixation. None maps to a currently
annotated human ABCC4 ligand-binding position in this alignment; that does not
rule out an indirect functional effect. Region assignments are inferred from
homology, with human ABCC4 at 47.5% identity across paired columns.

The three-residue length difference also includes a Q447 alignment column
absent from the other tested references. This lies in a variable region;
exact gap placement should not be interpreted as a validated structural change.
All differences, including shared states and the dorsata-variable site, are
retained in [`candidate_coding_differences.tsv`](results/candidate_coding_differences.tsv).

The complete AL transcript `AL|c41541_g1` contains a 4,008-nt ORF encoding all
1,336 reference amino acids exactly. Its sequence, orientation and ORF
coordinates are recorded. The AD raw reads demonstrate the gene's presence
in the RNA pool despite its failure to appear as a complete assembled unigene
in the earlier screen.

A family-reference comparison finds human ABCC4 closest among the tested
long reviewed human ABCC proteins. Several fly proteins are closer than
the fly gene named lethal(2)03659; CG7627 has the highest score in the retrieved
ABCC-domain panel. This is a similarity screen, not a resolved gene tree.
Specific fly Mrp4 or lethal(2)03659 functions must not be transferred solely
from the bee annotation name.

Source for projected features: [UniProt O15439](https://www.uniprot.org/uniprotkb/O15439/entry).
The full source records, comparison panel, alignment and mapped features are
preserved.

## Independent gut-protein evidence

[Arora, Mishra and Bonning, 2025](https://doi.org/10.1038/s41598-025-26662-1)
identified `XP_393750.5`, the exact mellifera member of our one-to-one bee
orthogroup, in midgut brush border membrane vesicle preparations. The authors
used mass spectrometry and a RefSeq database search, with a reported 1%
peptide-level false discovery rate and protein-retention criteria of at least
two unique peptides and two peptide spectral matches.

| Stage | Replicate 1 spectral count | Replicate 2 | Replicate 3 | Unique peptides listed in the stage summary |
| --- | ---: | ---: | ---: | ---: |
| Adult | 66 | 52 | 57 | 31 |
| Larva | 7 | 39 | 67 | 31 |

The adults were newly emerged, less than 24 hours old and unfed. Thus the
reported protein presence does not require prior adult foraging in that
experiment. The stage-level peptide totals are not 31 independently unique
peptides in each of the six replicates.

We extracted these observations from the complete supplementary workbook and
checked the orthogroup accession. Two other ABCC proteins are included as
comparison records. We did not re-search raw spectra: the public Figshare
peptide exports repeatedly redirected to expired download signatures. Their
file IDs, expected sizes, MD5 values and source URLs are preserved for later
retrieval. The analyzed workbook is archived in full.

The adult hybrid-normalized values disagree between Tables S4 and S6. We
preserve that discrepancy and use the reported spectral counts for the
presence conclusion. We make no new adult/larval abundance or significance
claim. This disagreement does not erase the candidate's reported detection
in all preparations.

Gut association was already published for mellifera. The connection to our
laboriosa candidate is the orthogroup match combined with the new genome,
RNA and coding checks. It does not demonstrate an exclusively gut-localized
protein, a laboriosa gut increase, transport direction, or GTX transport.

## What would discriminate the remaining hypotheses

| Hypothesis | Measurement that would support it | Outcome that would weaken it |
| --- | --- | --- |
| The protein transports GTX or a GTX conjugate | Direct, ATP-dependent transport measured chemically in a system expressing the candidate, with appropriate empty-system and activity controls | No transport despite verified protein localization and transport competence, within the tested chemical and concentration range |
| Laboriosa coding changes improve that transport | Different transport kinetics at matched protein abundance; a specified substitution changes the result | Comparable activity of the species proteins and tested variants under the same conditions |
| Regulation changes internal exposure | Replicated gut, excretory-tissue or neural-barrier protein measurements together with tissue GTX time courses | The RNA difference disappears in matched biological samples or does not track protein/internal exposure |
| The observed difference concerns another physiological role | Activity on another substrate or expression changes associated with age, site or viral burden, without a GTX effect | A replicated GTX-specific functional result persists after these alternatives are controlled |

[`untested_variant_panel.fna`](results/untested_variant_panel.fna) and the
matching protein FASTA specify AL and AD reference sequences, three single
AL-to-AD substitutions (L254F, L549F, F1134L), and their combination. Each
replacement uses the naturally occurring dorsata codon. These are sequence
designs only; no synthesis, expression or assay was performed. The CDS files
exclude the terminal stop codon, which a downstream expression design must
handle explicitly. No expression host or tagging arrangement is prescribed.

GTX tolerance in laboriosa itself still needs a controlled phenotype and
internal-exposure measurement. The existing data cannot answer the transport
question through docking scores, extra normalization or more technical
subsets of the same RNA pool.

## Reproduction

Run `make -C candidate_followup all`. The archive preserves 28 added files,
including the full two downloaded dorsata genomes, current Bombus protein
reference, complete proteomics workbook, source articles, sequence annotations,
saved archive queries and retained exploratory records. Reproduction also
uses the earlier checked base inputs and two fixed R1 prefixes.

All result tables, recovered coding sequences, protein alignments, untested
variant definitions and validation rules are generated by the scripts in this
directory. See [REPRODUCE.md](REPRODUCE.md) for exact versions and commands.
