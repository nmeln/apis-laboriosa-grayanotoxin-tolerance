# Population test of three ABCC coding candidates

Started 2026-09-05, before examining population alleles, from main commit
bdfd14ed58a8220effd0fe1ce7182d38cc3df7a2. This is an exploratory follow-up.

The primary question is whether laboriosa residues L254, L549 and F1134 in
XP_043798878.1 consistently distinguish the 29 deposited laboriosa workers
from the 28 dorsata workers in PRJNA931733. The paper describes one randomly
chosen worker per colony. These diploid individuals are the sampling units;
reads and reference assemblies are not additional individuals.

## Predictions and tests

| Hypothesis | Prediction | Test |
| --- | --- | --- |
| Each candidate state is consistent across the sampled laboriosa population | Callable workers carry the laboriosa amino-acid state at that site | Reconstruct codon coordinates, extract reads from all 57 complete paired libraries, and report per-worker diploid calls |
| The state separates the sampled species | Callable dorsata workers carry the other state | Compare amino-acid and codon counts with explicit callable denominators |
| The contrast is caused by reference-dependent read assignment | Genotypes change when the reference changes | Map the retained reads to each species' complete genome separately and require agreement for the primary summary |
| The apparent state is a technical error or insufficient depth | Low quality, duplicates, conflicting mates, strand imbalance or absent coverage drive the call | Preserve fragment counts, quality filters, genotype likelihoods and missing-data reasons; perform prespecified sensitivity checks |
| The method can detect variation within a species | Assembly-variable positions need not appear species-fixed | Include AL positions 452 and 822 as comparison sites, without assuming their allele frequencies in this cohort |

The other four stable substitutions from the previous alignment are secondary
descriptive comparisons. Gap columns are outside the primary SNP/codon test.

## Read collection

- Retrieve the run metadata and individual sample records. Reconcile aliases
  with the paper's actual supplementary sheets; do not guess locality from
  sample numbering. Record unresolved or conflicting metadata.
- Use the complete R1 and R2 FASTQs, validating their ENA byte counts and MD5s.
  Keep an exact record of input reads and retained pairs for every library.
- Extract pairs with a standard k-mer filter using baits from the complete
  candidate loci and flanks in all four already validated assemblies.
- Mask tested codons when constructing capture baits. Do not require a
  particular candidate allele for capture. Retain both mates when either
  matches. Sequence capture is followed by whole-genome mapping.
- Retain all selected read bases, qualities and identifiers. Archive these
  evidence files and the metadata; the full input download is approximately
  270 GB and remains independently retrievable by stable accessions.

## Genotypes and controls

Primary calls require at least ten independent fragments, mapping quality at
least 30, base quality at least 30 at all three codon bases, and genotype
quality at least 30. Exclude secondary/supplementary alignments and marked
duplicates. Overlapping mates count once; conflicting codon observations
within a fragment are excluded. A heterozygous call also requires at least
three supporting fragments per allele and a minor fraction of at least 0.20.

Preserve strand support and flag one-strand-only evidence. A primary result
must agree between the laboriosa and dorsata reference mappings. Repeat with
base quality 20 and minimum fragment depth 8 as a sensitivity analysis; do not
silently substitute a weaker call for a failed primary call.

Validate capture and calling on known synthetic genotypes, including both
alleles, heterozygotes, missing coverage and overlapping mates. These are
method checks and do not count as additional biological evidence.

Report observed codons as well as amino acids, since synonymous variation
can occur at a candidate position. Include unexpected states and failures.
No claim of species-wide fixation, positive selection, GTX transport or
causal tolerance follows from this population comparison. Geographic
sampling and relatedness constrain frequency estimates.
