# Comparative addendum: toxin tolerance across bees

This addendum asks whether public bee genomes contain a recurring molecular pattern that could support the proposed grayanotoxin tolerance mechanism in *Apis laboriosa*.

The scientific frame is comparative evolutionary genomics. More specifically, this is a phylogenetically informed screen for molecular convergence. The available phenotype matrix is too small for a valid formal convergence-association model: one bee lineage has directly measured oral tolerance to grayanotoxin, one has directly measured oral sensitivity, the focal *A. laboriosa* phenotype is inferred from ecology, and three *Apis* controls have no located grayanotoxin assay. The analysis therefore uses explicit exact-state rules and internal controls. It treats hits as leads.

## Result

| Question | Test | Result | Reading |
| --- | --- | --- | --- |
| Do *A. laboriosa* and orally tolerant *Bombus terrestris* share a known Para/Nav resistance pattern? | Six references at 16 experimental grayanotoxin-contact positions | All 16 positions are identical | Known-site target resistance has little sequence support |
| Do they share another exact Para protein state? | 2,050 callable Para residues | Zero strict focal sites | No exact whole-protein Para lead under the stated rule |
| Do they share broad expansions of detoxification, efflux, barrier, or excretion genes? | 10,003 orthogroups | One hit, an odorant-receptor group with a low-quality *A. laboriosa* model | No relevant broad expansion found |
| Are prespecified barrier and toxicokinetic genes enriched for exact shared states? | 4,902,144 callable amino-acid sites with matched permutations | Barrier/detox BH `q = 0.539`; core toxicokinetic BH `q = 0.539` | No enrichment |
| Did any narrow coding lead survive wider taxon checks? | NHE2/3-labelled exchanger orthogroup plus seven more *Bombus* and five other bees | S159 and I232 remain compatible with parallel state changes | Weak lead with no corrected enrichment or direct grayanotoxin evidence |

The cross-insect comparison leaves presystemic handling as the most useful working model: gut uptake, epithelial transport, tissue distribution, metabolism, or excretion may reduce the dose reaching excitable tissue. Public data do not establish that mechanism in *A. laboriosa*.

## Evidence boundary

| Taxon | Grayanotoxin evidence | Role here |
| --- | --- | --- |
| *Apis laboriosa* | Ecological exposure and mad-honey production | Focal hypothesis taxon |
| *Bombus terrestris audax* | Direct controlled oral tolerance | Measured tolerant comparator |
| *Apis mellifera mellifera* | Direct controlled oral sensitivity | Measured sensitive comparator |
| *Andrena carantonica* | Direct sublethal sensitivity | Additional phenotype evidence; no genome in the screen |
| *A. dorsata*, *A. cerana*, *A. florea* | No controlled assay located | Genomic controls with unknown grayanotoxin phenotype |

See [`phenotype_matrix.tsv`](phenotype_matrix.tsv) for the exact claims and sources.

## Files

| Path | Contents |
| --- | --- |
| [`REPORT.md`](REPORT.md) | Methods, results, limits, hypothesis ranking, and experiment design |
| [`REPRODUCE.md`](REPRODUCE.md) | Fresh-clone instructions and expected checks |
| [`AGENTS.md`](AGENTS.md) | Rules for an agent rerunning or extending the addendum |
| [`input_sources.tsv`](input_sources.tsv) | Six added inputs with URLs, accessions, sizes, and SHA-256 hashes |
| [`combined_input_manifest.tsv`](combined_input_manifest.tsv) | The 18 base inputs and six added inputs used by this analysis |
| [`hypothesis_matrix.tsv`](hypothesis_matrix.tsv) | Predictions, tests, results, and decisive follow-ups |
| [`literature_mechanism_table.tsv`](literature_mechanism_table.tsv) | Causal examples from other insect-toxin systems and their limits |
| [`results/`](results/) | Committed tables, alignments, summaries, and tree checks |
| [`scripts/`](scripts/) | Retrieval, orthology, screening, validation, and integrity checks |

## Reproduce

Requirements: Linux x86-64, micromamba, about 2 GB free disk space, and at least four CPU cores. From the repository root:

```bash
make -C comparative_addendum all
```

The command restores the pinned base and addendum input snapshots, verifies every file, creates the pinned bioinformatics environment, reconstructs the six primary proteomes, reruns OrthoFinder and all downstream tests, checks the stated claims directly from the result tables, and compares generated files with the committed SHA-256 manifests.

Large third-party inputs remain outside Git history. The exact added-input archive is published as the GitHub Release asset `comparative-addendum-inputs-v1.tar`; official NCBI and Europe PMC locations are the fallback. See [`REPRODUCE.md`](REPRODUCE.md) and [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md).

## Main sources

- Tiedeken et al. (2016), direct grayanotoxin feeding comparison: <https://doi.org/10.1111/1365-2435.12588>
- Lin et al. (2021), *A. laboriosa* genome: <https://doi.org/10.1093/gbe/evab227>
- Current *B. terrestris* reference `GCF_910591885.1`: <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_910591885.1/>
- Cross-insect mechanism sources are listed in [`literature_mechanism_table.tsv`](literature_mechanism_table.tsv).
