# Comparative addendum: cross-insect screen for grayanotoxin-tolerance mechanisms

Date: 2026-09-01

Repository: [nmeln/apis-laboriosa-grayanotoxin-tolerance](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance)

Status: reproducible addendum to the primary repository analysis

## Bottom line

The comparison is comparative evolutionary genomics. The narrower form is a phylogenetically informed screen for molecular convergence associated with toxin tolerance.

The current data support four conclusions.

1. A shared grayanotoxin-resistant sodium-channel sequence is now a weak explanation. All 16 experimentally implicated Para/Nav residues are identical across the six bee references. A separate whole-protein test found zero exact *Apis laboriosa* plus *Bombus terrestris* states absent from the other four *Apis* species across 2,050 callable Para residues.
2. Broad expansion of detoxification, efflux, gut-barrier, or excretion gene families is also weak. The only exact shared copy-number hit among 10,003 orthogroups was an odorant-receptor group containing one low-quality *A. laboriosa* model.
3. The genome-wide screen found no enrichment of exact *A. laboriosa* plus *B. terrestris* amino-acid sharing in the prespecified barrier, detoxification, or toxicokinetic gene sets. The observed densities were slightly below the matched genomic background.
4. One NHE2/3-labelled sodium/hydrogen-exchanger orthogroup is worth retaining as a narrow lead. It contains four exact *A. laboriosa* plus *B. terrestris* states. Two, S159 and I232 in the *A. laboriosa* protein, remain compatible with independent state changes after comparison with seven more *Bombus* proteins and five non-*Bombus* bees. The signal is not significant after category-wide correction, has no known connection to grayanotoxin, and is not accompanied by higher constitutive whole-worker RNA abundance in *A. laboriosa*.

The best current working model remains presystemic toxicokinetics: gut handling, epithelial transport, tissue localization, metabolism, or excretion that limits the grayanotoxin dose reaching excitable tissue. This is a testable model supported by causal examples involving other insect toxins. It has not been demonstrated in *A. laboriosa*.

## What the phenotype data actually establish

| Taxon | Status used here | Basis | Valid role | Main limit |
|---|---|---|---|---|
| *Apis laboriosa* | Tolerance suspected | Ecological exposure and mad-honey production | Focal hypothesis taxon | No controlled oral challenge, internal-dose measurement, or injection test was located |
| *Bombus terrestris audax* | Oral tolerance measured | Controlled non-choice feeding | Same-toxin tolerant comparator | Mechanism and internal exposure were not measured; the reference genome is not from the assayed population |
| *Apis mellifera mellifera* | Oral sensitivity measured | Controlled non-choice feeding | Same-toxin sensitive comparator | One subspecies and experimental population |
| *Andrena carantonica* | Sublethal sensitivity measured | Controlled non-choice feeding | Additional phenotype evidence | No clean tolerant phenotype and no matched genome in this analysis |
| *A. dorsata*, *A. cerana*, *A. florea* | Unknown for grayanotoxin | No controlled assay located | Genomic controls | They must not be labelled grayanotoxin-susceptible |

The available evidence does not justify calling *A. dorsata*, *A. cerana*, or *A. florea* grayanotoxin-susceptible. They are genomic controls with unknown grayanotoxin phenotypes.

## Why a formal phenotype-convergence model was not run

PCOC, RERconverge, CSUBST, and related convergence methods need several independent foreground lineages with defensible phenotypes and enough branches for calibration. The same-toxin matrix available here contains one directly measured tolerant lineage, one directly measured sensitive lineage, one ecologically inferred focal lineage, and three unphenotyped *Apis* controls. A branch-site model could test selection on one named lineage, but it would not establish association with grayanotoxin tolerance under this phenotype design.

That design cannot separate a grayanotoxin-associated signal from lineage history with a formal phenotype-convergence score. A larger model would still return numbers, but those numbers would not constitute a valid association test. The analysis therefore used a transparent exact-state screen, ran the identical rule on every other *Apis* species as an internal control, and treated the output as candidate generation.

## Data and tests

### Inputs

- Five RefSeq *Apis* genomes and annotations from the base repository
- Current *Bombus terrestris* reference `GCF_910591885.1`, assembly `iyBomTerr1.2`
- One longest protein per GFF gene: 9,783 to 10,310 proteins per species
- Older *B. terrestris* reference proteins for an independent sequence check
- Seven additional *Bombus* exchanger proteins and five non-*Bombus* bee exchanger proteins
- Existing pooled, untreated, whole-worker *A. laboriosa* and *A. dorsata* transcript assemblies and abundance tables

Every added URL, accession, byte count, and SHA-256 digest is recorded in [`input_sources.tsv`](input_sources.tsv). [`combined_input_manifest.tsv`](combined_input_manifest.tsv) records the base and addendum inputs used by this workflow.

### Orthology and strict-state rule

OrthoFinder 3.1.5 assigned 59,381 of 60,216 proteins to 10,003 orthogroups. The final set contained 8,683 orthogroups represented in all six species and 8,273 complete single-copy orthogroups.

For each callable aligned amino-acid column, the strict focal pattern required:

- *A. laboriosa* and *B. terrestris* to carry the same residue;
- all four other *Apis* references to carry one common residue;
- the focal residue to differ from the other-*Apis* residue.

The same rule was repeated with *A. dorsata*, *A. mellifera*, *A. cerana*, and *A. florea* as the focal *Apis*. A strict match can reflect parallel change, reversal, retained ancestral state, alignment error, annotation error, or chance. It is not a phenotype association by itself.

### Main results

| Test | Scale | Result | Interpretation |
|---|---:|---:|---|
| Known Para/Nav functional sites | 16 residues | 16 identical across all six references | Known-site target resistance disfavored |
| Whole Para strict-state screen | 2,050 callable residues | 0 focal strict sites | Exact shared Para coding mechanism disfavored |
| Complete single-copy genome screen | 8,273 orthogroups; 4,902,144 callable sites | 2,024 focal strict sites in 1,311 orthogroups | Candidate pool only; count is close to the *A. dorsata* control, 2,036 |
| Barrier or detox set | 244 orthogroups | 58 strict sites; matched density difference -0.0618 per 1,000 sites; permutation `p = 0.446`, BH `q = 0.539` | No enrichment |
| Core toxicokinetic set | 118 orthogroups | 32 strict sites; matched density difference -0.0690 per 1,000 sites; permutation `p = 0.485`, BH `q = 0.539` | No enrichment |
| Shared copy-number gain | 10,003 orthogroups | 1 hit, `OG0000196`, odorant receptor Or1-like | No shared detox, efflux, barrier, or excretion expansion |
| NHE2/3-labelled category | 6 orthogroups | 4 sites in one orthogroup; raw Fisher `p = 0.0759`, BH `q = 0.365` across 130 exploratory category/focal tests | Narrow lead without category-level support |

The matched permutation shuffled candidate labels within quartile bins of callable protein length and focal-to-*Bombus* divergence. It used 20,000 permutations with a fixed seed. This controls two obvious sources of excess hits, although it does not replace a phenotype-aware phylogenetic test.

The annotation-based family counts are rough screens. Several annotation terms are incomplete or inconsistent across species. Their useful result is the absence of an obvious *A. laboriosa* plus *B. terrestris* expansion, not an exact census of each biochemical family.

## The exchanger lead

The protein is labelled sodium/hydrogen exchanger 3 in four references and sodium/hydrogen exchanger 2 in two. The orthology result supports a shared exchanger group, but the NHE2 versus NHE3 name is not secure. It is therefore called the **NHE2/3-labelled exchanger orthogroup** here.

| *A. laboriosa* position | Focal state | Other four *Apis* | Additional *Bombus* matching focal | Five non-*Bombus* bees matching focal | Readout |
|---:|:---:|:---:|---:|---:|---|
| 44 | N | D | 3 of 8 | 0 of 5 | Variable in *Bombus*; weak |
| 159 | S | P | 7 of 8 | 0 of 5 | Candidate parallel state |
| 232 | I | V | 8 of 8 | 0 of 5 | Candidate parallel state |
| 353 | A | G | 8 of 8 | 5 of 5 | Broadly distributed state; weak as convergence evidence |

The eight additional *Bombus* checks include the older *B. terrestris* reference. At position 232, four of five external bees carry V and one carries L. Both *A. laboriosa* and *A. dorsata* assembled worker transcripts cover all four positions and match their reference proteins. This reduces the chance that the four focal residues are simple reference-annotation errors.

The existing whole-worker RNA data do not show constitutive elevation of this exchanger:

| Measure | *A. laboriosa* | *A. dorsata* | Ratio |
|---|---:|---:|---:|
| FPKM | 20.54 | 24.97 | 0.823 |
| CPM | 93.35 | 107.32 | 0.870 |

These values come from one pooled untreated whole-worker library per species. They cannot test tissue-specific expression, inducibility, isoform use, protein abundance, or transport activity.

Insect NHE proteins participate in epithelial ion and pH handling in the gut and Malpighian tubules, which gives this lead a defensible physiological context. No paper linking this exchanger to grayanotoxin transport or resistance was found. It may reflect unrelated osmoregulatory ecology, altitude adaptation, lineage history, or chance. It should be tested only in a design that also measures grayanotoxin pharmacokinetics.

## What other insect toxins tell us

Different toxins cannot confirm a bee grayanotoxin mechanism. They do establish useful priors about mechanisms that evolution repeatedly uses.

| System | Causal result | Mechanism class | Relevance here |
|---|---|---|---|
| Cardenolide-feeding insects | Repeated substitutions in Na,K-ATPase, including a functionally tested resistance change | Target-site insensitivity | Positive control showing that recurrent target-site adaptation can be visible across distant insects |
| *Pieris rapae* on Brassicales | Gut nitrile-specifier protein redirects plant-defense chemistry toward less-toxic nitriles that are excreted | Presystemic gut biotransformation | Direct example of protecting internal tissues before toxin exposure |
| *Manduca sexta* on tobacco | Midgut CYP6B46 knockdown changes nicotine transfer from gut to hemolymph and defensive exhalation | Tissue-specific toxicokinetics | A single regulated gut enzyme can alter systemic exposure without family expansion |
| *Drosophila sechellia* on toxic noni | RNAi implicates Osiris genes; important loci include tissue-specific regulation without fixed derived protein changes | Barrier or trafficking regulation | Protein sequence and gene-count screens can miss a real adaptation |
| *Apis cerana* and triptolide nectar | Controlled feeding shows greater tolerance than *A. mellifera* | Comparative toxic-nectar phenotype | A close pollinator precedent; the toxin and unresolved mechanism differ |

These examples strengthen the case for measuring internal dose, tissue distribution, and inducible regulation. They do not validate any particular *A. laboriosa* gene.

## Current hypothesis ranking

| Hypothesis | Current status | Reason |
|---|---|---|
| Known-site Para/Nav resistance | Low | All 16 implicated sites are identical |
| Other exact shared Para coding change | Low | Zero strict focal sites across the protein |
| Broad detox, efflux, or barrier expansion | Low | No broad count pattern and no relevant strict copy-number hit |
| Broad shared toxicokinetic coding convergence | Low | Prespecified sets are not enriched |
| NHE2/3-labelled exchanger coding effect | Interesting, weak | Two sequence states survive wider taxon checks; no corrected enrichment, functional evidence, or GTX link |
| Tissue-specific regulation or localization | Open | Existing whole-worker, untreated, unreplicated RNA cannot test it |
| Presystemic exclusion or rapid clearance | Highest-priority working model | Fits direct oral *Bombus* tolerance and causal cross-toxin precedents; no *A. laboriosa* pharmacokinetic data yet |
| Behavioral dose regulation | Open | Could contribute in *A. laboriosa*; cannot explain the non-choice *Bombus* result alone |
| Microbiome-mediated metabolism | Open | Current data do not test it |

## Decisive next experiment

The highest-information study is a species-by-route-by-dose experiment.

### Species

- *A. laboriosa*, the focal taxon
- *B. terrestris audax*, a directly measured orally tolerant comparator
- *A. mellifera mellifera*, a directly measured sensitive comparator
- *A. dorsata*, the closest genomic control

### Exposure

- purified GTX-I and GTX-III, alone and at the measured nectar ratio;
- matched vehicle controls;
- oral exposure and hemocoel injection;
- several doses spanning sublethal to lethal exposure;
- randomized allocation, blinded scoring, preregistered exclusions, and biological replication across colonies and collection sites.

### Measurements

- survival, motor function, feeding, and actual consumed dose;
- isotope-labelled or analytically validated LC-MS/MS time courses in crop, midgut, hemolymph, brain, Malpighian tubules, and frass;
- replicated tissue-resolved RNA-seq or proteomics at baseline and after exposure;
- Para electrophysiology if injection sensitivity differs among species;
- exchanger transport, pH, localization, and allele tests only if the pharmacokinetic data point toward epithelial handling.

The key discriminator is simple. Oral tolerance with injection sensitivity and low hemolymph or brain exposure supports presystemic restriction or clearance. Tolerance to both routes with similar internal exposure moves the priority toward target-tissue physiology, channel function, or downstream neural protection.

## Limits

- Controlled grayanotoxin tolerance has not been established in *A. laboriosa*.
- The genomes are species references and are not from the exact populations used in the feeding study.
- Six species with one measured tolerant lineage cannot support a formal phenotype-convergence scan.
- Exact residue sharing does not establish independent substitution, phenotype association, or causal function.
- Protein annotations and family counts contain false negatives, false positives, and inconsistent names.
- The transcript comparison has one pooled untreated whole-worker library per species and no biological replication.
- Ecological ingestion does not establish the grayanotoxin concentration reaching hemolymph, brain, or sodium channels.

## Reproduction and audit files

- [`README.md`](README.md): concise results and evidence boundary
- [`REPRODUCE.md`](REPRODUCE.md): environment, commands, expected checks, and interpretation boundary
- [`environment.yml`](environment.yml): pinned software environment
- [`run_analysis.sh`](run_analysis.sh): full analysis entry point
- [`AGENTS.md`](AGENTS.md): workflow and scientific guardrails for future agents
- [`scripts/`](scripts/): source code for fetching, orthology preparation, screens, validation, and manifests
- [`phenotype_matrix.tsv`](phenotype_matrix.tsv): measured and inferred phenotype boundary
- [`hypothesis_matrix.tsv`](hypothesis_matrix.tsv): predictions, tests, readouts, and decisive follow-ups
- [`literature_mechanism_table.tsv`](literature_mechanism_table.tsv): cross-insect evidence and limits
- [`results/`](results/): complete result tables, alignments, summaries, and tree quality-control files
- [`input_sources.tsv`](input_sources.tsv): added inputs with URLs and SHA-256 digests
- [`primary_proteomes.sha256`](primary_proteomes.sha256), [`orthofinder_key_outputs.sha256`](orthofinder_key_outputs.sha256), [`results.sha256`](results.sha256), [`scripts.sha256`](scripts.sha256): integrity checks

## Sources

- Tiedeken et al. 2016, direct grayanotoxin feeding comparison: <https://doi.org/10.1111/1365-2435.12588>
- Dobler et al. 2012, convergent cardenolide target-site resistance: <https://doi.org/10.1073/pnas.1202111109>
- Wittstock et al. 2004, *Pieris rapae* gut detoxification: <https://doi.org/10.1073/pnas.0308007101>
- Kumar et al. 2014, *Manduca sexta* nicotine toxicokinetics: <https://doi.org/10.1073/pnas.1314848111>
- Andrade Lopez et al. 2017, *Drosophila sechellia* Osiris loci: <https://doi.org/10.1111/mec.14001>
- Wang et al. 2022, *Apis cerana* and triptolide-containing nectar: <https://doi.org/10.1016/j.jinsphys.2022.104358>
- Pullikuth et al. 2006, insect NHE3 epithelial physiology: <https://doi.org/10.1242/jeb.02419>
- Rey et al. 2018, PCOC: <https://doi.org/10.1093/molbev/msy114>
- Kowalczyk et al. 2019, RERconverge: <https://doi.org/10.1093/bioinformatics/btz468>
- Fukushima and Pollock 2023, CSUBST: <https://doi.org/10.1038/s41559-022-01932-7>
- Lin et al. 2021, *A. laboriosa* genome: <https://doi.org/10.1093/gbe/evab227>
- NCBI current *B. terrestris* assembly: <https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_910591885.1/>
- Emms et al. 2025, OrthoFinder 3 methods: <https://doi.org/10.1101/2025.07.15.664860>
- Emms and Kelly 2019, OrthoFinder phylogenetic orthology inference: <https://doi.org/10.1186/s13059-019-1832-y>
- Buchfink et al. 2021, DIAMOND: <https://doi.org/10.1038/s41592-021-01101-x>
- Deorowicz et al. 2016, FAMSA: <https://doi.org/10.1038/srep33964>
- Minh et al. 2020, IQ-TREE 2: <https://doi.org/10.1093/molbev/msaa015>
- Kalyaanamoorthy et al. 2017, ModelFinder: <https://doi.org/10.1038/nmeth.4285>
- Hoang et al. 2018, UFBoot2: <https://doi.org/10.1093/molbev/msx281>
- Phys.org article that prompted the question: <https://phys.org/news/2026-08-mad-honey-heart-sold-online.html>
