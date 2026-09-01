# *Apis laboriosa* grayanotoxin-tolerance public-data analysis

Updated: 2026-09-01

## Conclusion

This analysis tests proposed genomic and physiological explanations for grayanotoxin tolerance against the *Apis laboriosa* reference genome, five comparative bee proteomes, an independent *A. laboriosa* assembly, a 2023 population-genomic supplement, and public pooled *A. laboriosa* and *A. dorsata* transcriptomes.

The current sequence data provide little support for a species-specific **Para sodium-channel resistance mechanism**. The known grayanotoxin-relevant residues, their surrounding sequence, every catalogued full-length Para isoform, the alternative-splice haplotypes in the only variable transmembrane block, and the abundant whole-worker transcript consensus show no laboriosa-specific protective channel state. The full channel is 99.90% identical to *A. dorsata*. The only two laboriosa-specific residues against the four-other-*Apis* consensus are outside mapped transmembrane helices and far from known grayanotoxin determinants.

If *A. laboriosa* has systemic tolerance, the current results prioritize **toxicokinetic protection** for further testing. Crop staging could limit the amount digested, while midgut transport, efflux, binding, or excretion could limit the free toxin reaching the brain. The only public pooled whole-worker transcriptome contains a high-abundance peritrophin-1-like transcript in *A. laboriosa* (FPKM 7,531.67; rank 12 of 106,176 unigenes), compared with combined peritrophin FPKM 15.29 in *A. dorsata*. Chitin-synthase, ABCC/MRP, and organic-anion-transporter signals are also higher. This result has one pooled sample per species, no toxin treatment, and no biological replication. The peritrophic matrix passes molecules much larger than 370 to 412 Da grayanotoxins, so simple molecular sieving by peritrophin is implausible.

A second lead comes from the population supplement. Its Tibetan “sodium ion transport” enrichment excludes Para: exact cross-assembly mapping places Para on eastern-assembly scaffold 105, while the selected gene models occur on scaffolds 8 and 25. Current RefSeq genes mapping to those two scaffolds suggest DSC1/NaCP60E, a Nach/DEG-ENaC channel, and a COMMD3-like protein. These loci could reflect altitude, sensation, or compensatory excitability. They provide no evidence for a grayanotoxin-resistant Para pocket. The 16 Nav-homologous grayanotoxin positions in DSC1/60E are identical across all six compared bees.

A controlled oral or injected grayanotoxin dose-response study in *A. laboriosa* was not found. Field observation of apparently unaffected foragers and production of toxic honey establish ecological exposure. They do not establish how much toxin reaches the bee's hemolymph or nervous system.

## Hypothesis scorecard

| Hypothesis | Practical test | Result | Current assessment |
| --- | --- | --- | --- |
| Special Para binding-pocket substitution | Map 16 experimentally implicated rat Nav1.4 positions across six bee species | All homologous residues identical among the compared bees | Unsupported at known sites |
| Distant/allosteric Para substitution | Whole-channel comparison against four other *Apis* species | Only N95 and V465 are laboriosa-specific; neither is transmembrane or near a known GTX site | Possible in principle, but little sequence support |
| Protective Para splice isoform | Compare every full-length RefSeq isoform and its variable transmembrane blocks | All three laboriosa DIII S3-S4 haplotypes occur verbatim in susceptible *A. dorsata*, *A. mellifera*, *A. cerana*, and *A. florea* | Catalogued special-isoform explanation rejected |
| Protective RNA editing | Align the complete 6,129-bp Para CDS in the pooled laboriosa unigene consensus | Four substitutions, all synonymous; one A→G; zero amino-acid changes | No fixed, abundant whole-worker protective edit; rare/tissue-specific editing remains untested |
| Para under local selection | Extract 2023 Table S6 and map Para into the same eastern assembly | Selected models are on scaffolds 8/25; all 29 independent Para probes map to scaffold 105 | Rejected for the published sodium-enrichment set |
| Another selected ion channel provides direct resistance | Map likely DSC1/60E and Nach candidates; compare DSC1 homologous sites | DSC1 is 99.37% identical to *A. dorsata* and all 16 mapped positions are invariant; Nach is a distinct DEG/ENaC family | Direct GTX target resistance unsupported; physiological compensation remains possible |
| More detox/efflux/barrier genes | Count parent genes in matched RefSeq annotations | No laboriosa expansion of P450, ABCB, ABCC, ABCG, OATP, MFS, aquaporin, or peritrophin families | Copy-number mechanism unsupported |
| Stronger constitutive gut/clearance program | Map species-specific RefSeq probes into the 2019 pooled transcriptomes and compare submitted FPKM/read counts | Large peritrophin/chitin signal; selective ABCC and OATP elevation; ABCB and ABCG not elevated | Positive descriptive lead; one pooled sample per species |
| Microbiome detoxification | Search available genomic/transcriptomic data | No direct microbiome or GTX-metabolite dataset | Open but unsupported |
| Colony routing/crop containment | Review anatomy and ecological observations | Anatomically plausible; no laboriosa-specific flux measurement | Possible exposure modifier; species-specific effect unmeasured |

## Tests and results

### 1. Para target-site comparison

Published site-directed mutagenesis in rat Nav1.4 implicates 16 positions in the DI and DIV S4-S5 linkers and the S6 helices of all four domains. Protein alignments mapped each position into bee Para.

| Rat Nav1.4 sites | Region | Six-bee result |
| --- | --- | --- |
| K237, L243, S246, K248, K249, L250, S251 | DI S4-S5 | Same homologous residues in every compared bee |
| I433, N434, L437 | DI S6 | Same residues in every compared bee |
| N784 | DII S6 | N in every compared bee |
| S1276 | DIII S6 | S in every compared bee |
| T1463 | DIV S4-S5 | T in every compared bee |
| I1575, F1579, Y1586 | DIV S6 | L/F/Y in every compared bee |

The comparison includes *A. laboriosa*, its close sister *A. dorsata*, susceptible *A. mellifera*, *A. cerana*, *A. florea*, and orally tolerant *Bombus terrestris*. Some bee residues differ from rat, but susceptible bees carry the same states; these shared differences cannot explain species-specific tolerance.

### 2. Para isoforms and distant sequence

Full-length isoforms recovered from the RefSeq proteomes:

| Species | Full-length Para isoforms | Variable known GTX positions |
| --- | ---: | ---: |
| *A. laboriosa* | 17 | 0 |
| *A. dorsata* | 32 | 0 |
| *A. mellifera* | 36 | 0 |
| *A. cerana* | 34 | 0 |
| *A. florea* | 23 | 0 |
| *B. terrestris* | 1 | 0 |

The only variable laboriosa block overlapping a mapped transmembrane region is DIII S3-S4. Laboriosa has three haplotypes there, and all three occur exactly in every other *Apis* reference. No unique catalogued laboriosa isoform was found.

The longest laboriosa Para is XP_043795192.1 (2,052 aa). It is 99.9025% identical to the longest *A. dorsata* model over aligned residues. Against the consensus of *A. dorsata*, *A. mellifera*, *A. cerana*, and *A. florea*, only two laboriosa-specific residues remain: N95 and V465. N95 maps to rat position 85 in a non-transmembrane region; V465 falls in an unmapped insertion/loop. Tolerant *Bombus* carries the other-*Apis* state at both.

### 3. Processed transcriptome test for fixed RNA editing

The public GSE130963 experiment contains one 150-bp paired-end library per species, each made by pooling equal RNA from three whole workers. The analysis uses its de novo unigene assemblies and submitted read-count/FPKM tables.

The top laboriosa Para contig, `AL|c47723_g1`, covers the complete 6,129-bp RefSeq CDS at 99.9347% identity. Its four single-base differences are synonymous. One has the A→G signature compatible with adenosine-to-inosine editing, but it also leaves glutamate unchanged. There are zero nonsynonymous consensus differences.

The top *A. dorsata* Para contig likewise covers its full CDS and has five single-base differences, including one A-to-G nonsynonymous I1681M difference. Editing-like consensus differences occur in both species. Raw read pileups from brain tissue would be needed to distinguish editing from population variants and Trinity assembly errors.

This test rules against a fixed, abundant, whole-worker protein-changing edit. It cannot exclude a rare, neuronal, caste-specific, seasonal, or exposure-induced edit.

### 4. Population-selection signal and scaffold mapping

Cao et al. (2023) resequenced workers from 29 *A. laboriosa* and 28 *A. dorsata* colonies. Table S6 reports three Tibetan-selected gene models in GO:0006814 (“sodium ion transport”; P = 0.01024): two on eastern scaffold 8 and one on scaffold 25.

Using 29 independent exact 51-bp probes sampled from RefSeq Para CDS exons, all 29 map uniquely and concordantly to eastern accession `GWHAOTM00000105` (`OriSeqID=scaffold_105`). None maps to scaffold 8 or 25.

The eastern submission does not publish its MAKER GFF, so the three Table S6 model IDs cannot be assigned directly. RefSeq genes were mapped into scaffolds 8 and 25 with independent exact 31-bp CDS probes and cross-referenced against current RefSeq sodium-related GO annotations. The likely set is:

| Eastern scaffold | Likely current RefSeq gene | Product | Probe support | Interpretation |
| --- | --- | --- | ---: | --- |
| 8 | LOC122714529 | sodium channel protein 60E-like / DSC1 | 16/16 | Nav-like calcium-selective cation-channel family; distinct from canonical Para |
| 8 | LOC122714475 | COMMD3-like | 3/3 | Intracellular protein with a broad, electronically inferred sodium-transport GO term |
| 25 | LOC122718769 | Nach-like | 12/12 | DEG/ENaC/pickpocket-family channel; lacks the four-domain Nav architecture |

The candidate identities are scaffold-and-GO inferences because the MAKER GFF is unavailable. The independent mapping still excludes Para from that enrichment.

The population caveat matters. Chinese laboriosa populations are structured: published FST values between Tibetan and other populations are 0.173-0.203, versus 0.044 between western and eastern Yunnan. The available reference and pooled transcriptome therefore may not represent Nepalese mad-honey populations.

### 5. DSC1/NaCP60E stress test

DSC1 is the only selected-channel candidate with the same broad four-domain architecture as Nav, although functional work identifies this family as a voltage-gated calcium-selective cation channel. The longest laboriosa DSC1/60E model is 99.3723% identical to *A. dorsata* over aligned residues.

Mapping the 16 rat Nav1.4 GTX determinants into DSC1/60E gives exactly the same homologous residues in laboriosa, dorsata, mellifera, cerana, florea, and *Bombus*. Seven laboriosa-specific DSC1 positions remain against the four-other-*Apis* consensus; the nearest is 40 aa from a mapped GTX determinant. There is no direct evidence that grayanotoxin binds DSC1 in the first place.

This leaves a possible indirect role-altered excitability, calcium entry, sensory physiology, or compensation for Para perturbation-but provides no target-site-resistance signal.

### 6. Gene-family copy number and chitin-locus QC

Selected parent-gene counts:

| Family | *A. laboriosa* | *A. dorsata* | *A. mellifera* | *B. terrestris* |
| --- | ---: | ---: | ---: | ---: |
| Cytochrome P450-related | 34 | 36 | 42 | 50 |
| UDP glycosyltransferase | 4 | 9 | 12 | 9 |
| Glutathione S-transferase | 7 | 8 | 7 | 13 |
| ABC/multidrug transporter | 33 | 38 | 35 | 38 |
| ABCB/P-glycoprotein | 3 | 5 | 5 | 5 |
| ABCC/MRP | 9 | 9 | 8 | 9 |
| ABCG | 8 | 9 | 9 | 10 |
| Organic-anion transporter | 6 | 6 | 5 | 6 |
| Major facilitator | 8 | 9 | 9 | 8 |
| Aquaporin | 5 | 5 | 5 | 6 |

No broad laboriosa expansion appears. An initial count of three “chitin synthase” annotations was also misleading: two are credible 1,486-aa and 1,573-aa proteins, while the third is a 114-aa fragment wholly nested inside the span of the 1,573-aa locus. The credible full-length count is therefore two, matching *A. dorsata* and *A. mellifera*.

Copy number cannot test altered promoter activity, tissue localization, transport direction, or substrate specificity.

### 7. Exploratory constitutive-expression screen

Species-specific exact probes mapped annotated genes into each species’ Trinity unigenes. Family totals were deduplicated by unigene. FPKM and counts-per-million (CPM) ratios are descriptive only because there is one untreated pooled whole-worker library per species.

| Category | Laboriosa FPKM | Dorsata FPKM | FPKM ratio | CPM ratio | Reading |
| --- | ---: | ---: | ---: | ---: | --- |
| Peritrophin/chitin-binding | 7,531.67 | 15.29 | 492.6× | 383.6× | Strongest positive lead; one laboriosa transcript is rank 12 genome-wide |
| Chitin synthase annotations | 64.89 | 3.34 | 19.4× | 46.5× | Mainly full-length chs-2; raw category also includes the 114-aa fragment |
| ABCC/MRP | 411.10 | 129.02 | 3.19× | 3.02× | Compatible with selective efflux/clearance |
| Organic-anion transport | 185.31 | 47.75 | 3.88× | 5.09× | Compatible with selective transport, substrate unknown |
| Major facilitator | 188.44 | 117.61 | 1.60× | 1.71× | Modest elevation |
| ABCB/P-glycoprotein | 221.39 | 298.25 | 0.74× | 0.60× | Not elevated |
| ABCG | 101.56 | 110.25 | 0.92× | 0.81× | Not elevated |
| Para | 5.57 | 6.36 | 0.88× | 1.08× | Similar expression |

The P450 total is much higher in laboriosa and is dominated by CYP4G15 (FPKM 5,658.96). CYP4G enzymes are associated with cuticular-hydrocarbon synthesis, so this total is not evidence of GTX metabolism.

The peritrophin result is biologically coherent with a gut-barrier hypothesis: peritrophins determine bee peritrophic-matrix structure and permeability. But it cannot be the whole answer. In *A. mellifera*, 40-kDa dextran crosses the matrix while ≥70-kDa dextran is retained; GTX-I and GTX-III are only about 0.412 and 0.371 kDa. A thicker or compositionally different matrix might alter binding, residence time, microbiota, or epithelial signaling, but simple size exclusion is implausible.

## Working mechanism ranking

1. **Restricted internal exposure plus selective toxicokinetic handling.** Crop storage limits the fraction entering the midgut; epithelial uptake/efflux, binding, and Malpighian clearance can then keep free brain exposure low. The ABCC/OATP and gut-matrix transcript signals support testing this mechanism. Causality is untested.
2. **A specialized gut physiological state.** Peritrophin-1 and chs-2 are candidates for replicated, tissue-specific testing. Altered epithelial physiology or toxin residence is physically more plausible than passive size exclusion of a 0.4-kDa molecule.
3. **Secondary ion-channel compensation.** Population-selected DSC1/60E and Nach-like loci could buffer sensory or excitability changes. The selection could also reflect altitude or foraging ecology. Direct GTX binding has not been shown for either channel.
4. **Specific metabolism or sequestration - open.** No family expansion is required for one enzyme, binding protein, or transporter to matter. Parent GTX persists in honey, so complete detoxification of the crop cargo is unlikely.
5. **Microbiome and colony routing - possible modifiers.** There are no laboriosa-specific data sufficient to rank them higher.
6. **Para structural resistance, catalogued special splicing, or fixed abundant editing - low probability.** Multiple independent sequence tests are negative.

## Revised prediction about physiological cost

A neural-gating trade-off has low support because a resistant Para channel was not identified. If toxicokinetic handling is important, plausible costs include:

- ATP, ion, and water-balance costs of epithelial efflux and Malpighian excretion;
- protein/chitin synthesis and turnover in the midgut barrier;
- reduced immediate access to toxin-rich crop contents as flight fuel;
- transporter competition with pesticides or normal metabolites;
- condition-dependent gut, oxidative, immune, or longevity costs without obvious paralysis.

If DSC1/Nach compensation proves causal, a sensory or neural-performance trade-off could return to the model, but there is currently no evidence that it does.

## Experiments that would resolve the mechanism

The highest-value experiment remains a paired **oral-versus-hemocoel-injection** GTX-I/GTX-III dose-response in *A. laboriosa*, *A. dorsata*, *A. mellifera*, and *B. terrestris*.

| Outcome | Interpretation |
| --- | --- |
| Oral tolerance, injection sensitivity | Crop/gut absorption, binding, efflux, or rapid presystemic clearance is primary |
| Oral and injection tolerance | Systemic metabolism/sequestration, physiological compensation, or an undiscovered target mechanism |
| Similar injection sensitivity plus similar expressed-Para sensitivity | Strong confirmation of the toxicokinetic model |

The same animals should receive isotope-labelled toxin with LC-MS/MS time courses in crop, midgut lumen, gut wall, hemolymph, brain, Malpighian tubules, feces, regurgitated nectar, and honey.

The new gut lead makes the following follow-ups especially valuable:

1. replicated midgut, fat-body, Malpighian-tubule, glial, and brain RNA-seq/proteomics before and after field-realistic GTX exposure;
2. histology and permeability measurements of the laboriosa peritrophic matrix, plus direct GTX flux across isolated midgut;
3. knockdown or inhibition of peritrophin-1, chs-2, the leading ABCCs, and OATPs followed by tissue GTX measurements;
4. long-read brain/antenna RNA sequencing for rare or cell-specific Para edits;
5. heterologous Para electrophysiology, still the definitive check for a cryptic gating effect despite the negative sequence evidence;
6. DSC1/60E and Nach functional assays testing physiological compensation, with no assumption of direct GTX binding.

## Interpretation limits

- No controlled *A. laboriosa* toxin challenge, tissue pharmacokinetic dataset, or expressed-channel assay was found.
- The expression comparison has no biological replicates, uses different collection locations, pools three whole workers per species, and has no toxin treatment. It cannot support a formal differential-expression claim.
- The de novo consensus can miss low-frequency and tissue-specific RNA editing; raw SRA data total roughly 20 GB compressed and 52 Gbp.
- Current Chinese references may not represent Nepalese populations; substantial population structure is documented.
- Annotation releases differ among species, and name-based gene-family counts are screening tests.
- The eastern assembly publishes genomic FASTA without the MAKER GFF, so the three population candidates are inferred from exact scaffold mapping plus current GO annotations.
- Peritrophin abundance is a mechanistic lead. Direct grayanotoxin blocking has not been shown.

## Reproducibility outputs

Key scripts:

- `analyze_nav.py` - Para target residues, whole-channel identities, and species-specific substitutions
- `analyze_para_isoforms.py` and `compare_para_isoforms.py` - isoform variability and exact cross-species haplotypes
- `analyze_para_transcriptomes.py` - complete-CDS transcript-consensus and mismatch analysis
- `extract_population_table_s6.py` - exact extraction of the published sodium-transport candidates
- `map_para_to_eastern.py` - independent Para probes into the eastern assembly
- `map_refseq_genes_to_selected_scaffolds.py` - RefSeq cross-map for eastern scaffolds 8 and 25
- `analyze_selected_ion_channels.py` - DSC1/60E and Nach candidate tests
- `analyze_detox_families.py` and `analyze_barrier_clearance.py` - parent-gene family screens
- `qc_chitin_synthase_loci.py` - overlap and protein-length audit of the apparent chitin-synthase surplus
- `analyze_constitutive_transcript_expression.py` - exact-probe whole-worker expression screen

The `results/` directory contains all residue tables, scaffold mappings, transcript mismatches, gene-family counts, expression mappings/contrasts, and locus QC.

## Principal sources

1. NCBI RefSeq assembly GCF_014066325.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_014066325.1/
2. Lin D. et al. (2021), *Comparative Genomics Reveals Recent Adaptive Evolution in Himalayan Giant Honeybee Apis laboriosa*: https://doi.org/10.1093/gbe/evab227
3. Cao L. et al. (2023), *Population Structure, Demographic History, and Adaptation of Giant Honeybees in China Revealed by Population Genomic Data*: https://doi.org/10.1093/gbe/evad025
4. NCBI GEO/SRA GSE130963 / PRJNA542114: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA542114
5. Tiedeken E. J. et al. (2016), differential bee toxicity of grayanotoxin-containing nectar: https://doi.org/10.1111/1365-2435.12588
6. Maejima H. et al. (2002), DI/DIV S4-S5 determinants: https://pubmed.ncbi.nlm.nih.gov/12150970/
7. Ishii H. et al. (1999), DI S6 determinants: https://pubmed.ncbi.nlm.nih.gov/10603430/
8. Yamaoka K. et al. (2003), DII/DIII/DIV S6 determinants: https://pubmed.ncbi.nlm.nih.gov/12524436/
9. Zhang T. et al. (2013), DSC1 and neuronal excitability: https://pmc.ncbi.nlm.nih.gov/articles/PMC3591268/
10. Rinkevich F. D. et al. (2015), distinct Para and DSC1 roles: https://pmc.ncbi.nlm.nih.gov/articles/PMC4486006/
11. Oliveira A. H. et al. (2019), bee peritrophin and peritrophic-matrix permeability: https://pubmed.ncbi.nlm.nih.gov/31614307/
12. Ahn S. Y. et al. (2022), GTX-I/III in Nepalese mad honey: https://doi.org/10.5806/AST.2022.35.2.82
13. EFSA CONTAM Panel (2023), grayanotoxins in honey: https://doi.org/10.2903/j.efsa.2023.7866
14. PubChem GTX-I and GTX-III records: https://pubchem.ncbi.nlm.nih.gov/compound/9548612 and https://pubchem.ncbi.nlm.nih.gov/compound/11057730
