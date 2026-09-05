# *Apis laboriosa* grayanotoxin tolerance

Public-data tests of proposed tolerance mechanisms.

Status: preliminary computational analysis. No controlled *A. laboriosa* grayanotoxin challenge data were found. The repository contains hypotheses, tests, generated tables, source accessions, checksums, and deterministic scripts.

## Result

Current public sequence data provide little support for a species-specific grayanotoxin-resistant Para sodium channel in *A. laboriosa*.

**2026-09-05 correction:** The earlier gut-barrier interpretation is superseded by the [transcriptomic addendum](transcriptomic_addendum/). Both source samples are described as “whole body without belly.” Direct analysis of two million R1 reads per species finds black queen cell virus matches in 29.24% of the dorsata subset and 0.0094% of laboriosa. Large chitin-related RNA contrasts survive, but the elevated chitin synthase is closer to a cuticle/trachea reference enzyme. These libraries cannot identify an enhanced midgut barrier.

Reduced internal exposure, transport, metabolism, and excretion remain possible mechanisms. One specific ABCC/MRP group remains an expression candidate, with 1,161 versus 54 raw reads. Species, collection conditions, tissue composition, and physiological state cannot be separated with one pool per species. No GTX transport or causal protection has been demonstrated.

## Tests at a glance

| Hypothesis | Prediction | Test | Result | Assessment |
| --- | --- | --- | --- | --- |
| Altered Para binding site | Unique laboriosa residues at experimentally implicated grayanotoxin positions | Map 16 rat Nav1.4 mutagenesis positions into six bee Para proteins | All 16 homologous residues are identical across the compared bees | No support at known sites |
| Distant Para substitution | Laboriosa-specific changes in or near transmembrane regions | Whole-channel comparison against four other *Apis* species | Para is 99.9025% identical to *A. dorsata*; only N95 and V465 are laboriosa-specific; neither is transmembrane or close to a mapped site | Weak sequence support |
| Protective Para splice form | A laboriosa-specific transmembrane haplotype | Compare all full-length RefSeq isoforms | 17 laboriosa isoforms; zero variable known sites; all three laboriosa DIII S3-S4 haplotypes also occur in four other *Apis* species | No unique catalogued splice form |
| Fixed abundant Para RNA edit | Protein-changing consensus difference in processed laboriosa Para transcript | Align the complete 6,129-base CDS to the pooled transcriptome | Four mismatches, all synonymous | No fixed abundant protein-changing edit detected |
| Para under reported sodium selection | Para maps to selected eastern scaffolds 8 or 25 | Map 29 independent exact Para probes into the eastern assembly | 29 of 29 probes map uniquely to scaffold 105 | Para excluded from that selected set |
| Broad detoxification or barrier expansion | More detoxification, efflux, or barrier genes | Count parent genes in matched annotations | No broad laboriosa expansion | Copy-number explanation unsupported |
| Constitutive gut or clearance program | Higher barrier or transporter transcripts | Original exact-probe screen, followed by source audit and symmetric raw-read mapping | Original family ratios reproduce; large focal RNA differences survive | Gut interpretation withdrawn: belly excluded, viral imbalance, tissue and function unresolved |

## Transcriptomic addendum

The [new report](transcriptomic_addendum/REPORT.md) tests mapping bias, library composition, viral sequence identity, individual transporter candidates, and chitin-synthase class. It also audits an independent GTX-I transcriptome for cross-species reuse. All counts, sensitivity analyses, scripts, input hashes, raw-read subsets, and agent instructions are preserved.

```bash
make -C transcriptomic_addendum all
```

## Comparative addendum

The [`comparative_addendum/`](comparative_addendum/) directory tests whether *A. laboriosa* and the directly measured orally tolerant *Bombus terrestris audax* share a broader coding pattern.

The screen found zero exact shared focal states across 2,050 callable Para residues, no relevant shared gene-family expansion, and no enrichment in prespecified barrier, detoxification, or toxicokinetic gene sets. A weak NHE2/3-labelled exchanger lead remains after wider bee comparisons, with no corrected category-level enrichment and no direct grayanotoxin evidence. Presystemic handling remains open. Oral versus injection exposure with tissue LC-MS/MS would help distinguish it from systemic mechanisms.

The addendum has its own pinned environment, input snapshot, scripts, validation assertions, result manifests, and agent instructions. Run it with:

```bash
make -C comparative_addendum all
```

## Key limitations

- No controlled oral or injected grayanotoxin dose-response was found for *A. laboriosa*.
- The transcriptome comparison has one untreated pool of three individuals per species, described by GEO as “whole body without belly,” with a large viral RNA imbalance.
- FPKM and CPM ratios are descriptive. They are not a replicated differential-expression result.
- The Chinese reference genomes and transcriptomes may not represent Nepalese mad-honey populations.
- The eastern assembly has no published MAKER GFF. Candidate identities on scaffolds 8 and 25 are inferred from exact cross-assembly probes and current RefSeq GO annotations.
- Processed pooled-worker transcripts cannot exclude rare, tissue-specific, seasonal, or exposure-induced RNA editing.
- Peritrophin cannot plausibly exclude a 0.37 to 0.41 kDa molecule by size alone. Honeybee peritrophic matrix permeability has been demonstrated for much larger molecules.

## Most discriminating unresolved test

A paired oral and hemocoel-injection GTX-I/GTX-III dose-response would separate presystemic handling from systemic tolerance.

| Observation | Interpretation |
| --- | --- |
| Oral tolerance with injection sensitivity | Crop handling, gut binding, restricted uptake, epithelial efflux, or rapid presystemic clearance |
| Oral and injection tolerance | Systemic sequestration, metabolism, compensation, or an unidentified target mechanism |

Tissue-resolved LC-MS/MS measurements would show whether toxin reaches the hemolymph and brain.

## Repository map

| Path | Contents |
| --- | --- |
| [`REPORT.md`](REPORT.md) | Full analysis, hypothesis ranking, methods, limitations, and source list |
| [`REPRODUCE.md`](REPRODUCE.md) | Fresh-clone reproduction instructions |
| [`AGENTS.md`](AGENTS.md) | Rules and workflow for coding agents extending or rerunning the analysis |
| [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md) | Input-source terms, citations, and reuse notices |
| [`data_sources.tsv`](data_sources.tsv) | Input paths, accessions, source URLs, sizes, and SHA-256 hashes for the analyzed snapshot |
| [`input_snapshot.sha256`](input_snapshot.sha256) | Hash of the versioned GitHub Release input archive |
| [`results/`](results/) | Generated residue, isoform, scaffold, transcript, copy-number, and expression tables |
| [`scripts/`](scripts/) | Input retrieval, analyses, claim validation, and checksum verification |
| [`results.sha256`](results.sha256) | Expected hashes for every generated result file |
| [`comparative_addendum/`](comparative_addendum/) | Reproducible cross-bee convergence screen, controls, wider taxon checks, and result tables |
| [`transcriptomic_addendum/`](transcriptomic_addendum/) | Tissue correction, raw-read validation, viral burden, focused transport candidates, and functional comparison |

## Key result tables

- [Mapped grayanotoxin positions in bee Para](results/gtx_target_residues.tsv)
- [Whole-channel laboriosa comparisons](results/laboriosa_pairwise_identity.tsv)
- [Laboriosa-specific Para residues](results/laboriosa_unique_genomewide.tsv)
- [Para isoform summary](results/para_isoform_species_summary.tsv)
- [DIII S3-S4 haplotype overlap](results/para_diii_s3_s4_haplotype_overlap.tsv)
- [Processed Para transcript comparison](results/para_transcriptome_alignment_summary.tsv)
- [Para mapping into the eastern assembly](results/para_eastern_scaffold_probe_hits.tsv)
- [Selected sodium-candidate cross-map](results/population_sodium_candidate_crossmap.tsv)
- [DSC1/60E mapped positions](results/selected_60e_gtx_site_residues.tsv)
- [Gene-family copy-number screen](results/detox_family_counts.tsv)
- [Barrier and clearance copy-number screen](results/barrier_clearance_gene_counts.tsv)
- [Chitin-synthase locus QC](results/chitin_synthase_locus_qc.tsv)
- [Original pooled-worker expression contrasts, with interpretation corrected above](results/constitutive_expression_contrasts.tsv)

## Reproduce

Linux and macOS:

```bash
git clone https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance.git
cd apis-laboriosa-grayanotoxin-tolerance
make all
```

`make all` creates a virtual environment, downloads the pinned public inputs, verifies their SHA-256 hashes, reruns every analysis, validates the headline claims against the generated TSV files, and compares all generated outputs with `results.sha256`.

Individual stages:

```bash
make setup          # create .venv and install pinned Python packages
make fetch          # retrieve public input files
make verify-inputs  # compare all 25 inputs with the analyzed snapshot
make results        # rerun the analysis pipeline
make validate       # assert the headline claims from generated tables
make verify-results # compare generated results with the committed snapshot
make check          # compile scripts, validate claims, verify committed results
```

See [`REPRODUCE.md`](REPRODUCE.md) for storage requirements, manual alternatives, and the exact command sequence.

## Data policy

Large third-party genome and transcriptome files are kept out of Git history. A
versioned GitHub Release asset preserves the exact analyzed input snapshot.
`scripts/fetch_inputs.py` tries that pinned archive first and uses NCBI, GWH,
GEO, OUP, and UniProt as a fallback. Every file is checked against
`data_sources.tsv`. See [`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md) for the
source terms. Generated result tables are committed.

## License

Repository code is licensed under the [MIT License](LICENSE). Downloaded and
archived source datasets keep their original terms. See
[`THIRD_PARTY_DATA.md`](THIRD_PARTY_DATA.md).

## Primary sources

1. Lin D. et al. (2021), *Comparative Genomics Reveals Recent Adaptive Evolution in Himalayan Giant Honeybee Apis laboriosa*: https://doi.org/10.1093/gbe/evab227
2. Cao L. et al. (2023), *Population Structure, Demographic History, and Adaptation of Giant Honeybees in China Revealed by Population Genomic Data*: https://doi.org/10.1093/gbe/evad025
3. NCBI RefSeq assembly GCF_014066325.1: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_014066325.1/
4. NCBI GEO/SRA GSE130963 / PRJNA542114: https://www.ncbi.nlm.nih.gov/bioproject/PRJNA542114
5. Tiedeken E. J. et al. (2016), differential grayanotoxin toxicity among bees: https://doi.org/10.1111/1365-2435.12588
6. Maejima H. et al. (2002), DI/DIV S4-S5 determinants: https://pubmed.ncbi.nlm.nih.gov/12150970/
7. Ishii H. et al. (1999), DI S6 determinants: https://pubmed.ncbi.nlm.nih.gov/10603430/
8. Yamaoka K. et al. (2003), DII/DIII/DIV S6 determinants: https://pubmed.ncbi.nlm.nih.gov/12524436/
9. Oliveira A. H. et al. (2019), honeybee peritrophin and matrix permeability: https://pubmed.ncbi.nlm.nih.gov/31614307/

## Scope

This repository reports a computational screen of public data. It contains no live-bee challenge, toxin pharmacokinetics, tissue-specific expression, electrophysiology, or direct transporter assays.
