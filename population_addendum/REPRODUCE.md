# Reproduce the population analysis

This pipeline analyzes all 57 deposited workers. Its input release preserves
296 files in a 337,141,760-byte archive, including all four complete genomes
and the 58,049,708 bytes of selected original FASTQ records.

A fresh clone with a new Python environment restored all 296 files and
reproduced all 29 result files byte for byte. The separately named
`population-inputs-v1.tar` asset is hosted under the existing
[`candidate-inputs-v1` release](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/candidate-inputs-v1);
the GitHub new-tag form was unable to validate the proposed population tag.
The archive contents and hashes are independent of that hosting choice.

## Routine reproduction from the preserved evidence

Use Linux x86-64, Python 3.12 and Java 17. Allow about 4 GB RAM and 3 GB of free
disk space for the preserved inputs, whole-genome references and working files.
The large original source libraries are unnecessary for routine reproduction.

```sh
git clone https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance.git
cd apis-laboriosa-grayanotoxin-tolerance
make -C population_addendum all
```

This creates the Python environment, installs checksum-pinned minimap2 and
BBDuk distributions, restores the analyzed inputs, checks source and dependency
hashes, reconstructs genomic codon positions and masked capture baits, runs the
synthetic controls, verifies all selected reads, maps them to each complete
species genome, marks duplicates within each worker, calls diploid codons,
performs the independent BCFtools check, and produces sample, species and
locality tables. The final stage compares every output with its saved SHA-256.

The archived selected FASTQs preserve the original read identifiers, bases,
qualities and both mates. Genome FASTAs are complete assemblies. The prior
candidate gene models and protein alignment are committed dependency files
checked by `dependencies.sha256`; their reconstruction is documented in
`../candidate_followup/REPRODUCE.md`.

## Recollect a library from the sequencing archive

The source inventory contains all 114 FASTQ URLs, exact byte counts and MD5s.
All source files together are 269,858,255,055 compressed bytes. Each individual
requires roughly 4 to 7 GB of compressed input, plus temporary space.

```sh
make -C population_addendum setup
.venv/bin/python population_addendum/scripts/download_parallel.py --run SRR23343494 --connections 4
.venv/bin/python population_addendum/scripts/collect_reads.py --run SRR23343494 --output-root population_addendum/work/recollection --delete-raw
```

The optional first command downloads each file in at most four disjoint HTTP
ranges, checks every response range, assembles the complete gzip file in order,
and verifies its original MD5. The collector can also download complete files
by itself. It checks the source MD5s again before extraction. Transfer scheduling
does not affect the accepted source bytes.

The separate output directory preserves the pinned evidence. Compare its two
FASTQs and `collection.json` with the archived worker directory. Those files
should match exactly. The BBDuk log contains local paths and run times, so its
bytes can differ between independent collections.

BBDuk finds read names with an exact 21-base bait match. A second pass through
the verified original files copies the complete selected FASTQ records. This
step is required because the tested BBDuk version changed quality scores at N
bases even when quality modification was disabled. The synthetic control tests
exact restoration, retention of unmatched mates and all 64 possible codons.

Downloads reside in the system temporary directory under `bee_population_raw`.
`--delete-raw` removes only the verified temporary source and capture files after
the selected records and collection manifest have been written successfully.

The GitHub workflow `Collect ABCC population reads` performs this stage for all
57 deposited accessions or a supplied comma-separated subset. Its optional
`prior_run` input merges completed libraries from an earlier collection run.
Repeated accessions must have identical selected records and manifests. The
combined artifact must contain all 57 workers when an earlier run is supplied.

## Inspect one worker

```sh
.venv/bin/python population_addendum/scripts/map_population.py --runs SRR23343494 --label one_worker
.venv/bin/python population_addendum/scripts/call_codons.py --label one_worker --output-prefix one_worker_
.venv/bin/python population_addendum/scripts/crosscheck_bcftools.py --label one_worker --prefix one_worker_
.venv/bin/python population_addendum/scripts/summarize_population.py --prefix one_worker_
```

These additional result files are intentionally outside the final checksum
snapshot. Move them into a separate working directory before running the strict
whole-output verification target again.

## Interpret the files

| File | Meaning |
| --- | --- |
| `results/sample_inventory.tsv` | One deposited worker per row, with archive locality, complete source URLs and checksums |
| `results/target_coordinates.tsv` | Codon coordinates in four assemblies and the mapped protein numbering |
| `results/collection_audit.tsv` | Full-source validation and selected-read totals for each worker |
| `results/sample_codon_calls.tsv` | Every sample, site, reference and quality-threshold combination, including failures |
| `results/fragment_evidence.tsv` | The individual fragment observations that support each call |
| `results/genotype_filter_audit.tsv` | Excluded alignments and conflicting-mate counts |
| `results/consensus_genotypes.tsv` | Calls that agree and pass quality checks against both whole-genome references |
| `results/species_site_summary.tsv` | Species counts with explicit callable and missing denominators |
| `results/location_site_summary.tsv` | The same counts separated by collection locality |
| `results/individual_primary_profiles.tsv` | Three-site individual genotypes; these are unphased |
| `results/bcftools_base_crosscheck.tsv` | Independent SNP calls at every tested codon base, including reference sites and missing positions |
| `results/locus_coverage_summary.tsv` | Fragment depth across coding sequence, introns and flanks in both reference mappings |
| `results/locus_coverage_100bp.tsv` | The same coverage described in consecutive 100-base bins |
| `results/collection_batch_merge.json` | Complete cohort membership and the identical overlapping collection between batches |
| `input_sources.tsv` | Exact analyzed files and hashes within the preserved input archive |
| `pilot_expected_reads.sha256` | Read-file hashes from independent local collections of SRR23343443, SRR23343450 and SRR23343494 |
| `results.sha256` | Exact expected generated outputs |

Never replace a failed checksum with a newly calculated value merely to make
the check pass. Preserve the difference, identify its cause and version any
intentional change with an explanation.
