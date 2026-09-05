# Population follow-up instructions

Read the root AGENTS.md and the candidate_followup instructions and report,
then this directory's ANALYSIS_PLAN.md, README.md, REPORT.md, REPRODUCE.md and
relevant generating scripts before extending the work.

- The primary sites are laboriosa protein positions 254, 549 and 1134.
  Position numbering must be mapped through the saved protein alignment and
  genomic coding exons, including reverse-strand coordinates and split codons.
- Use the 57 deposited accessions and reconcile sample aliases. A row in a
  supplementary table is not evidence that its raw reads were deposited.
- The paper describes one worker per colony. Use diploid individual calls.
  Do not count paired mates, PCR duplicates or reference genomes as extra bees.
- Mask tested codons in capture baits and retain both mates. Keep whole-genome
  competitors in downstream mapping. A capture match is not a genotype.
- Check complete source-file byte counts and MD5s before accepting extraction.
  Failed or partial downloads do not supply a completed sample.
- Preserve exact selected read sequences, qualities, names and source checksums.
- Missing or ambiguous calls must remain missing. Report callable denominators
  at every site and quality threshold. Do not infer reference alleles from
  absent VCF rows or zero coverage.
- Keep reference swaps, duplicate/overlap handling, strand evidence, unexpected
  codons and prespecified quality sensitivity results.
- Do not label observed cohort consistency as species-wide fixation or toxin
  adaptation. These Chinese samples do not establish the Nepalese population
  state. No GTX phenotype or direct transport assay is available.
- Never weaken a validation assertion to obtain a passing run. Inspect and
  explain differences before deliberately versioning new result checksums.
- Preserve full analyzed evidence in a versioned release. Downloaded FASTQs,
  indexes, BAMs, intermediate files and input archives stay out of Git history.

New analyses and results must be reproducible from the archived evidence.
Also retain a separately runnable full-FASTQ collection stage so the extraction
can be independently checked against the original sequencing archives.

The source collection used two batches: the first eight workers were retained
from run 33993933619, and run 33994894621 collected SRR23343446 through
SRR23343494 with faster checksum-verified range downloads. The first queue was
cancelled after its required workers completed. The second workflow combines
completed artifacts and requires all 57 accessions. Its repeated SRR23343446
collection must have identical selected records. The independent local
SRR23343443 extraction is pinned in `pilot_expected_reads.sha256`.

A low-coverage worker, SRR23343455 (HD17), prompted the descriptive exon,
intron and flank coverage check. Do not turn a depth deficit at this captured
locus into a gene-deletion claim. Keep the original quality thresholds and
retain missing calls. The complete-source files were checked before capture.

The completed cohort has 28 callable laboriosa workers at each primary site;
dorsata has 27, 26 and 28 at positions 254, 549 and 1134. Keep the explicit
exceptions in REPORT.md, including D16's single TTG fragment among 18 and
D31's single synonymous TTC fragment among 16. The original quality rules
exclude these calls. Do not promote a minority fragment into a diploid allele.

`pilot_expected_reads.sha256` preserves the initial SRR23343443 collection and
two additional independent local collections, SRR23343450 and SRR23343494.
All six selected FASTQs match the remote evidence. The batch merge record
also preserves the independent repeated SRR23343446 collection.

For a complete rerun, use `make -C population_addendum all` from the repository
root. It restores the 296-file input archive, verifies dependencies, produces
the full result set and checks the observed findings and every result hash.
The separate collection option `--output-root` lets an independent check write
under `work/` without replacing archived evidence or its timing-dependent log.
