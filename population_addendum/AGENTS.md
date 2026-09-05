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
