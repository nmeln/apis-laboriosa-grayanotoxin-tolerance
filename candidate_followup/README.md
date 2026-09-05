# ABCC candidate follow-up

2026-09-05. Follow-up to the [transcriptomic audit](../transcriptomic_addendum/REPORT.md).

The specific transporter candidate survives independent genome and RNA-mapping
checks. A separate study detected its *Apis mellifera* counterpart in midgut
membrane preparations. This supports studying the gene's transport activity.
Grayanotoxin transport and protection remain unmeasured.

| Question | Result |
| --- | --- |
| Is there another dorsata genome? | Yes. Thai assembly `GCA_009792835.1` is independent of the Malaysian reference `GCF_000469605.1`. Both contain the complete 14-exon candidate. |
| Is the gene missing or broken in dorsata? | Both recovered dorsata proteins are 1,333 amino acids long and differ at one residue. Neither model has an internal stop or ambiguous coding base. |
| Does the RNA contrast depend on the reference? | Host-normalized laboriosa/dorsata ratios are 17.43 using the laboriosa reference and 17.60 using dorsata. |
| Can an aligner-free control reproduce the direction? | Identical, gene-specific 100-nt markers count 258 laboriosa versus 15 dorsata reads, a raw ratio of 17.2. These are a subset of the same libraries. |
| Is there independent gut evidence? | The exact mellifera counterpart, `XP_393750.5`, is reported in all three adult and all three larval midgut membrane preparations. |
| Are there specific coding candidates? | L254, L549 and F1134 differ from the other five bee references and agree across two laboriosa assemblies. The full assembled laboriosa transcript encodes the reference protein exactly. |
| Is there a new focal-species gut RNA dataset? | None found in the saved ENA/SRA searches. Each species has one indexed transcriptomic run, the original belly-excluded pool. Related-bee gut and whole-individual datasets exist. |

The candidate is `OG0000499`: laboriosa `LOC122718161 / XP_043798878.1`,
dorsata `LOC102671280 / XP_006621542.1`, and mellifera
`LOC410269 / XP_393750.5`. The RefSeq product label is probable multidrug
resistance-associated protein lethal(2)03659. Functional conclusions should
use the measured evidence rather than the product name.

The RNA comparison still has one untreated pool per species, different
collection conditions, excluded bellies, and a large viral RNA imbalance.
The gut-protein result concerns mellifera. It does not locate the original
laboriosa RNA signal in the gut.

## Run

```bash
make -C candidate_followup all
```

Linux x86-64, Python 3.11+, two parallel jobs with two threads each, and about
5 GB free disk space are sufficient. The workflow restores checked inputs,
rebuilds the pooled host reference, runs all tests, validates headline values,
and compares every result byte with its committed checksum.

- [Full report](REPORT.md)
- [Reproduction instructions](REPRODUCE.md)
- [Analysis plan and later additions](ANALYSIS_PLAN.md)
- [Agent instructions](AGENTS.md)
- [Generated evidence tables and sequences](results/)
- [Input sources and SHA-256 hashes](input_sources.tsv)
- [Archived inputs](https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/tag/candidate-inputs-v1)
- [Third-party data terms](THIRD_PARTY_DATA.md)

The variant FASTAs are **untested sequence designs**. They have not been
synthesized, expressed or assayed.
