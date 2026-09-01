#!/usr/bin/env python3
"""Search public A. laboriosa/A. dorsata transcriptome assemblies for Para edits.

The GEO data are pooled-worker de novo consensus assemblies, not read pileups.  They
can therefore rule out a fixed, abundant edited protein form, but cannot establish
low-frequency or tissue-specific RNA editing.  For A. laboriosa, mismatch contexts
are also checked against the independently assembled eastern-Yunnan genome.
"""

from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / ".deps"))

from Bio import SeqIO  # type: ignore  # noqa: E402
from Bio.Align import PairwiseAligner  # type: ignore  # noqa: E402
from Bio.Seq import Seq  # type: ignore  # noqa: E402


CONFIG = {
    "Apis_laboriosa": {
        "protein": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
        "rna": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz",
        "gff": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
        "unigene": ROOT / "references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz",
        "fpkm": ROOT / "references/transcriptome_2019/GSM3757258_AL.Readcount_FPKM.txt.gz",
    },
    "Apis_dorsata": {
        "protein": ROOT / "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
        "rna": ROOT / "genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz",
        "gff": ROOT / "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
        "unigene": ROOT / "references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz",
        "fpkm": ROOT / "references/transcriptome_2019/GSM3757259_AD.Readcount_FPKM.txt.gz",
    },
}
PARA_RE = re.compile(r"product=sodium channel protein (?:para|paralytic)", re.I)


def records(path: Path):
    with gzip.open(path, "rt") as handle:
        return list(SeqIO.parse(handle, "fasta"))


def revcomp(seq: str):
    return str(Seq(seq).reverse_complement())


def para_pairs(gff: Path):
    pairs = {}
    with gzip.open(gff, "rt") as handle:
        for line in handle:
            if line.startswith("#") or "\tCDS\t" not in line or not PARA_RE.search(line):
                continue
            attributes = line.rstrip().split("\t", 8)[8]
            parent = re.search(r"(?:^|;)Parent=rna-([^;]+)", attributes)
            protein = re.search(r"(?:^|;)protein_id=([^;]+)", attributes)
            if parent and protein:
                pairs[parent.group(1)] = protein.group(1)
    return pairs


def extract_cds(mrna: str, protein: str):
    for frame in range(3):
        translated = str(Seq(mrna[frame:]).translate())
        aa_start = translated.find(protein)
        if aa_start >= 0:
            nt_start = frame + aa_start * 3
            return mrna[nt_start : nt_start + len(protein) * 3]
    raise RuntimeError("Protein not found exactly in its annotated mRNA")


def load_expression(path: Path):
    result = {}
    with gzip.open(path, "rt") as handle:
        next(handle)
        for line in handle:
            gene_id, _sample, read_count, fpkm = line.rstrip().split("\t")
            result[gene_id] = (float(read_count), float(fpkm))
    return result


def local_aligner():
    a = PairwiseAligner()
    a.mode = "local"
    a.match_score = 2.0
    a.mismatch_score = -2.0
    a.open_gap_score = -10.0
    a.extend_gap_score = -1.0
    return a


def find_candidates(unigenes, canonical_cds, k=31, stride=19):
    probes = [(i, canonical_cds[i : i + k]) for i in range(0, len(canonical_cds) - k + 1, stride)]
    direct = {probe: i for i, probe in probes}
    reverse = {revcomp(probe): i for i, probe in probes}
    hits = []
    for record in unigenes:
        seq = str(record.seq).upper()
        d = set()
        r = set()
        for i in range(len(seq) - k + 1):
            word = seq[i : i + k]
            if word in direct:
                d.add(direct[word])
            if word in reverse:
                r.add(reverse[word])
        count = max(len(d), len(r))
        if count >= 5:
            hits.append((record, count, "+" if len(d) >= len(r) else "-"))
    return sorted(hits, key=lambda x: x[1], reverse=True)


def best_isoform(oriented_seq, isoforms, k=31, stride=13):
    scores = []
    for transcript_id, protein_id, cds, protein in isoforms:
        probes = {cds[i : i + k] for i in range(0, len(cds) - k + 1, stride)}
        scores.append((sum(probe in oriented_seq for probe in probes), transcript_id, protein_id, cds, protein))
    return max(scores, key=lambda x: x[0])


def load_eastern_para_scaffold():
    path = ROOT / "genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz"
    for record in records(path):
        if "OriSeqID=scaffold_105" in record.description:
            return str(record.seq).upper()
    raise RuntimeError("Eastern-Yunnan scaffold_105 not found")


def main():
    out = ROOT / "results"
    out.mkdir(exist_ok=True)
    eastern = load_eastern_para_scaffold()
    summary_rows = []
    mismatch_rows = []
    candidate_rows = []
    aligner = local_aligner()

    for species, cfg in CONFIG.items():
        proteins = {r.id: str(r.seq).replace("*", "") for r in records(cfg["protein"])}
        mrnas = {r.id: str(r.seq).upper() for r in records(cfg["rna"])}
        pairs = para_pairs(cfg["gff"])
        isoforms = []
        for transcript_id, protein_id in pairs.items():
            if transcript_id in mrnas and protein_id in proteins:
                cds = extract_cds(mrnas[transcript_id], proteins[protein_id])
                isoforms.append((transcript_id, protein_id, cds, proteins[protein_id]))
        if not isoforms:
            raise RuntimeError(f"No Para transcript/protein pairs for {species}")
        canonical = max(isoforms, key=lambda x: len(x[2]))
        expression = load_expression(cfg["fpkm"])
        unigenes = records(cfg["unigene"])
        candidates = find_candidates(unigenes, canonical[2])
        if not candidates:
            raise RuntimeError(f"No Para unigene candidate for {species}")

        for rank, (record, probe_hits, orientation) in enumerate(candidates[:10], start=1):
            read_count, fpkm = expression.get(record.id, (0.0, 0.0))
            candidate_rows.append(
                {
                    "species": species,
                    "rank": rank,
                    "unigene": record.id,
                    "length_nt": len(record.seq),
                    "probe_hits": probe_hits,
                    "orientation": orientation,
                    "read_count": read_count,
                    "fpkm": fpkm,
                }
            )

        record, probe_hits, orientation = candidates[0]
        query = str(record.seq).upper()
        if orientation == "-":
            query = revcomp(query)
        exact_probe_score, transcript_id, protein_id, ref_cds, ref_protein = best_isoform(query, isoforms)
        alignment = aligner.align(ref_cds, query)[0]
        ref_blocks, query_blocks = alignment.aligned
        aligned_bases = 0
        matches = 0
        mismatches = []
        for (r0, r1), (q0, q1) in zip(ref_blocks, query_blocks):
            rs = ref_cds[r0:r1]
            qs = query[q0:q1]
            aligned_bases += len(rs)
            for offset, (rb, qb) in enumerate(zip(rs, qs)):
                if rb == qb:
                    matches += 1
                else:
                    mismatches.append((int(r0 + offset), int(q0 + offset), rb, qb))

        for ref_i, query_i, ref_base, query_base in mismatches:
            codon_start = (ref_i // 3) * 3
            ref_codon = ref_cds[codon_start : codon_start + 3]
            alt_codon_chars = list(ref_codon)
            alt_codon_chars[ref_i % 3] = query_base
            alt_codon = "".join(alt_codon_chars)
            ref_aa = str(Seq(ref_codon).translate()) if len(ref_codon) == 3 else "-"
            alt_aa = str(Seq(alt_codon).translate()) if len(alt_codon) == 3 else "-"
            radius = 15
            ref_context = ref_cds[max(0, ref_i - radius) : ref_i + radius + 1]
            query_context = query[max(0, query_i - radius) : query_i + radius + 1]
            context_complete = len(ref_context) == 31 and len(query_context) == 31
            if species == "Apis_laboriosa" and context_complete:
                ref_in_east = ref_context in eastern or revcomp(ref_context) in eastern
                query_in_east = query_context in eastern or revcomp(query_context) in eastern
            else:
                ref_in_east = "not_tested"
                query_in_east = "not_tested"
            mismatch_rows.append(
                {
                    "species": species,
                    "unigene": record.id,
                    "reference_transcript": transcript_id,
                    "reference_protein": protein_id,
                    "cds_position_1based": ref_i + 1,
                    "protein_position_1based": ref_i // 3 + 1,
                    "ref_base": ref_base,
                    "transcriptome_base": query_base,
                    "a_to_g_signature": ref_base == "A" and query_base == "G",
                    "ref_codon": ref_codon,
                    "transcriptome_codon": alt_codon,
                    "ref_aa": ref_aa,
                    "transcriptome_aa": alt_aa,
                    "nonsynonymous": ref_aa != alt_aa,
                    "ref_context_in_eastern_genome": ref_in_east,
                    "transcriptome_context_in_eastern_genome": query_in_east,
                }
            )

        summary_rows.append(
            {
                "species": species,
                "para_refseq_isoforms": len(isoforms),
                "top_unigene": record.id,
                "top_unigene_length_nt": len(query),
                "read_count": expression.get(record.id, (0.0, 0.0))[0],
                "fpkm": expression.get(record.id, (0.0, 0.0))[1],
                "best_reference_transcript": transcript_id,
                "best_reference_protein": protein_id,
                "exact_probe_score": exact_probe_score,
                "reference_cds_length_nt": len(ref_cds),
                "aligned_bases": aligned_bases,
                "reference_coverage": f"{aligned_bases / len(ref_cds):.6f}",
                "identity_over_aligned": f"{matches / aligned_bases:.6f}",
                "single_base_mismatches": len(mismatches),
                "a_to_g_mismatches": sum(rb == "A" and qb == "G" for _, _, rb, qb in mismatches),
                "nonsynonymous_mismatches": sum(
                    row["species"] == species and row["unigene"] == record.id and row["nonsynonymous"]
                    for row in mismatch_rows
                ),
            }
        )

    for filename, rows in [
        ("para_transcriptome_candidates.tsv", candidate_rows),
        ("para_transcriptome_alignment_summary.tsv", summary_rows),
    ]:
        with (out / filename).open("w", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)

    mismatch_fields = [
        "species", "unigene", "reference_transcript", "reference_protein",
        "cds_position_1based", "protein_position_1based", "ref_base", "transcriptome_base",
        "a_to_g_signature", "ref_codon", "transcriptome_codon", "ref_aa", "transcriptome_aa",
        "nonsynonymous", "ref_context_in_eastern_genome", "transcriptome_context_in_eastern_genome",
    ]
    with (out / "para_transcriptome_mismatches.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=mismatch_fields, delimiter="\t")
        writer.writeheader()
        writer.writerows(mismatch_rows)

    for row in summary_rows:
        print(row)


if __name__ == "__main__":
    main()
