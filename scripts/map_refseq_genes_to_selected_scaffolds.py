#!/usr/bin/env python3
"""Map RefSeq A. laboriosa gene models to population-study scaffolds 8 and 25.

The population paper's Table S6 lists three GO:0006814 candidates only by
MAKER IDs on eastern-assembly scaffolds 8 and 25.  No public GFF accompanies
that assembly.  This script uses independent exact 31-bp probes from RefSeq
CDS exons to identify RefSeq gene models on those two scaffolds and makes the
candidate neighborhood interpretable without assuming scaffold-number
equivalence between assemblies.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
REF_GENOME = ROOT / "genomes/apis_laboriosa/GCF_014066325.1_genomic.fna.gz"
REF_GFF = ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz"
EASTERN_GENOME = ROOT / "genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz"
REF_GAF = ROOT / "references/GCF_014066325.1_gene_ontology.gaf.gz"
OUT = ROOT / "results/refseq_genes_on_population_sodium_scaffolds.tsv"
CANDIDATE_OUT = ROOT / "results/population_sodium_candidate_crossmap.tsv"
TARGET_ORIGINAL_IDS = {"scaffold_8", "scaffold_25"}
K = 31
MAX_PROBES_PER_GENE = 16
SODIUM_RELATED_GO_TERMS = {
    "GO:0006814",  # sodium ion transport
    "GO:0035725",  # sodium ion transmembrane transport
    "GO:0005248",  # voltage-gated sodium channel activity
    "GO:0015280",  # ligand-gated sodium channel activity
}


def parse_attrs(text: str) -> dict[str, str]:
    result = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def fasta_records(path: Path):
    with gzip.open(path, "rt") as handle:
        header = None
        chunks = []
        for line in handle:
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks).upper()
                header = line[1:].strip()
                chunks = []
            else:
                chunks.append(line.strip())
        if header is not None:
            yield header, "".join(chunks).upper()


def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def main() -> None:
    genes = {}
    products = defaultdict(set)
    cds_intervals = defaultdict(list)
    with gzip.open(REF_GFF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attr_text = fields
            attrs = parse_attrs(attr_text)
            if feature == "gene" and attrs.get("ID"):
                genes[attrs["ID"]] = {
                    "refseq_gene_parent": attrs["ID"],
                    "refseq_gene_name": attrs.get("gene", attrs.get("Name", attrs["ID"])),
                    "refseq_scaffold": seqid,
                    "refseq_start": int(start),
                    "refseq_end": int(end),
                    "strand": strand,
                }
            elif feature in {"mRNA", "transcript"} and attrs.get("Parent") and attrs.get("product"):
                products[attrs["Parent"]].add(attrs["product"])
            elif feature == "CDS" and attrs.get("gene"):
                cds_intervals[f"gene-{attrs['gene']}"] .append((seqid, int(start), int(end)))

    needed_scaffolds = {interval[0] for rows in cds_intervals.values() for interval in rows}
    ref_sequences = {}
    for header, sequence in fasta_records(REF_GENOME):
        accession = header.split()[0]
        if accession in needed_scaffolds:
            ref_sequences[accession] = sequence

    gene_probes = defaultdict(set)
    for gene_id, intervals in cds_intervals.items():
        # One central probe per distinct CDS interval, sampled evenly when there
        # are many isoforms/exons.
        unique_intervals = sorted(set(intervals))
        candidates = []
        for seqid, start, end in unique_intervals:
            segment = ref_sequences[seqid][start - 1 : end]
            if len(segment) >= K and "N" not in segment:
                offset = (len(segment) - K) // 2
                candidates.append(segment[offset : offset + K])
        if len(candidates) > MAX_PROBES_PER_GENE:
            indexes = [round(i * (len(candidates) - 1) / (MAX_PROBES_PER_GENE - 1)) for i in range(MAX_PROBES_PER_GENE)]
            candidates = [candidates[i] for i in indexes]
        gene_probes[gene_id].update(candidates)

    probe_to_genes = defaultdict(set)
    for gene_id, probes in gene_probes.items():
        for probe in probes:
            probe_to_genes[probe].add(gene_id)
            probe_to_genes[revcomp(probe)].add(gene_id)

    hits = defaultdict(lambda: defaultdict(set))
    hit_positions = defaultdict(lambda: defaultdict(list))
    target_accessions = {}
    for header, sequence in fasta_records(EASTERN_GENOME):
        fields = header.split()
        metadata = dict(part.split("=", 1) for part in fields[1:] if "=" in part)
        original_id = metadata.get("OriSeqID")
        if original_id not in TARGET_ORIGINAL_IDS:
            continue
        accession = fields[0]
        target_accessions[original_id] = accession
        for position in range(len(sequence) - K + 1):
            kmer = sequence[position : position + K]
            for gene_id in probe_to_genes.get(kmer, ()):
                # The probe itself, not every genomic occurrence, is the
                # independent evidence unit.
                hits[original_id][gene_id].add(kmer if kmer in gene_probes[gene_id] else revcomp(kmer))
                hit_positions[original_id][gene_id].append(position + 1)

    go_terms = defaultdict(set)
    with gzip.open(REF_GAF, "rt") as handle:
        for line in handle:
            if line.startswith("!"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) >= 9:
                go_terms[f"gene-LOC{fields[1]}"] .add(fields[4])

    rows = []
    for original_id in sorted(TARGET_ORIGINAL_IDS, key=lambda value: int(value.split("_")[-1])):
        for gene_id, matched in hits[original_id].items():
            total = len(gene_probes[gene_id])
            if len(matched) < 2:
                continue
            row = {
                "eastern_original_scaffold": original_id,
                "eastern_accession": target_accessions[original_id],
                "eastern_probe_min_position": min(hit_positions[original_id][gene_id]),
                "eastern_probe_max_position": max(hit_positions[original_id][gene_id]),
                **genes[gene_id],
                "product": min(
                    products.get(gene_id, {"unlabelled"}),
                    key=lambda label: (len(label), label.casefold(), label),
                ),
                "exact_probe_hits": len(matched),
                "total_probes": total,
                "probe_fraction": f"{len(matched) / total:.6f}" if total else "NA",
                "refseq_sodium_related_go_terms": ",".join(
                    sorted(go_terms.get(gene_id, set()) & SODIUM_RELATED_GO_TERMS)
                ),
                "refseq_sodium_go_match": bool(
                    go_terms.get(gene_id, set()) & SODIUM_RELATED_GO_TERMS
                ),
                "sodium_or_ion_keyword": bool(
                    re.search(
                        r"\bsodium\b|\bcation\b|\bion channel\b|\bion transport|\batpase\b",
                        " ".join(products.get(gene_id, set())),
                        re.I,
                    )
                ),
            }
            rows.append(row)

    rows.sort(
        key=lambda row: (
            int(str(row["eastern_original_scaffold"]).split("_")[-1]),
            not bool(row["refseq_sodium_go_match"]),
            not bool(row["sodium_or_ion_keyword"]),
            -int(row["exact_probe_hits"]),
            str(row["product"]),
        )
    )
    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    candidate_rows = []
    for row in rows:
        if not row["refseq_sodium_go_match"]:
            continue
        candidate_rows.append(
            {
                "eastern_original_scaffold": row["eastern_original_scaffold"],
                "eastern_accession": row["eastern_accession"],
                "eastern_probe_min_position": row["eastern_probe_min_position"],
                "eastern_probe_max_position": row["eastern_probe_max_position"],
                "likely_refseq_gene": row["refseq_gene_name"],
                "likely_product": row["product"],
                "refseq_sodium_related_go_terms": row["refseq_sodium_related_go_terms"],
                "exact_probe_hits": row["exact_probe_hits"],
                "total_probes": row["total_probes"],
                "assignment_status": "scaffold-and-GO candidate; eastern MAKER GFF unavailable, so not a direct maker-ID assignment",
            }
        )
    with CANDIDATE_OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=candidate_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(candidate_rows)

    for original_id in sorted(TARGET_ORIGINAL_IDS):
        scaffold_rows = [row for row in rows if row["eastern_original_scaffold"] == original_id]
        go_rows = [row for row in scaffold_rows if row["refseq_sodium_go_match"]]
        print(f"{original_id}: {len(scaffold_rows)} RefSeq gene mappings; {len(go_rows)} sodium-GO candidates")
        for row in go_rows:
            print(
                f"  {row['refseq_gene_name']} ({row['exact_probe_hits']}/{row['total_probes']} probes; "
                f"{row['refseq_sodium_related_go_terms']}): {row['product']}"
            )


if __name__ == "__main__":
    main()
