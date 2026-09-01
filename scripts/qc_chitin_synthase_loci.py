#!/usr/bin/env python3
"""Audit apparent chitin-synthase copy number in the A. laboriosa GFF.

The annotation-level family screen counts parent gene IDs.  This script checks
whether any of those gene models occupy the same genomic locus, which would make
the raw count unsafe to interpret as a lineage-specific duplication.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote

from Bio import SeqIO


ROOT = Path(__file__).resolve().parents[1]
GFF = ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz"
PROTEINS = ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz"
OUT = ROOT / "results/chitin_synthase_locus_qc.tsv"
PATTERN = re.compile(r"chitin synthase", re.I)


def parse_attrs(text: str) -> dict[str, str]:
    attrs = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            attrs[key] = unquote(value)
    return attrs


def main() -> None:
    genes: dict[str, dict[str, object]] = {}
    products: dict[str, set[str]] = defaultdict(set)
    protein_ids: dict[str, set[str]] = defaultdict(set)

    with gzip.open(GFF, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9:
                continue
            seqid, _source, feature, start, end, _score, strand, _phase, attr_text = fields
            attrs = parse_attrs(attr_text)
            if feature == "gene":
                gene_id = attrs.get("ID")
                if gene_id:
                    genes[gene_id] = {
                        "gene_parent": gene_id,
                        "gene_name": attrs.get("gene", attrs.get("Name", gene_id)),
                        "scaffold": seqid,
                        "start": int(start),
                        "end": int(end),
                        "strand": strand,
                    }
            elif feature in {"mRNA", "transcript"}:
                parent = attrs.get("Parent")
                product = attrs.get("product")
                if parent and product and PATTERN.search(product):
                    products[parent].add(product)
            elif feature == "CDS" and attrs.get("gene") and attrs.get("protein_id"):
                protein_ids[f"gene-{attrs['gene']}"] .add(attrs["protein_id"])

    protein_lengths = {}
    with gzip.open(PROTEINS, "rt") as handle:
        for record in SeqIO.parse(handle, "fasta"):
            protein_lengths[record.id] = len(record.seq)

    hits = []
    for parent, labels in products.items():
        row = dict(genes[parent])
        row["product"] = min(
            labels,
            key=lambda label: (len(label), label.casefold(), label),
        )
        accessions = sorted(protein_ids.get(parent, set()))
        lengths = [protein_lengths[accession] for accession in accessions if accession in protein_lengths]
        row["protein_accessions"] = ",".join(accessions)
        row["protein_lengths_aa"] = ",".join(str(length) for length in lengths)
        row["max_protein_length_aa"] = max(lengths) if lengths else ""
        row["full_length_chitin_synthase_like"] = bool(lengths and max(lengths) >= 1000)
        hits.append(row)
    hits.sort(key=lambda row: (row["scaffold"], row["start"], row["end"]))

    # Build connected components of overlapping gene spans on the same scaffold.
    adjacency: dict[str, set[str]] = {str(row["gene_parent"]): set() for row in hits}
    pair_details: dict[tuple[str, str], tuple[int, float, str]] = {}
    for i, left in enumerate(hits):
        for right in hits[i + 1 :]:
            if left["scaffold"] != right["scaffold"]:
                continue
            overlap = max(0, min(int(left["end"]), int(right["end"])) - max(int(left["start"]), int(right["start"])) + 1)
            if not overlap:
                continue
            left_id, right_id = str(left["gene_parent"]), str(right["gene_parent"])
            adjacency[left_id].add(right_id)
            adjacency[right_id].add(left_id)
            shorter = min(int(left["end"]) - int(left["start"]) + 1, int(right["end"]) - int(right["start"]) + 1)
            if int(left["start"]) <= int(right["start"]) and int(left["end"]) >= int(right["end"]):
                relation = f"{right_id} nested in {left_id}"
            elif int(right["start"]) <= int(left["start"]) and int(right["end"]) >= int(left["end"]):
                relation = f"{left_id} nested in {right_id}"
            else:
                relation = "partial overlap"
            pair_details[(left_id, right_id)] = (overlap, overlap / shorter, relation)

    component_id: dict[str, int] = {}
    component = 0
    for node in adjacency:
        if node in component_id:
            continue
        component += 1
        stack = [node]
        component_id[node] = component
        while stack:
            current = stack.pop()
            for neighbor in adjacency[current]:
                if neighbor not in component_id:
                    component_id[neighbor] = component
                    stack.append(neighbor)

    output_rows = []
    for row in hits:
        gene_id = str(row["gene_parent"])
        overlaps = []
        relations = []
        fractions = []
        for (left_id, right_id), (_bases, fraction, relation) in pair_details.items():
            if gene_id in {left_id, right_id}:
                overlaps.append(right_id if gene_id == left_id else left_id)
                relations.append(relation)
                fractions.append(f"{fraction:.6f}")
        output_rows.append(
            {
                **row,
                "overlap_component": component_id[gene_id],
                "overlaps_gene": ",".join(overlaps),
                "overlap_fraction_of_shorter": ",".join(fractions),
                "overlap_relation": "; ".join(relations),
            }
        )

    OUT.parent.mkdir(exist_ok=True)
    with OUT.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=output_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(output_rows)

    print(f"annotated_chitin_synthase_genes={len(hits)}")
    print(f"nonoverlapping_locus_components={len(set(component_id.values()))}")
    for row in output_rows:
        if row["overlap_relation"]:
            print(row["overlap_relation"])


if __name__ == "__main__":
    main()
