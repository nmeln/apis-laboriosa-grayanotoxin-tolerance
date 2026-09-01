#!/usr/bin/env python3
"""Create one-protein-per-gene FASTA files for an exploratory OrthoFinder run.

The public RefSeq proteomes contain many alternative isoforms.  Counting those as
independent genes would inflate lineage-specific copy numbers, so this script maps
each protein accession to the GFF `gene` attribute and retains the longest protein
per gene.  Ties are resolved lexicographically by accession for determinism.
"""

from __future__ import annotations

import argparse
import csv
import gzip
from pathlib import Path
from urllib.parse import unquote


ADDENDUM_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ADDENDUM_ROOT.parent
DEFAULT_OUT_DIR = ADDENDUM_ROOT / "work/primary_proteomes"

SOURCES = {
    "Apis_laboriosa": (
        PROJECT_ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
        PROJECT_ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
    ),
    "Apis_dorsata": (
        PROJECT_ROOT / "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
        PROJECT_ROOT / "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
    ),
    "Apis_mellifera": (
        PROJECT_ROOT / "genomes/apis_mellifera/GCF_003254395.2_genomic.gff.gz",
        PROJECT_ROOT / "genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz",
    ),
    "Apis_cerana": (
        PROJECT_ROOT / "genomes/apis_cerana/GCF_029169275.2_genomic.gff.gz",
        PROJECT_ROOT / "genomes/apis_cerana/GCF_029169275.2_protein.faa.gz",
    ),
    "Apis_florea": (
        PROJECT_ROOT / "genomes/apis_florea/GCF_048593485.1_genomic.gff.gz",
        PROJECT_ROOT / "genomes/apis_florea/GCF_048593485.1_protein.faa.gz",
    ),
    "Bombus_terrestris": (
        PROJECT_ROOT / "genomes/bombus_terrestris/GCF_000214255.1_genomic.gff.gz",
        PROJECT_ROOT / "genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz",
    ),
}

CURRENT_BOMBUS = (
    ADDENDUM_ROOT / "inputs/bombus_terrestris_current/GCF_910591885.1_iyBomTerr1.2_genomic.gff.gz",
    ADDENDUM_ROOT / "inputs/bombus_terrestris_current/GCF_910591885.1_iyBomTerr1.2_protein.faa.gz",
)


def attrs(text: str) -> dict[str, str]:
    parsed: dict[str, str] = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            parsed[key] = unquote(value)
    return parsed


def protein_to_gene(gff_path: Path) -> dict[str, str]:
    mapping: dict[str, str] = {}
    with gzip.open(gff_path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9 or fields[2] != "CDS":
                continue
            at = attrs(fields[8])
            protein_id = at.get("protein_id")
            gene = at.get("gene")
            if not protein_id or not gene:
                continue
            previous = mapping.setdefault(protein_id, gene)
            if previous != gene:
                raise RuntimeError(f"Conflicting gene mapping for {protein_id}: {previous} vs {gene}")
    return mapping


def read_fasta(path: Path):
    accession: str | None = None
    description = ""
    chunks: list[str] = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if accession is not None:
                    yield accession, description, "".join(chunks).replace("*", "")
                description = line[1:].strip()
                accession = description.split(maxsplit=1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if accession is not None:
        yield accession, description, "".join(chunks).replace("*", "")


def wrap(sequence: str, width: int = 80):
    for start in range(0, len(sequence), width):
        yield sequence[start : start + width]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    parser.add_argument("--current-bombus", action="store_true")
    args = parser.parse_args()

    output_dir = args.output_dir
    output_dir.mkdir(parents=True, exist_ok=True)
    sources = dict(SOURCES)
    if args.current_bombus:
        sources["Bombus_terrestris"] = CURRENT_BOMBUS
    all_rows: list[dict[str, object]] = []
    summary_rows: list[dict[str, object]] = []

    for species, (gff_path, proteome_path) in sources.items():
        mapping = protein_to_gene(gff_path)
        records = list(read_fasta(proteome_path))
        by_gene: dict[str, tuple[str, str, str]] = {}
        unmapped: list[str] = []
        for accession, description, sequence in records:
            gene = mapping.get(accession)
            if gene is None:
                unmapped.append(accession)
                continue
            candidate = (accession, description, sequence)
            current = by_gene.get(gene)
            if current is None or (-len(sequence), accession) < (-len(current[2]), current[0]):
                by_gene[gene] = candidate

        output_path = output_dir / f"{species}.faa"
        with output_path.open("w") as handle:
            for gene in sorted(by_gene):
                accession, description, sequence = by_gene[gene]
                handle.write(f">{accession} gene={gene}\n")
                for line in wrap(sequence):
                    handle.write(line + "\n")
                all_rows.append(
                    {
                        "species": species,
                        "gene": gene,
                        "accession": accession,
                        "length_aa": len(sequence),
                        "description": description,
                    }
                )

        summary_rows.append(
            {
                "species": species,
                "proteins_in_refseq_fasta": len(records),
                "protein_accessions_mapped_to_gene": sum(a in mapping for a, _, _ in records),
                "primary_gene_representatives": len(by_gene),
                "unmapped_protein_accessions": len(unmapped),
                "output_fasta": output_path.name,
            }
        )

    member_name = "primary_proteome_members_current.tsv" if args.current_bombus else "primary_proteome_members.tsv"
    summary_name = "primary_proteome_summary_current.tsv" if args.current_bombus else "primary_proteome_summary.tsv"
    results_dir = ADDENDUM_ROOT / "results"
    results_dir.mkdir(parents=True, exist_ok=True)
    with (results_dir / member_name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(all_rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(all_rows)

    with (results_dir / summary_name).open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summary_rows[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summary_rows)

    for row in summary_rows:
        print("\t".join(str(row[key]) for key in row))


if __name__ == "__main__":
    main()
