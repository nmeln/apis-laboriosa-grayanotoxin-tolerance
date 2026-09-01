#!/usr/bin/env python3
"""Audit local inputs without changing the canonical snapshot manifest."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
DEFAULT_OUT = ROOT / "local_input_manifest.tsv"

SOURCES = [
    ("genomes/apis_laboriosa/GCF_014066325.1_genomic.fna.gz", "GCF_014066325.1", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/GCF_014066325.1_ASM1406632v1/GCF_014066325.1_ASM1406632v1_genomic.fna.gz"),
    ("genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz", "GCF_014066325.1", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/GCF_014066325.1_ASM1406632v1/GCF_014066325.1_ASM1406632v1_genomic.gff.gz"),
    ("genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz", "GCF_014066325.1", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/GCF_014066325.1_ASM1406632v1/GCF_014066325.1_ASM1406632v1_protein.faa.gz"),
    ("genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz", "GCF_014066325.1", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/GCF_014066325.1_ASM1406632v1/GCF_014066325.1_ASM1406632v1_rna.fna.gz"),
    ("references/GCF_014066325.1_gene_ontology.gaf.gz", "GCF_014066325.1 GO annotation", "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/014/066/325/GCF_014066325.1_ASM1406632v1/GCF_014066325.1_ASM1406632v1_gene_ontology.gaf.gz"),
    ("genomes/apis_laboriosa/eastern_yunnan/GWHAOTM00000000.genome.fasta.gz", "GWHAOTM00000000", "https://download.cncb.ac.cn/gwh/Animals/Apis_laboriosa_A.laboriosa_scaffold_GWHAOTM00000000/GWHAOTM00000000.genome.fasta.gz"),
    ("genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz", "GCF_000469605.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000469605.1/"),
    ("genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz", "GCF_000469605.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000469605.1/"),
    ("genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz", "GCF_000469605.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000469605.1/"),
    ("genomes/apis_mellifera/GCF_003254395.2_genomic.gff.gz", "GCF_003254395.2", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003254395.2/"),
    ("genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz", "GCF_003254395.2", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_003254395.2/"),
    ("genomes/apis_cerana/GCF_029169275.2_genomic.gff.gz", "GCF_029169275.2", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029169275.2/"),
    ("genomes/apis_cerana/GCF_029169275.2_protein.faa.gz", "GCF_029169275.2", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_029169275.2/"),
    ("genomes/apis_florea/GCF_048593485.1_genomic.gff.gz", "GCF_048593485.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048593485.1/"),
    ("genomes/apis_florea/GCF_048593485.1_protein.faa.gz", "GCF_048593485.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_048593485.1/"),
    ("genomes/bombus_terrestris/GCF_000214255.1_genomic.gff.gz", "GCF_000214255.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000214255.1/"),
    ("genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz", "GCF_000214255.1", "https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000214255.1/"),
    ("references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz", "GSM3757258 / GSE130963", "https://www.ncbi.nlm.nih.gov/gds/303757258"),
    ("references/transcriptome_2019/GSM3757258_AL.Readcount_FPKM.txt.gz", "GSM3757258 / GSE130963", "https://www.ncbi.nlm.nih.gov/gds/303757258"),
    ("references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz", "GSM3757259 / GSE130963", "https://www.ncbi.nlm.nih.gov/gds/303757259"),
    ("references/transcriptome_2019/GSM3757259_AD.Readcount_FPKM.txt.gz", "GSM3757259 / GSE130963", "https://www.ncbi.nlm.nih.gov/gds/303757259"),
    ("references/population_2023/extracted/Supplemental_Tables.xlsx", "Cao et al. 2023 supplement", "https://doi.org/10.1093/gbe/evad025"),
    ("references/rat_Nav1.4_NP_037310.2.faa", "NP_037310.2", "https://www.ncbi.nlm.nih.gov/protein/NP_037310.2"),
    ("references/uniprot_P15390.json", "UniProt P15390 topology", "https://rest.uniprot.org/uniprotkb/P15390.json"),
    ("references/bee_nav_functional_panel.faa", "compiled RefSeq/GenBank accessions", "https://www.ncbi.nlm.nih.gov/protein/"),
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_OUT,
        help="audit output path; defaults to ignored local_input_manifest.tsv",
    )
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output
    if output.resolve() == (ROOT / "data_sources.tsv").resolve():
        raise SystemExit(
            "Refusing to overwrite data_sources.tsv. Update the canonical snapshot "
            "only after reviewing upstream and result changes."
        )

    rows = []
    for relative_path, accession, source_url in SOURCES:
        path = ROOT / relative_path
        rows.append(
            {
                "relative_path": relative_path,
                "accession_or_dataset": accession,
                "source_url": source_url,
                "size_bytes": path.stat().st_size if path.exists() else "missing",
                "sha256": sha256(path) if path.exists() else "missing",
            }
        )
    output.parent.mkdir(parents=True, exist_ok=True)
    with output.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=rows[0].keys(),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)
    missing = [row["relative_path"] for row in rows if row["size_bytes"] == "missing"]
    print(f"Wrote {len(rows)} input records to {output}; missing={len(missing)}")
    for path in missing:
        print(f"  {path}")


if __name__ == "__main__":
    main()
