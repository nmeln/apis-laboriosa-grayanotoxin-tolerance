#!/usr/bin/env python3
"""Maintain checksum manifests for the comparative addendum."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent
WORK = ROOT / "work"

BASE_INPUTS = [
    "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
    "genomes/apis_laboriosa/refseq/GCF_014066325.1_protein.faa.gz",
    "genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz",
    "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
    "genomes/apis_dorsata/GCF_000469605.1_protein.faa.gz",
    "genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz",
    "genomes/apis_mellifera/GCF_003254395.2_genomic.gff.gz",
    "genomes/apis_mellifera/GCF_003254395.2_protein.faa.gz",
    "genomes/apis_cerana/GCF_029169275.2_genomic.gff.gz",
    "genomes/apis_cerana/GCF_029169275.2_protein.faa.gz",
    "genomes/apis_florea/GCF_048593485.1_genomic.gff.gz",
    "genomes/apis_florea/GCF_048593485.1_protein.faa.gz",
    "genomes/bombus_terrestris/GCF_000214255.1_protein.faa.gz",
    "references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz",
    "references/transcriptome_2019/GSM3757258_AL.Readcount_FPKM.txt.gz",
    "references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz",
    "references/transcriptome_2019/GSM3757259_AD.Readcount_FPKM.txt.gz",
    "references/rat_Nav1.4_NP_037310.2.faa",
]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def write_sha_manifest(paths: list[Path], output: Path) -> None:
    with output.open("w") as handle:
        for path in sorted(paths):
            if not path.is_file():
                raise FileNotFoundError(path)
            handle.write(f"{sha256(path)}  {path.relative_to(PROJECT)}\n")


def write_combined_input_manifest() -> int:
    with (PROJECT / "data_sources.tsv").open(newline="") as handle:
        base = {row["relative_path"]: row for row in csv.DictReader(handle, delimiter="\t")}
    with (ROOT / "input_sources.tsv").open(newline="") as handle:
        extra = list(csv.DictReader(handle, delimiter="\t"))
    rows = []
    for relative in BASE_INPUTS:
        row = base[relative]
        path = PROJECT / relative
        if path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
            raise RuntimeError(f"Pinned base input mismatch: {relative}")
        rows.append({"role": "base_repository_input", **row})
    for row in extra:
        path = ROOT / row["relative_path"]
        if path.stat().st_size != int(row["size_bytes"]) or sha256(path) != row["sha256"]:
            raise RuntimeError(f"Pinned addendum input mismatch: {row['relative_path']}")
        rows.append({"role": "addendum_input", **row})
    with (ROOT / "combined_input_manifest.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)
    return len(extra)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--refresh-results", action="store_true")
    parser.add_argument("--refresh-scripts", action="store_true")
    parser.add_argument("--refresh-work", action="store_true")
    parser.add_argument("--refresh-all", action="store_true")
    args = parser.parse_args()
    addendum_input_count = write_combined_input_manifest()

    primary = sorted((WORK / "primary_proteomes").glob("*.faa"))
    orthofinder = WORK / "primary_proteomes/OrthoFinder/Results_comparative_addendum_v1"
    orthofinder_key = [
        orthofinder / "Orthogroups/Orthogroups.tsv",
        orthofinder / "Orthogroups/Orthogroups_SingleCopyOrthologues.txt",
        orthofinder / "Species_Tree/SpeciesTree_rooted.txt",
    ]
    scripts = sorted((ROOT / "scripts").glob("*.py")) + [ROOT / "run_analysis.sh"]

    if args.refresh_all or args.refresh_work:
        write_sha_manifest(primary, ROOT / "primary_proteomes.sha256")
        write_sha_manifest(orthofinder_key, ROOT / "orthofinder_key_outputs.sha256")
    if args.refresh_all or args.refresh_scripts:
        write_sha_manifest(scripts, ROOT / "scripts.sha256")
    if args.refresh_all or args.refresh_results:
        results = [path for path in (ROOT / "results").rglob("*") if path.is_file()]
        write_sha_manifest(results, ROOT / "results.sha256")
    print(
        f"base_inputs={len(BASE_INPUTS)} addendum_inputs={addendum_input_count} "
        f"primary_proteomes={len(primary)} scripts={len(scripts)}"
    )


if __name__ == "__main__":
    main()
