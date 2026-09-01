#!/usr/bin/env python3
"""Re-run the repository annotation screens with the current Bombus reference."""

from __future__ import annotations

import csv
import importlib.util
from pathlib import Path


ADDENDUM_ROOT = Path(__file__).resolve().parents[1]
PROJECT_ROOT = ADDENDUM_ROOT.parent
RESULTS = ADDENDUM_ROOT / "results"
CURRENT_BOMBUS_GFF = (
    ADDENDUM_ROOT
    / "inputs/bombus_terrestris_current/GCF_910591885.1_iyBomTerr1.2_genomic.gff.gz"
)


def load_module(name: str, path: Path):
    spec = importlib.util.spec_from_file_location(name, path)
    if spec is None or spec.loader is None:
        raise RuntimeError(f"Cannot load {path}")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def write_table(path: Path, rows: list[dict[str, object]]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(rows)


def run_detox() -> None:
    module = load_module("repo_detox", PROJECT_ROOT / "scripts/analyze_detox_families.py")
    gffs = dict(module.GFFS)
    gffs["Bombus_terrestris"] = CURRENT_BOMBUS_GFF
    products = {species: module.load_products(path) for species, path in gffs.items()}
    counts: list[dict[str, object]] = []
    members: list[dict[str, object]] = []
    for category, pattern in module.CATEGORIES.items():
        row: dict[str, object] = {"category": category}
        for species, gene_products in products.items():
            found = sorted(
                ((gene, product) for gene, product in gene_products.items() if pattern.search(product)),
                key=lambda pair: (pair[1].lower(), pair[0]),
            )
            row[species] = len(found)
            members.extend(
                {"category": category, "species": species, "gene_parent": gene, "product": product}
                for gene, product in found
            )
        counts.append(row)
    RESULTS.mkdir(parents=True, exist_ok=True)
    write_table(RESULTS / "detox_family_counts_current.tsv", counts)
    write_table(RESULTS / "detox_family_members_current.tsv", members)


def run_barrier() -> None:
    module = load_module("repo_barrier", PROJECT_ROOT / "scripts/analyze_barrier_clearance.py")
    gffs = dict(module.GFFS)
    gffs["Bombus_terrestris"] = CURRENT_BOMBUS_GFF
    products = {species: module.load_gene_products(path) for species, path in gffs.items()}
    counts: list[dict[str, object]] = []
    members: list[dict[str, object]] = []
    for category, pattern in module.PATTERNS.items():
        row: dict[str, object] = {"category": category}
        for species, gene_products in products.items():
            found = sorted(
                ((gene, product) for gene, product in gene_products.items() if pattern.search(product)),
                key=lambda pair: (pair[1].lower(), pair[0]),
            )
            row[species] = len(found)
            members.extend(
                {"category": category, "species": species, "gene_parent": gene, "product": product}
                for gene, product in found
            )
        controls = [int(row[species]) for species in ("Apis_dorsata", "Apis_mellifera", "Apis_cerana", "Apis_florea")]
        row["laboriosa_vs_dorsata_ratio"] = (
            f"{int(row['Apis_laboriosa']) / int(row['Apis_dorsata']):.3f}"
            if int(row["Apis_dorsata"])
            else "NA"
        )
        row["laboriosa_exceeds_other_apis_max"] = int(row["Apis_laboriosa"]) > max(controls)
        counts.append(row)
    write_table(RESULTS / "barrier_clearance_gene_counts_current.tsv", counts)
    write_table(RESULTS / "barrier_clearance_gene_members_current.tsv", members)


def main() -> None:
    run_detox()
    run_barrier()


if __name__ == "__main__":
    main()
