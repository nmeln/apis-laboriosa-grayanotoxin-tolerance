#!/usr/bin/env python3
"""Screen gene-level annotations for gut/BBB barrier and clearance systems.

This is a copy-number screen, not an expression assay.  It counts parent genes rather
than protein isoforms and is most informative for large lineage-specific expansions.
"""

from __future__ import annotations

import csv
import gzip
import re
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
GFFS = {
    "Apis_laboriosa": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
    "Apis_dorsata": ROOT / "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
    "Apis_mellifera": ROOT / "genomes/apis_mellifera/GCF_003254395.2_genomic.gff.gz",
    "Apis_cerana": ROOT / "genomes/apis_cerana/GCF_029169275.2_genomic.gff.gz",
    "Apis_florea": ROOT / "genomes/apis_florea/GCF_048593485.1_genomic.gff.gz",
    "Bombus_terrestris": ROOT / "genomes/bombus_terrestris/GCF_000214255.1_genomic.gff.gz",
}

CATEGORIES = {
    "ABCB_P_glycoprotein": r"ATP-binding cassette sub-family B|P-glycoprotein|multidrug resistance protein(?!-associated)",
    "ABCC_MRP": r"ATP-binding cassette sub-family C|multidrug resistance-associated protein",
    "ABCG_half_transporter": r"ATP-binding cassette sub-family G|white protein|scarlet protein|brown protein",
    "organic_anion_transport": r"organic anion transport(?:er|ing polypeptide)",
    "major_facilitator": r"major facilitator superfamily",
    "aquaporin": r"aquaporin",
    "sodium_hydrogen_exchange": r"sodium.?hydrogen exchanger|Na\(\+\)/H\(\+\) exchange",
    "sodium_potassium_ATPase": r"sodium/potassium-transporting ATPase",
    "V_type_proton_ATPase": r"V-type proton ATPase|V-type sodium ATPase",
    "neurexin_IV": r"neurexin[- ]IV|neurexin 4",
    "contactin": r"\bcontactin\b",
    "neuroglian": r"\bneuroglian\b",
    "coracle": r"\bcoracle\b",
    "gliotactin": r"\bgliotactin\b",
    "mesh": r"\bmesh\b",
    "snakeskin": r"\bsnakeskin\b",
    "claudin_septate": r"claudin|megatrachea|\bsinuous\b",
    "chitin_synthase": r"chitin synthase",
    "peritrophic_matrix_chitin_binding": r"peritrophin|peritrophic matrix|chitin-binding protein",
    "mucin": r"\bmucin\b",
    "lipophorin": r"\blipophorin\b",
}
PATTERNS = {name: re.compile(pattern, re.I) for name, pattern in CATEGORIES.items()}


def parse_attrs(text):
    result = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def load_gene_products(path):
    labels = defaultdict(set)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = parse_attrs(fields[8])
            if attrs.get("Parent") and attrs.get("product"):
                labels[attrs["Parent"]].add(attrs["product"])
    return {
        gene: min(products, key=lambda label: (len(label), label.casefold(), label))
        for gene, products in labels.items()
    }


def main():
    products = {species: load_gene_products(path) for species, path in GFFS.items()}
    counts = []
    members = []
    for category, pattern in PATTERNS.items():
        row = {"category": category}
        for species, gene_products in products.items():
            found = sorted(
                [(gene, product) for gene, product in gene_products.items() if pattern.search(product)],
                key=lambda pair: (pair[1].lower(), pair[0]),
            )
            row[species] = len(found)
            for gene, product in found:
                members.append({"category": category, "species": species, "gene_parent": gene, "product": product})
        close_counts = [row[name] for name in ("Apis_dorsata", "Apis_mellifera", "Apis_cerana", "Apis_florea")]
        apis_max = max(close_counts)
        row["laboriosa_vs_dorsata_ratio"] = f"{row['Apis_laboriosa'] / row['Apis_dorsata']:.3f}" if row["Apis_dorsata"] else "NA"
        row["laboriosa_exceeds_other_apis_max"] = row["Apis_laboriosa"] > apis_max
        counts.append(row)

    with (OUT / "barrier_clearance_gene_counts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=counts[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(counts)
    with (OUT / "barrier_clearance_gene_members.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=members[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(members)

    for row in counts:
        print(row)


if __name__ == "__main__":
    main()
