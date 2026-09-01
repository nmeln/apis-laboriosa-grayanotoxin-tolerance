#!/usr/bin/env python3
"""Count annotated xenobiotic-handling gene families across bee genomes.

Counts are unique GFF parent genes, not transcript or protein isoforms.  They
are a screening tool, not proof of activity: annotation versions and gene-model
quality differ among assemblies.
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
OUT.mkdir(parents=True, exist_ok=True)

GFFS = {
    "Apis_laboriosa": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
    "Apis_dorsata": ROOT / "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
    "Apis_mellifera": ROOT / "genomes/apis_mellifera/GCF_003254395.2_genomic.gff.gz",
    "Apis_cerana": ROOT / "genomes/apis_cerana/GCF_029169275.2_genomic.gff.gz",
    "Apis_florea": ROOT / "genomes/apis_florea/GCF_048593485.1_genomic.gff.gz",
    "Bombus_terrestris": ROOT / "genomes/bombus_terrestris/GCF_000214255.1_genomic.gff.gz",
}

CATEGORIES = {
    "cytochrome_P450": re.compile(r"cytochrome P450", re.I),
    "UDP_glycosyltransferase": re.compile(r"UDP-(?:glycosyl|glucuronosyl)transferase", re.I),
    "glutathione_S_transferase": re.compile(r"glutathione S-transferase", re.I),
    "ABC_or_multidrug_transporter": re.compile(
        r"ATP-binding cassette|ABC transporter|multidrug resistance(?:-associated)? protein|P-glycoprotein",
        re.I,
    ),
    "sulfotransferase": re.compile(r"sulfotransferase", re.I),
    "carboxylesterase_or_esterase": re.compile(r"carboxylesterase|(?<!chol)esterase", re.I),
    "epoxide_hydrolase": re.compile(r"epoxide hydrolase", re.I),
    "aldo_keto_reductase": re.compile(r"aldo-keto reductase", re.I),
    "short_chain_dehydrogenase": re.compile(r"short-chain dehydrogenase|short chain dehydrogenase", re.I),
    "organic_anion_transporter": re.compile(r"organic anion transporter", re.I),
    "major_facilitator_transporter": re.compile(r"major facilitator superfamily", re.I),
    "aquaporin": re.compile(r"aquaporin", re.I),
    "lipocalin": re.compile(r"lipocalin", re.I),
    "lipophorin": re.compile(r"lipophorin", re.I),
}


def attrs(text: str):
    result = {}
    for item in text.rstrip().split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = unquote(value)
    return result


def load_products(path: Path):
    """Return one representative product string per annotated parent gene."""
    products = defaultdict(set)
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            at = attrs(fields[8])
            parent = at.get("Parent")
            product = at.get("product")
            if parent and product:
                products[parent].add(product)
    # Prefer the shortest product label, which usually omits isoform suffixes.
    return {
        gene: min(labels, key=lambda label: (len(label), label.casefold(), label))
        for gene, labels in products.items()
    }


def main():
    by_species = {species: load_products(path) for species, path in GFFS.items()}

    count_rows = []
    detail_rows = []
    for category, pattern in CATEGORIES.items():
        row = {"category": category}
        for species, gene_products in by_species.items():
            matches = sorted(
                ((gene, product) for gene, product in gene_products.items() if pattern.search(product)),
                key=lambda x: (x[1].lower(), x[0]),
            )
            row[species] = len(matches)
            for gene, product in matches:
                detail_rows.append(
                    {
                        "category": category,
                        "species": species,
                        "gene_parent": gene,
                        "product": product,
                    }
                )
        count_rows.append(row)

    with (OUT / "detox_family_counts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=count_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(count_rows)

    with (OUT / "detox_family_members.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=detail_rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(detail_rows)

    for row in count_rows:
        print("\t".join(str(row[key]) for key in row))


if __name__ == "__main__":
    main()
