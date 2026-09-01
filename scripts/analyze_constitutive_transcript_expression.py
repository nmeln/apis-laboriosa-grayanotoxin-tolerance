#!/usr/bin/env python3
"""Exploratory whole-worker expression screen for barrier/detox hypotheses.

The public 2019 transcriptomes contain one pooled, untreated whole-worker sample
per species.  They cannot establish differential expression, but they can test
for a gross constitutive family-wide signal.  Species-specific exact probes map
RefSeq genes to Trinity unigenes, after which the submitted FPKM values are
summarized by annotation category.
"""

from __future__ import annotations

import csv
import gzip
import re
import statistics
from bisect import bisect_right
from collections import defaultdict
from pathlib import Path
from urllib.parse import unquote


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results"
K = 31
MAX_SAMPLED_PER_TRANSCRIPT = 48
MAX_UNIQUE_PROBES_PER_GENE = 120
MIN_HITS = 3

SPECIES = {
    "Apis_laboriosa": {
        "gff": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_genomic.gff.gz",
        "rna": ROOT / "genomes/apis_laboriosa/refseq/GCF_014066325.1_rna.fna.gz",
        "unigenes": ROOT / "references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz",
        "fpkm": ROOT / "references/transcriptome_2019/GSM3757258_AL.Readcount_FPKM.txt.gz",
    },
    "Apis_dorsata": {
        "gff": ROOT / "genomes/apis_dorsata/GCF_000469605.1_genomic.gff.gz",
        "rna": ROOT / "genomes/apis_dorsata/GCF_000469605.1_rna.fna.gz",
        "unigenes": ROOT / "references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz",
        "fpkm": ROOT / "references/transcriptome_2019/GSM3757259_AD.Readcount_FPKM.txt.gz",
    },
}

CATEGORIES = {
    "cytochrome_P450": r"cytochrome P450",
    "UDP_glycosyltransferase": r"UDP[- ](?:glycosyl|glucuronosyl)transferase",
    "glutathione_S_transferase": r"glutathione S-transferase",
    "ABCB_P_glycoprotein": r"ATP-binding cassette sub-family B|P-glycoprotein|multidrug resistance protein(?!-associated)",
    "ABCC_MRP": r"ATP-binding cassette sub-family C|multidrug resistance-associated protein",
    "ABCG_half_transporter": r"ATP-binding cassette sub-family G|white protein|scarlet protein|brown protein",
    "organic_anion_transport": r"organic anion transport(?:er|ing polypeptide)",
    "major_facilitator": r"major facilitator superfamily",
    "aquaporin": r"aquaporin",
    "V_type_proton_ATPase": r"V-type proton ATPase|V-type sodium ATPase",
    "chitin_synthase": r"chitin synthase",
    "peritrophic_matrix_chitin_binding": r"peritrophin|peritrophic matrix|chitin-binding protein",
    "mucin": r"\bmucin\b",
    "Para": r"sodium (?:voltage-gated channel paralytic|channel protein (?:para|paralytic))",
    "DSC1_60E": r"sodium channel protein 60E",
    "Nach_DEG_ENaC": r"sodium channel protein Nach|pickpocket protein 28",
    "COMMD3": r"COMM domain-containing protein 3",
}
PATTERNS = {name: re.compile(pattern, re.I) for name, pattern in CATEGORIES.items()}


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


def revcomp(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def sample_probes(sequence: str):
    if len(sequence) < K:
        return set()
    count = min(MAX_SAMPLED_PER_TRANSCRIPT, len(sequence) - K + 1)
    if count == 1:
        starts = [0]
    else:
        starts = [round(i * (len(sequence) - K) / (count - 1)) for i in range(count)]
    return {
        sequence[start : start + K]
        for start in starts
        if "N" not in sequence[start : start + K] and len(set(sequence[start : start + K])) >= 3
    }


def read_expression(path: Path):
    values = {}
    with gzip.open(path, "rt") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            values[row["gene_id"]] = (float(row["Read_count"]), float(row["FPKM"]))
    sorted_fpkms = sorted(fpkm for _read_count, fpkm in values.values())
    return values, sum(read_count for read_count, _fpkm in values.values()), sorted_fpkms


def analyze_species(species: str, paths: dict[str, Path]):
    transcript_meta = {}
    gene_products = defaultdict(set)
    gene_categories = defaultdict(set)
    with gzip.open(paths["gff"], "rt") as handle:
        for line in handle:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] not in {"mRNA", "transcript"}:
                continue
            attrs = parse_attrs(fields[8])
            parent = attrs.get("Parent")
            transcript = attrs.get("transcript_id", attrs.get("Name"))
            product = attrs.get("product")
            if not parent or not transcript or not product:
                continue
            categories = {name for name, pattern in PATTERNS.items() if pattern.search(product)}
            if not categories:
                continue
            transcript_meta[transcript] = (parent, product)
            gene_products[parent].add(product)
            gene_categories[parent].update(categories)

    raw_gene_probes = defaultdict(set)
    for header, sequence in fasta_records(paths["rna"]):
        transcript = header.split()[0]
        if transcript in transcript_meta:
            gene, _product = transcript_meta[transcript]
            raw_gene_probes[gene].update(sample_probes(sequence))

    owners = defaultdict(set)
    for gene, probes in raw_gene_probes.items():
        for probe in probes:
            owners[min(probe, revcomp(probe))].add(gene)
    gene_probes = defaultdict(set)
    for canonical, genes in owners.items():
        if len(genes) == 1:
            gene_probes[next(iter(genes))].add(canonical)
    for gene, probes in list(gene_probes.items()):
        if len(probes) > MAX_UNIQUE_PROBES_PER_GENE:
            ordered = sorted(probes)
            indexes = [
                round(i * (len(ordered) - 1) / (MAX_UNIQUE_PROBES_PER_GENE - 1))
                for i in range(MAX_UNIQUE_PROBES_PER_GENE)
            ]
            gene_probes[gene] = {ordered[index] for index in indexes}

    lookup = defaultdict(list)
    for gene, probes in gene_probes.items():
        for canonical in probes:
            lookup[canonical].append((gene, canonical))
            reverse = revcomp(canonical)
            if reverse != canonical:
                lookup[reverse].append((gene, canonical))

    best = {}
    for header, sequence in fasta_records(paths["unigenes"]):
        unigene = header.split()[0]
        local = defaultdict(set)
        for position in range(len(sequence) - K + 1):
            for gene, canonical in lookup.get(sequence[position : position + K], ()):
                local[gene].add(canonical)
        for gene, matched in local.items():
            candidate = (len(matched), unigene)
            if gene not in best or candidate > (len(best[gene][1]), best[gene][0]):
                best[gene] = (unigene, matched)

    expression, total_read_count, sorted_fpkms = read_expression(paths["fpkm"])
    rows = []
    for gene in sorted(gene_categories):
        unigene, matched = best.get(gene, ("", set()))
        mapped = len(matched) >= MIN_HITS
        read_count, fpkm = expression.get(unigene, (0.0, 0.0)) if mapped else (0.0, 0.0)
        for category in sorted(gene_categories[gene]):
            rows.append(
                {
                    "species": species,
                    "category": category,
                    "gene_parent": gene,
                    "product": min(
                        gene_products[gene],
                        key=lambda label: (len(label), label.casefold(), label),
                    ),
                    "best_unigene": unigene if mapped else "",
                    "exact_probe_hits": len(matched),
                    "unique_gene_probes": len(gene_probes.get(gene, set())),
                    "probe_fraction": f"{len(matched) / len(gene_probes[gene]):.6f}" if gene_probes.get(gene) else "NA",
                    "mapped": mapped,
                    "read_count": f"{read_count:.2f}",
                    "cpm": f"{read_count / total_read_count * 1_000_000:.6f}" if total_read_count else "NA",
                    "fpkm": f"{fpkm:.6f}",
                    "fpkm_percentile": f"{bisect_right(sorted_fpkms, fpkm) / len(sorted_fpkms) * 100:.6f}" if mapped else "NA",
                    "fpkm_rank_desc": len(sorted_fpkms) - bisect_right(sorted_fpkms, fpkm) + 1 if mapped else "NA",
                }
            )
    return rows


def main() -> None:
    rows = []
    for species, paths in SPECIES.items():
        print(f"Mapping {species} target genes to its pooled whole-worker transcriptome...")
        rows.extend(analyze_species(species, paths))

    with (OUT / "constitutive_expression_gene_mappings.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=rows[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    summary = []
    for category in CATEGORIES:
        for species in SPECIES:
            selected = [row for row in rows if row["species"] == species and row["category"] == category]
            mapped = [row for row in selected if row["mapped"]]
            fpkms = [float(row["fpkm"]) for row in mapped]
            # A Trinity unigene can attract probes from more than one closely
            # related annotation model.  Deduplicate aggregate expression by
            # unigene so family sums are not artificially inflated.
            unique_by_unigene = {}
            for row in mapped:
                unigene = row["best_unigene"]
                if unigene not in unique_by_unigene or int(row["exact_probe_hits"]) > int(unique_by_unigene[unigene]["exact_probe_hits"]):
                    unique_by_unigene[unigene] = row
            summary.append(
                {
                    "category": category,
                    "species": species,
                    "annotated_genes": len(selected),
                    "mapped_genes": len(mapped),
                    "mapped_unique_unigenes": len(unique_by_unigene),
                    "mapping_fraction": f"{len(mapped) / len(selected):.6f}" if selected else "NA",
                    "sum_fpkm_unique": f"{sum(float(row['fpkm']) for row in unique_by_unigene.values()):.6f}",
                    "sum_cpm_unique": f"{sum(float(row['cpm']) for row in unique_by_unigene.values()):.6f}",
                    "median_fpkm": f"{statistics.median(fpkms):.6f}" if fpkms else "NA",
                    "max_fpkm": f"{max(fpkms):.6f}" if fpkms else "NA",
                }
            )
    with (OUT / "constitutive_expression_family_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(summary)

    contrasts = []
    for category in CATEGORIES:
        laboriosa = next(
            row for row in summary if row["category"] == category and row["species"] == "Apis_laboriosa"
        )
        dorsata = next(
            row for row in summary if row["category"] == category and row["species"] == "Apis_dorsata"
        )
        lab_fpkm, dor_fpkm = float(laboriosa["sum_fpkm_unique"]), float(dorsata["sum_fpkm_unique"])
        lab_cpm, dor_cpm = float(laboriosa["sum_cpm_unique"]), float(dorsata["sum_cpm_unique"])
        contrasts.append(
            {
                "category": category,
                "laboriosa_sum_fpkm_unique": f"{lab_fpkm:.6f}",
                "dorsata_sum_fpkm_unique": f"{dor_fpkm:.6f}",
                "laboriosa_to_dorsata_fpkm_ratio": f"{lab_fpkm / dor_fpkm:.6f}" if dor_fpkm else "NA",
                "laboriosa_sum_cpm_unique": f"{lab_cpm:.6f}",
                "dorsata_sum_cpm_unique": f"{dor_cpm:.6f}",
                "laboriosa_to_dorsata_cpm_ratio": f"{lab_cpm / dor_cpm:.6f}" if dor_cpm else "NA",
                "design_warning": "one pooled untreated whole-worker library per species; descriptive only",
            }
        )
    with (OUT / "constitutive_expression_contrasts.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=contrasts[0].keys(), delimiter="\t")
        writer.writeheader()
        writer.writerows(contrasts)

    for category in CATEGORIES:
        values = [row for row in summary if row["category"] == category]
        print(
            category
            + ": "
            + "; ".join(
                f"{row['species']} {row['mapped_genes']}/{row['annotated_genes']} genes, "
                f"unique-transcript FPKM {row['sum_fpkm_unique']}"
                for row in values
            )
        )


if __name__ == "__main__":
    main()
