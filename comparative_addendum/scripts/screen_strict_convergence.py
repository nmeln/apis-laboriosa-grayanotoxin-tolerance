#!/usr/bin/env python3
"""Exploratory six-bee screen for strict amino-acid sharing and copy number.

This deliberately narrow analysis asks whether Apis laboriosa and the directly
phenotyped grayanotoxin-tolerant Bombus terrestris share protein states that are
absent from four other Apis references.  It also runs the identical screen with
each other Apis species as a negative-control focal species.

The site pattern is descriptive.  With these taxa alone it cannot distinguish
parallel substitution from reversal to an ancestral state, and the focal
A. laboriosa phenotype has not been established in a controlled challenge.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from pyfamsa import Aligner, Sequence
from scipy.stats import fisher_exact


SPECIES = [
    "Apis_laboriosa",
    "Apis_dorsata",
    "Apis_mellifera",
    "Apis_cerana",
    "Apis_florea",
    "Bombus_terrestris",
]
APIS = SPECIES[:-1]
BOMBUS = "Bombus_terrestris"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")

CORE_CATEGORIES = {
    "ABCB_P_glycoprotein",
    "ABCC_MRP",
    "ABCG_half_transporter",
    "V_type_proton_ATPase",
    "aquaporin",
    "chitin_synthase",
    "claudin_septate",
    "contactin",
    "coracle",
    "lipophorin",
    "major_facilitator",
    "major_facilitator_transporter",
    "mesh",
    "mucin",
    "neuroglian",
    "organic_anion_transport",
    "organic_anion_transporter",
    "peritrophic_matrix_chitin_binding",
    "snakeskin",
    "sodium_hydrogen_exchange",
}


def read_fasta(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    accession: str | None = None
    chunks: list[str] = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = "".join(chunks)
                accession = line[1:].split(maxsplit=1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if accession is not None:
        records[accession] = "".join(chunks)
    return records


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=fields,
            delimiter="\t",
            extrasaction="ignore",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(rows)


def normalize_gene(gene: str) -> str:
    return gene.removeprefix("gene-")


def load_candidate_categories(project: Path, session_root: Path) -> dict[tuple[str, str], set[str]]:
    categories: dict[tuple[str, str], set[str]] = defaultdict(set)
    paths = [
        project / "results/barrier_clearance_gene_members.tsv",
        project / "results/detox_family_members.tsv",
        session_root / "results/barrier_clearance_gene_members_current.tsv",
        session_root / "results/detox_family_members_current.tsv",
    ]
    for path in paths:
        if not path.exists():
            continue
        for row in read_tsv(path):
            key = (row["species"], normalize_gene(row["gene_parent"]))
            categories[key].add(row["category"])
    return categories


def parse_orthogroups(
    path: Path, accession_species: dict[str, str] | None = None
) -> tuple[list[str], list[dict[str, list[str]]]]:
    groups: list[dict[str, list[str]]] = []
    if path.suffix == ".txt":
        if accession_species is None:
            raise RuntimeError("Accession-to-species mapping is required for Orthogroups.txt")
        with path.open() as handle:
            for line in handle:
                name, member_text = line.rstrip().split(":", 1)
                parsed: dict[str, list[str]] = {"Orthogroup": [name]}
                parsed.update({species: [] for species in SPECIES})
                for accession in member_text.split():
                    species = accession_species.get(accession)
                    if species is None:
                        raise RuntimeError(f"Unknown accession in {name}: {accession}")
                    parsed[species].append(accession)
                groups.append(parsed)
        return SPECIES, groups

    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        species_columns = [x for x in reader.fieldnames or [] if x != "Orthogroup"]
        missing = set(SPECIES) - set(species_columns)
        if missing:
            raise RuntimeError(f"Orthogroups table lacks species: {sorted(missing)}")
        for row in reader:
            parsed: dict[str, list[str]] = {"Orthogroup": [row["Orthogroup"]]}
            for species in SPECIES:
                parsed[species] = [x.strip() for x in row[species].split(",") if x.strip()]
            groups.append(parsed)
    return species_columns, groups


def position_maps(aligned: dict[str, str]) -> dict[str, list[int | None]]:
    maps: dict[str, list[int | None]] = {}
    for species, sequence in aligned.items():
        position = 0
        values: list[int | None] = []
        for aa in sequence:
            if aa == "-":
                values.append(None)
            else:
                position += 1
                values.append(position)
        maps[species] = values
    return maps


def bh_adjust(pvalues: list[float]) -> list[float]:
    if not pvalues:
        return []
    order = sorted(range(len(pvalues)), key=pvalues.__getitem__)
    adjusted = [1.0] * len(pvalues)
    running = 1.0
    for reverse_rank, index in enumerate(reversed(order), start=1):
        rank = len(pvalues) - reverse_rank + 1
        running = min(running, pvalues[index] * len(pvalues) / rank)
        adjusted[index] = min(1.0, running)
    return adjusted


def matched_permutation(
    summaries: list[dict[str, object]],
    candidate_key: str,
    hit_key: str,
    divergence_key: str,
    permutations: int = 20_000,
    seed: int = 240801,
) -> tuple[float, float, int]:
    """Shuffle labels within callable-length and divergence strata.

    The statistic is candidate minus background mean strict-hit density per
    1,000 callable sites.  Strata reduce, but do not eliminate, rate and length
    confounding.  This is an exploratory calibration rather than a phylogenetic
    association test.
    """

    eligible = [row for row in summaries if int(row["callable_sites"]) > 0]
    if not eligible:
        return math.nan, math.nan, 0

    lengths = sorted(int(row["callable_sites"]) for row in eligible)
    divergences = sorted(float(row[divergence_key]) for row in eligible)

    def quantile(values: list[float], fraction: float) -> float:
        index = round((len(values) - 1) * fraction)
        return values[index]

    length_cuts = [quantile(lengths, q) for q in (0.25, 0.5, 0.75)]
    divergence_cuts = [quantile(divergences, q) for q in (0.25, 0.5, 0.75)]

    def bin_of(value: float, cuts: list[float]) -> int:
        return sum(value > cut for cut in cuts)

    strata: dict[tuple[int, int], list[int]] = defaultdict(list)
    labels: list[bool] = []
    hits: list[int] = []
    sites: list[int] = []
    for index, row in enumerate(eligible):
        strata[(
            bin_of(float(row["callable_sites"]), length_cuts),
            bin_of(float(row[divergence_key]), divergence_cuts),
        )].append(index)
        labels.append(bool(row[candidate_key]))
        hits.append(int(row[hit_key]))
        sites.append(int(row["callable_sites"]))

    hit_array = np.asarray(hits, dtype=np.int64)
    site_array = np.asarray(sites, dtype=np.int64)
    label_array = np.asarray(labels, dtype=bool)
    total_hits = int(hit_array.sum())
    total_sites = int(site_array.sum())

    def statistic(candidate_hits: int, candidate_sites: int) -> float:
        background_hits = total_hits - candidate_hits
        background_sites = total_sites - candidate_sites
        if candidate_sites == 0 or background_sites == 0:
            return math.nan
        return 1000 * (candidate_hits / candidate_sites - background_hits / background_sites)

    observed_candidate_hits = int(hit_array[label_array].sum())
    observed_candidate_sites = int(site_array[label_array].sum())
    observed = statistic(observed_candidate_hits, observed_candidate_sites)
    rng = np.random.default_rng(seed)
    stratum_arrays = []
    for values in strata.values():
        indices = np.asarray(values, dtype=np.int64)
        stratum_arrays.append((indices, int(label_array[indices].sum())))
    exceed = 0
    for _ in range(permutations):
        candidate_hits = 0
        candidate_sites = 0
        for indices, count in stratum_arrays:
            if count == 0:
                continue
            chosen = indices if count == len(indices) else rng.choice(indices, size=count, replace=False)
            candidate_hits += int(hit_array[chosen].sum())
            candidate_sites += int(site_array[chosen].sum())
        permuted = statistic(candidate_hits, candidate_sites)
        if abs(permuted) >= abs(observed):
            exceed += 1
    return observed, (exceed + 1) / (permutations + 1), len(eligible)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthofinder", type=Path, required=True)
    parser.add_argument("--session-root", type=Path, default=Path(__file__).resolve().parents[1])
    parser.add_argument("--project", type=Path, default=Path(__file__).resolve().parents[2])
    parser.add_argument("--proteomes", type=Path)
    parser.add_argument("--members", type=Path)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--threads", type=int, default=4)
    args = parser.parse_args()

    output = args.output or args.session_root / "results"
    output.mkdir(parents=True, exist_ok=True)

    proteomes = args.proteomes or args.session_root / "primary_proteomes"
    members_path = args.members or args.session_root / "primary_proteome_members.tsv"
    members = read_tsv(members_path)
    metadata = {row["accession"]: row for row in members}
    accession_species = {row["accession"]: row["species"] for row in members}
    orthogroups_path = args.orthofinder / "Orthogroups" / "Orthogroups.tsv"
    if not orthogroups_path.exists():
        orthogroups_path = args.orthofinder / "Orthogroups" / "Orthogroups.txt"
    if not orthogroups_path.exists():
        raise FileNotFoundError(orthogroups_path)
    sequences: dict[str, dict[str, str]] = {
        species: read_fasta(proteomes / f"{species}.faa")
        for species in SPECIES
    }
    categories_by_gene = load_candidate_categories(args.project, args.session_root)
    _, orthogroups = parse_orthogroups(orthogroups_path, accession_species)

    aligner = Aligner(threads=max(1, args.threads))
    orthofinder_alignments = args.orthofinder / "MultipleSequenceAlignments"
    use_saved_alignments = orthofinder_alignments.is_dir()
    site_hits: list[dict[str, object]] = []
    gene_summaries: list[dict[str, object]] = []
    single_copy_count = 0

    for group in orthogroups:
        if not all(len(group[species]) == 1 for species in SPECIES):
            continue
        single_copy_count += 1
        og = group["Orthogroup"][0]
        accessions = {species: group[species][0] for species in SPECIES}
        saved_alignment = orthofinder_alignments / f"{og}.fa"
        if use_saved_alignments:
            if not saved_alignment.exists():
                raise FileNotFoundError(saved_alignment)
            raw_alignment = read_fasta(saved_alignment)
            aligned = {}
            for species in SPECIES:
                expected_header = f"{species}_{accessions[species]}"
                if expected_header not in raw_alignment:
                    raise RuntimeError(
                        f"Saved alignment {saved_alignment} lacks {expected_header}"
                    )
                aligned[species] = raw_alignment[expected_header]
        else:
            aligned_records = aligner.align([
                Sequence(species.encode(), sequences[species][accessions[species]].encode())
                for species in SPECIES
            ])
            aligned = {record.id.decode(): record.sequence.decode() for record in aligned_records}
        if set(aligned) != set(SPECIES):
            raise RuntimeError(f"Alignment species mismatch in {og}")
        maps = position_maps(aligned)

        category_union: set[str] = set()
        for species, accession in accessions.items():
            gene = metadata[accession]["gene"]
            category_union.update(categories_by_gene.get((species, gene), set()))
        description = metadata[accessions["Apis_laboriosa"]]["description"]
        is_para = "sodium channel protein para" in description.lower() or "sodium voltage-gated channel paralytic" in description.lower()
        if is_para:
            category_union.add("Para_voltage_gated_sodium_channel")
        is_barrier_or_detox = bool(category_union - {"Para_voltage_gated_sodium_channel"})
        is_core = bool(category_union & CORE_CATEGORIES)

        hit_counts = Counter()
        callable_sites = 0
        bombus_differences = Counter()
        length = len(next(iter(aligned.values())))
        for column in range(length):
            residues = {species: aligned[species][column] for species in SPECIES}
            if not all(aa in VALID_AA for aa in residues.values()):
                continue
            callable_sites += 1
            for species in APIS:
                if residues[BOMBUS] != residues[species]:
                    bombus_differences[species] += 1
            for focal in APIS:
                controls = [species for species in APIS if species != focal]
                control_states = {residues[species] for species in controls}
                if len(control_states) != 1:
                    continue
                control_state = next(iter(control_states))
                if residues[focal] == residues[BOMBUS] and residues[focal] != control_state:
                    hit_counts[focal] += 1
                    site_hits.append({
                        "orthogroup": og,
                        "alignment_column_1based": column + 1,
                        "focal_apis": focal,
                        "focal_state": residues[focal],
                        "other_apis_state": control_state,
                        "bombus_state": residues[BOMBUS],
                        "focal_position_1based": maps[focal][column],
                        "laboriosa_position_1based": maps["Apis_laboriosa"][column],
                        "bombus_position_1based": maps[BOMBUS][column],
                        "candidate_any": is_barrier_or_detox,
                        "candidate_core": is_core,
                        "candidate_categories": ";".join(sorted(category_union)),
                        "laboriosa_accession": accessions["Apis_laboriosa"],
                        "laboriosa_gene": metadata[accessions["Apis_laboriosa"]]["gene"],
                        "laboriosa_description": description,
                    })

        summary: dict[str, object] = {
            "orthogroup": og,
            "callable_sites": callable_sites,
            "bombus_difference_fraction": bombus_differences["Apis_laboriosa"] / callable_sites if callable_sites else math.nan,
            "candidate_any": is_barrier_or_detox,
            "candidate_core": is_core,
            "is_para": is_para,
            "candidate_categories": ";".join(sorted(category_union)),
            "laboriosa_accession": accessions["Apis_laboriosa"],
            "laboriosa_gene": metadata[accessions["Apis_laboriosa"]]["gene"],
            "laboriosa_description": description,
        }
        for focal in APIS:
            summary[f"hits_{focal}"] = hit_counts[focal]
            summary[f"bombus_difference_fraction_{focal}"] = (
                bombus_differences[focal] / callable_sites if callable_sites else math.nan
            )
        gene_summaries.append(summary)

    site_fields = [
        "orthogroup", "alignment_column_1based", "focal_apis", "focal_state",
        "other_apis_state", "bombus_state", "focal_position_1based",
        "laboriosa_position_1based", "bombus_position_1based", "candidate_any",
        "candidate_core", "candidate_categories", "laboriosa_accession",
        "laboriosa_gene", "laboriosa_description",
    ]
    write_tsv(output / "strict_site_hits.tsv", site_hits, site_fields)

    summary_fields = [
        "orthogroup", "callable_sites", "bombus_difference_fraction", "candidate_any",
        "candidate_core", "is_para", "candidate_categories", "laboriosa_accession",
        "laboriosa_gene", "laboriosa_description",
    ] + [f"hits_{species}" for species in APIS] + [f"bombus_difference_fraction_{species}" for species in APIS]
    write_tsv(output / "orthogroup_site_summary.tsv", gene_summaries, summary_fields)

    aggregates: list[dict[str, object]] = []
    total_callable = sum(int(row["callable_sites"]) for row in gene_summaries)
    for focal in APIS:
        key = f"hits_{focal}"
        hits = sum(int(row[key]) for row in gene_summaries)
        genes = sum(int(row[key]) > 0 for row in gene_summaries)
        aggregates.append({
            "focal_apis": focal,
            "single_copy_orthogroups": len(gene_summaries),
            "callable_sites": total_callable,
            "strict_sites": hits,
            "strict_sites_per_million_callable": hits * 1_000_000 / total_callable if total_callable else math.nan,
            "orthogroups_with_strict_site": genes,
        })
    write_tsv(
        output / "strict_site_aggregate.tsv",
        aggregates,
        ["focal_apis", "single_copy_orthogroups", "callable_sites", "strict_sites", "strict_sites_per_million_callable", "orthogroups_with_strict_site"],
    )

    enrichment_rows: list[dict[str, object]] = []
    for focal in APIS:
        hit_key = f"hits_{focal}"
        for label_key, label_name in (("candidate_any", "barrier_or_detox"), ("candidate_core", "core_toxicokinetic")):
            candidate = [row for row in gene_summaries if bool(row[label_key])]
            background = [row for row in gene_summaries if not bool(row[label_key])]
            candidate_hits = sum(int(row[hit_key]) for row in candidate)
            background_hits = sum(int(row[hit_key]) for row in background)
            candidate_sites = sum(int(row["callable_sites"]) for row in candidate)
            background_sites = sum(int(row["callable_sites"]) for row in background)
            site_table = [[candidate_hits, candidate_sites - candidate_hits], [background_hits, background_sites - background_hits]]
            site_or, site_p = fisher_exact(site_table, alternative="two-sided")
            candidate_genes_hit = sum(int(row[hit_key]) > 0 for row in candidate)
            background_genes_hit = sum(int(row[hit_key]) > 0 for row in background)
            gene_table = [
                [candidate_genes_hit, len(candidate) - candidate_genes_hit],
                [background_genes_hit, len(background) - background_genes_hit],
            ]
            gene_or, gene_p = fisher_exact(gene_table, alternative="two-sided")
            observed, permutation_p, permuted_genes = matched_permutation(
                gene_summaries,
                label_key,
                hit_key,
                f"bombus_difference_fraction_{focal}",
            )
            enrichment_rows.append({
                "focal_apis": focal,
                "candidate_set": label_name,
                "candidate_orthogroups": len(candidate),
                "background_orthogroups": len(background),
                "candidate_callable_sites": candidate_sites,
                "background_callable_sites": background_sites,
                "candidate_strict_sites": candidate_hits,
                "background_strict_sites": background_hits,
                "candidate_orthogroups_with_hit": candidate_genes_hit,
                "background_orthogroups_with_hit": background_genes_hit,
                "site_fisher_odds_ratio": site_or,
                "site_fisher_p_two_sided": site_p,
                "gene_fisher_odds_ratio": gene_or,
                "gene_fisher_p_two_sided": gene_p,
                "matched_density_difference_per_1000_sites": observed,
                "matched_permutation_p_two_sided": permutation_p,
                "permutation_orthogroups": permuted_genes,
                "permutations": 20_000,
            })
    pvalues = [float(row["matched_permutation_p_two_sided"]) for row in enrichment_rows]
    for row, adjusted in zip(enrichment_rows, bh_adjust(pvalues)):
        row["matched_permutation_bh_q"] = adjusted
    enrichment_fields = list(enrichment_rows[0])
    write_tsv(output / "candidate_site_enrichment.tsv", enrichment_rows, enrichment_fields)

    copy_hits: list[dict[str, object]] = []
    copy_aggregate = Counter()
    for group in orthogroups:
        og = group["Orthogroup"][0]
        counts = {species: len(group[species]) for species in SPECIES}
        group_categories: set[str] = set()
        laboriosa_descriptions: list[str] = []
        for species in SPECIES:
            for accession in group[species]:
                gene = metadata.get(accession, {}).get("gene")
                if gene:
                    group_categories.update(categories_by_gene.get((species, gene), set()))
                if species == "Apis_laboriosa" and accession in metadata:
                    laboriosa_descriptions.append(metadata[accession]["description"])
        is_candidate = bool(group_categories)
        for focal in APIS:
            controls = [species for species in APIS if species != focal]
            control_max = max(counts[species] for species in controls)
            strict = counts[focal] >= 2 and counts[BOMBUS] >= 2 and counts[focal] > control_max and counts[BOMBUS] > control_max
            if not strict:
                continue
            copy_aggregate[focal] += 1
            copy_hits.append({
                "orthogroup": og,
                "focal_apis": focal,
                **{f"count_{species}": counts[species] for species in SPECIES},
                "other_apis_max": control_max,
                "candidate_any": is_candidate,
                "candidate_categories": ";".join(sorted(group_categories)),
                "laboriosa_descriptions": " | ".join(laboriosa_descriptions),
            })
    copy_fields = ["orthogroup", "focal_apis"] + [f"count_{species}" for species in SPECIES] + [
        "other_apis_max", "candidate_any", "candidate_categories", "laboriosa_descriptions"
    ]
    write_tsv(output / "strict_copy_number_hits.tsv", copy_hits, copy_fields)
    copy_rows = [{"focal_apis": species, "strict_shared_expansion_orthogroups": copy_aggregate[species]} for species in APIS]
    write_tsv(output / "strict_copy_number_aggregate.tsv", copy_rows, ["focal_apis", "strict_shared_expansion_orthogroups"])

    para_rows = [row for row in gene_summaries if bool(row["is_para"])]
    result = {
        "orthogroups_total": len(orthogroups),
        "complete_single_copy_orthogroups": single_copy_count,
        "callable_amino_acid_sites": total_callable,
        "laboriosa_bombus_strict_sites": sum(int(row["hits_Apis_laboriosa"]) for row in gene_summaries),
        "laboriosa_bombus_strict_site_orthogroups": sum(int(row["hits_Apis_laboriosa"]) > 0 for row in gene_summaries),
        "para_complete_single_copy_orthogroups": len(para_rows),
        "para_laboriosa_bombus_strict_sites": sum(int(row["hits_Apis_laboriosa"]) for row in para_rows),
        "strict_copy_number_hits_laboriosa_bombus": copy_aggregate["Apis_laboriosa"],
        "alignment_source": (
            "OrthoFinder final MultipleSequenceAlignments"
            if use_saved_alignments
            else "pyfamsa recomputation"
        ),
        "interpretation_guardrail": "Descriptive foreground-sharing screen; not a phylogenetic association test or evidence of causality.",
    }
    (output / "analysis_summary.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
