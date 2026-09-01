#!/usr/bin/env python3
"""Test exact Bombus-plus-one-Apis amino-acid sharing across Para.

OrthoFinder groups Para with the homologous 60E channel in this six-bee
proteome set.  That makes Para ineligible for a complete one-copy orthogroup
screen.  This focused analysis identifies the primary Para model by product
annotation, aligns exactly one Para protein per species, and applies the same
strict state-sharing rule used by screen_strict_convergence.py.

The result is descriptive.  It does not identify ancestral states, test a
phenotype association, or establish toxin tolerance in Apis laboriosa.
"""

from __future__ import annotations

import argparse
import csv
import json
from collections import Counter
from pathlib import Path

from pyfamsa import Aligner, Sequence


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


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--proteomes", type=Path, required=True)
    parser.add_argument("--members", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    para_rows: dict[str, dict[str, str]] = {}
    with args.members.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            description = row["description"].lower()
            is_para = (
                "sodium channel protein para" in description
                or "sodium channel protein paralytic" in description
                or "sodium voltage-gated channel paralytic" in description
            )
            if is_para:
                species = row["species"]
                if species in para_rows:
                    raise RuntimeError(f"Multiple primary Para annotations for {species}")
                para_rows[species] = row
    if set(para_rows) != set(SPECIES):
        raise RuntimeError(f"Para species mismatch: {sorted(para_rows)}")

    sequences = {}
    for species in SPECIES:
        records = read_fasta(args.proteomes / f"{species}.faa")
        accession = para_rows[species]["accession"]
        sequences[species] = records[accession]

    aligned_records = Aligner(threads=1, refine=False).align([
        Sequence(species.encode(), sequences[species].encode())
        for species in SPECIES
    ])
    aligned = {record.id.decode(): record.sequence.decode() for record in aligned_records}
    if set(aligned) != set(SPECIES):
        raise RuntimeError("Alignment does not contain exactly the six expected species")
    if len({len(sequence) for sequence in aligned.values()}) != 1:
        raise RuntimeError("Alignment sequences have unequal lengths")
    maps = position_maps(aligned)

    hits: list[dict[str, object]] = []
    hit_counts = Counter()
    callable_sites = 0
    for column in range(len(next(iter(aligned.values())))):
        residues = {species: aligned[species][column] for species in SPECIES}
        if not all(aa in VALID_AA for aa in residues.values()):
            continue
        callable_sites += 1
        for focal in APIS:
            controls = [species for species in APIS if species != focal]
            control_states = {residues[species] for species in controls}
            if len(control_states) != 1:
                continue
            control_state = next(iter(control_states))
            if residues[focal] == residues[BOMBUS] and residues[focal] != control_state:
                hit_counts[focal] += 1
                hits.append({
                    "alignment_column_1based": column + 1,
                    "focal_apis": focal,
                    "focal_position_1based": maps[focal][column],
                    "laboriosa_position_1based": maps["Apis_laboriosa"][column],
                    "bombus_position_1based": maps[BOMBUS][column],
                    "focal_and_bombus_state": residues[focal],
                    "other_apis_state": control_state,
                    **{f"state_{species}": residues[species] for species in SPECIES},
                })

    hit_fields = [
        "alignment_column_1based",
        "focal_apis",
        "focal_position_1based",
        "laboriosa_position_1based",
        "bombus_position_1based",
        "focal_and_bombus_state",
        "other_apis_state",
    ] + [f"state_{species}" for species in SPECIES]
    with (args.output / "para_strict_shared_sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=hit_fields, delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(hits)

    aggregate = [
        {
            "focal_apis": focal,
            "callable_sites": callable_sites,
            "strict_sites_shared_with_bombus": hit_counts[focal],
            "para_accession": para_rows[focal]["accession"],
            "para_length_aa": len(sequences[focal]),
        }
        for focal in APIS
    ]
    with (args.output / "para_strict_sharing_aggregate.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(aggregate[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(aggregate)

    summary = {
        "alignment_columns": len(next(iter(aligned.values()))),
        "callable_amino_acid_sites": callable_sites,
        "laboriosa_bombus_strict_shared_sites": hit_counts["Apis_laboriosa"],
        "interpretation_guardrail": (
            "Exact state-sharing screen only; no ancestral-state, phenotype-association, "
            "or causal inference."
        ),
    }
    (args.output / "para_strict_sharing_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n"
    )
    print(json.dumps(summary, indent=2))


if __name__ == "__main__":
    main()
