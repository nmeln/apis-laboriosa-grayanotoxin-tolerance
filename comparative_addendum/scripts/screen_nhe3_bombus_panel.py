#!/usr/bin/env python3
"""Evaluate the four NHE3 sharing sites in an expanded Bombus panel."""

from __future__ import annotations

import argparse
import csv
import gzip
import re
from pathlib import Path

from pyfamsa import Aligner, Sequence


APIS = [
    "Apis_laboriosa",
    "Apis_dorsata",
    "Apis_mellifera",
    "Apis_cerana",
    "Apis_florea",
]
CURRENT_BOMBUS = "Bombus_terrestris_current"
OLD_BOMBUS = "Bombus_terrestris_old"
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path, compressed: bool = False) -> dict[str, tuple[str, str]]:
    records: dict[str, tuple[str, str]] = {}
    accession: str | None = None
    description = ""
    chunks: list[str] = []
    opener = gzip.open if compressed else open
    with opener(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if accession is not None:
                    records[accession] = (description, "".join(chunks).replace("*", ""))
                description = line[1:].strip()
                accession = description.split(maxsplit=1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if accession is not None:
        records[accession] = (description, "".join(chunks).replace("*", ""))
    return records


def position_maps(aligned: dict[str, str]) -> dict[str, list[int | None]]:
    result = {}
    for key, sequence in aligned.items():
        position = 0
        values = []
        for aa in sequence:
            if aa == "-":
                values.append(None)
            else:
                position += 1
                values.append(position)
        result[key] = values
    return result


def taxon_from_description(description: str) -> str:
    match = re.search(r"\[([^]]+)\]$", description)
    if not match:
        raise RuntimeError(f"Cannot parse taxon from {description}")
    return match.group(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--proteomes", type=Path, required=True)
    parser.add_argument("--old-bombus-proteome", type=Path, required=True)
    parser.add_argument("--bombus-panel", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--orthogroup", default="OG0006407")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.orthogroups.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["Orthogroup"] == args.orthogroup]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one {args.orthogroup} row, found {len(rows)}")
    group = rows[0]

    sequences = {}
    accession_by_key = {}
    taxon_by_key = {}
    for species in APIS + ["Bombus_terrestris"]:
        accession = group[species].strip()
        if not accession or "," in accession:
            raise RuntimeError(f"{args.orthogroup} is not single-copy for {species}")
        proteome = read_fasta(args.proteomes / f"{species}.faa")
        key = CURRENT_BOMBUS if species == "Bombus_terrestris" else species
        sequences[key] = proteome[accession][1]
        accession_by_key[key] = accession
        taxon_by_key[key] = species.replace("_", " ")

    old_records = read_fasta(args.old_bombus_proteome, compressed=True)
    old_candidates = [
        (accession, description, sequence)
        for accession, (description, sequence) in old_records.items()
        if "sodium/hydrogen exchanger 3 [Bombus terrestris]" in description
    ]
    if len(old_candidates) != 1:
        raise RuntimeError(f"Expected one old Bombus NHE3 model, found {len(old_candidates)}")
    old_accession, old_description, old_sequence = old_candidates[0]
    sequences[OLD_BOMBUS] = old_sequence
    accession_by_key[OLD_BOMBUS] = old_accession
    taxon_by_key[OLD_BOMBUS] = taxon_from_description(old_description)

    panel_records = read_fasta(args.bombus_panel)
    for accession, (description, sequence) in panel_records.items():
        key = f"panel_{accession}"
        sequences[key] = sequence
        accession_by_key[key] = accession
        taxon_by_key[key] = taxon_from_description(description)

    aligned_records = Aligner(threads=1).align([
        Sequence(key.encode(), sequence.encode()) for key, sequence in sequences.items()
    ])
    aligned = {record.id.decode(): record.sequence.decode() for record in aligned_records}
    maps = position_maps(aligned)

    strict_columns = []
    controls = [species for species in APIS if species != "Apis_laboriosa"]
    for column in range(len(next(iter(aligned.values())))):
        residues = {species: aligned[species][column] for species in APIS}
        residues[CURRENT_BOMBUS] = aligned[CURRENT_BOMBUS][column]
        if not all(aa in VALID_AA for aa in residues.values()):
            continue
        control_states = {residues[species] for species in controls}
        if len(control_states) != 1:
            continue
        control_state = next(iter(control_states))
        if (
            residues["Apis_laboriosa"] == residues[CURRENT_BOMBUS]
            and residues["Apis_laboriosa"] != control_state
        ):
            strict_columns.append((column, control_state))
    if len(strict_columns) != 4:
        raise RuntimeError(f"Expected four strict NHE3 sites, found {len(strict_columns)}")

    panel_keys = [OLD_BOMBUS] + sorted(key for key in sequences if key.startswith("panel_"))
    detail = []
    site_summary = []
    for site_number, (column, control_state) in enumerate(strict_columns, start=1):
        focal_state = aligned["Apis_laboriosa"][column]
        matches = 0
        valid = 0
        for key in panel_keys:
            state = aligned[key][column]
            if state in VALID_AA:
                valid += 1
                matches += state == focal_state
            detail.append({
                "site_number": site_number,
                "laboriosa_position_1based": maps["Apis_laboriosa"][column],
                "current_bombus_position_1based": maps[CURRENT_BOMBUS][column],
                "laboriosa_and_current_bombus_state": focal_state,
                "other_apis_state": control_state,
                "panel_taxon": taxon_by_key[key],
                "panel_accession": accession_by_key[key],
                "panel_position_1based": maps[key][column],
                "panel_state": state,
                "matches_laboriosa_and_current_bombus": state == focal_state,
            })
        site_summary.append({
            "site_number": site_number,
            "laboriosa_position_1based": maps["Apis_laboriosa"][column],
            "current_bombus_position_1based": maps[CURRENT_BOMBUS][column],
            "laboriosa_and_current_bombus_state": focal_state,
            "other_apis_state": control_state,
            "additional_bombus_sequences_with_valid_state": valid,
            "additional_bombus_sequences_matching_focal_state": matches,
        })

    with (args.output / "nhe3_bombus_panel_sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(detail[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(detail)
    with (args.output / "nhe3_bombus_panel_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(site_summary[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(site_summary)
    print(f"sites={len(site_summary)} panel_sequences={len(panel_keys)}")


if __name__ == "__main__":
    main()
