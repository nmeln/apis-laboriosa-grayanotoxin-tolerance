#!/usr/bin/env python3
"""Evaluate the NHE2/3-labelled orthogroup sites in other sequenced bees."""

from __future__ import annotations

import argparse
import csv
import re
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
CONTROLS = [species for species in SPECIES[:-1] if species != "Apis_laboriosa"]
VALID_AA = set("ACDEFGHIKLMNPQRSTVWY")


def read_fasta(path: Path) -> dict[str, tuple[str, str]]:
    records = {}
    accession = None
    description = ""
    chunks = []
    with path.open() as handle:
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


def taxon(description: str) -> str:
    match = re.search(r"\[([^]]+)\]$", description)
    if not match:
        raise RuntimeError(f"Cannot parse taxon from {description}")
    return match.group(1)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--proteomes", type=Path, required=True)
    parser.add_argument("--panel", type=Path, action="append", required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--orthogroup", default="OG0006407")
    parser.add_argument("--prefix", default="nhe3_external_bee")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.orthogroups.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["Orthogroup"] == args.orthogroup]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one {args.orthogroup} row, found {len(rows)}")
    group = rows[0]

    sequences = {}
    accessions = {}
    taxa = {}
    for species in SPECIES:
        accession = group[species].strip()
        if not accession or "," in accession:
            raise RuntimeError(f"{args.orthogroup} is not single-copy for {species}")
        records = read_fasta(args.proteomes / f"{species}.faa")
        sequences[species] = records[accession][1]
        accessions[species] = accession
        taxa[species] = species.replace("_", " ")

    excluded = []
    for panel_path in args.panel:
        for accession, (description, sequence) in read_fasta(panel_path).items():
            panel_taxon = taxon(description)
            if not 900 <= len(sequence) <= 1200:
                excluded.append((accession, panel_taxon, len(sequence)))
                continue
            key = f"external_{accession}"
            if key in sequences:
                raise RuntimeError(f"Duplicate panel accession {accession}")
            sequences[key] = sequence
            accessions[key] = accession
            taxa[key] = panel_taxon

    external_keys = sorted(key for key in sequences if key.startswith("external_"))
    if len(external_keys) < 3:
        raise RuntimeError("Too few length-compatible external bee sequences")
    if len({taxa[key] for key in external_keys}) != len(external_keys):
        raise RuntimeError("Expected one length-compatible sequence per external bee taxon")

    aligned_records = Aligner(threads=1).align([
        Sequence(key.encode(), sequence.encode()) for key, sequence in sequences.items()
    ])
    aligned = {record.id.decode(): record.sequence.decode() for record in aligned_records}
    maps = position_maps(aligned)

    with (args.output / f"{args.prefix}_alignment.faa").open("w") as handle:
        for key in sorted(aligned):
            handle.write(f">{key}\n")
            sequence = aligned[key]
            for start in range(0, len(sequence), 80):
                handle.write(sequence[start : start + 80] + "\n")
    metadata = [
        {
            "alignment_id": key,
            "taxon": taxa[key],
            "accession": accessions[key],
            "unaligned_length_aa": len(sequences[key]),
        }
        for key in sorted(sequences)
    ]
    with (args.output / f"{args.prefix}_sequence_metadata.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(metadata[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(metadata)

    strict_columns = []
    for column in range(len(next(iter(aligned.values())))):
        states = {species: aligned[species][column] for species in SPECIES}
        if not all(state in VALID_AA for state in states.values()):
            continue
        control_states = {states[species] for species in CONTROLS}
        if len(control_states) != 1:
            continue
        control_state = next(iter(control_states))
        if (
            states["Apis_laboriosa"] == states["Bombus_terrestris"]
            and states["Apis_laboriosa"] != control_state
        ):
            strict_columns.append((column, control_state))
    if len(strict_columns) != 4:
        raise RuntimeError(f"Expected four strict sites, found {len(strict_columns)}")

    detail = []
    summaries = []
    for site_number, (column, control_state) in enumerate(strict_columns, start=1):
        focal_state = aligned["Apis_laboriosa"][column]
        for key in external_keys:
            state = aligned[key][column]
            detail.append({
                "site_number": site_number,
                "laboriosa_position_1based": maps["Apis_laboriosa"][column],
                "laboriosa_and_bombus_terrestris_state": focal_state,
                "other_four_apis_state": control_state,
                "external_taxon": taxa[key],
                "external_accession": accessions[key],
                "external_position_1based": maps[key][column],
                "external_state": state,
                "matches_focal_state": state == focal_state,
                "matches_other_four_apis_state": state == control_state,
            })
        site_rows = [row for row in detail if row["site_number"] == site_number]
        summaries.append({
            "site_number": site_number,
            "laboriosa_position_1based": maps["Apis_laboriosa"][column],
            "laboriosa_and_bombus_terrestris_state": focal_state,
            "other_four_apis_state": control_state,
            "external_bee_sequences": len(site_rows),
            "external_bees_matching_focal_state": sum(row["matches_focal_state"] for row in site_rows),
            "external_bees_matching_other_apis_state": sum(row["matches_other_four_apis_state"] for row in site_rows),
        })

    with (args.output / f"{args.prefix}_sites.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(detail[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(detail)
    with (args.output / f"{args.prefix}_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(summaries[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(summaries)
    with (args.output / f"{args.prefix}_excluded.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["accession", "taxon", "length_aa", "reason"])
        for accession, panel_taxon, length in excluded:
            writer.writerow([accession, panel_taxon, length, "outside 900-1200 aa orthogroup length window"])
    print(f"sites={len(summaries)} external_taxa={len(external_keys)} excluded={len(excluded)}")


if __name__ == "__main__":
    main()
