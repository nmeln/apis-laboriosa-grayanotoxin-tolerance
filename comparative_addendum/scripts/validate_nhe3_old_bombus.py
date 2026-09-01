#!/usr/bin/env python3
"""Check the NHE3 foreground-sharing sites against an older Bombus reference."""

from __future__ import annotations

import argparse
import csv
import gzip
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
    for species, sequence in aligned.items():
        position = 0
        values = []
        for aa in sequence:
            if aa == "-":
                values.append(None)
            else:
                position += 1
                values.append(position)
        result[species] = values
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--orthogroups", type=Path, required=True)
    parser.add_argument("--proteomes", type=Path, required=True)
    parser.add_argument("--old-bombus-proteome", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--orthogroup", default="OG0006407")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    with args.orthogroups.open(newline="") as handle:
        rows = [row for row in csv.DictReader(handle, delimiter="\t") if row["Orthogroup"] == args.orthogroup]
    if len(rows) != 1:
        raise RuntimeError(f"Expected one {args.orthogroup} row, found {len(rows)}")
    group = rows[0]
    if any("," in group[species] or not group[species] for species in SPECIES):
        raise RuntimeError(f"{args.orthogroup} is not complete single-copy")

    records = {}
    accessions = {}
    for species in SPECIES:
        accession = group[species].strip()
        proteome = read_fasta(args.proteomes / f"{species}.faa")
        records[species] = proteome[accession][1]
        accessions[species] = accession

    old_records = read_fasta(args.old_bombus_proteome, compressed=True)
    old_candidates = [
        (accession, sequence)
        for accession, (description, sequence) in old_records.items()
        if "sodium/hydrogen exchanger 3 [Bombus terrestris]" in description
    ]
    if len(old_candidates) != 1:
        raise RuntimeError(f"Expected one old Bombus NHE3 model, found {len(old_candidates)}")
    old_accession, old_sequence = old_candidates[0]
    records[OLD_BOMBUS] = old_sequence

    aligned_records = Aligner(threads=1).align([
        Sequence(species.encode(), sequence.encode()) for species, sequence in records.items()
    ])
    aligned = {record.id.decode(): record.sequence.decode() for record in aligned_records}
    maps = position_maps(aligned)

    hits = []
    for column in range(len(next(iter(aligned.values())))):
        residues = {species: aligned[species][column] for species in SPECIES}
        if not all(aa in VALID_AA for aa in residues.values()):
            continue
        controls = [species for species in APIS if species != "Apis_laboriosa"]
        control_states = {residues[species] for species in controls}
        if len(control_states) != 1:
            continue
        control_state = next(iter(control_states))
        if residues["Apis_laboriosa"] != residues[BOMBUS] or residues[BOMBUS] == control_state:
            continue
        hits.append({
            "alignment_column_1based": column + 1,
            "laboriosa_position_1based": maps["Apis_laboriosa"][column],
            "current_bombus_position_1based": maps[BOMBUS][column],
            "old_bombus_position_1based": maps[OLD_BOMBUS][column],
            "laboriosa_state": residues["Apis_laboriosa"],
            "current_bombus_state": residues[BOMBUS],
            "old_bombus_state": aligned[OLD_BOMBUS][column],
            "other_apis_state": control_state,
            "old_matches_current_bombus": aligned[OLD_BOMBUS][column] == residues[BOMBUS],
            "old_bombus_accession": old_accession,
            "current_bombus_accession": accessions[BOMBUS],
        })

    if len(hits) != 4:
        raise RuntimeError(f"Expected four strict NHE3 sites after seven-sequence realignment, found {len(hits)}")
    with (args.output / "nhe3_old_bombus_validation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(hits[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(hits)
    print(f"strict_sites={len(hits)} old_matches={sum(row['old_matches_current_bombus'] for row in hits)}")


if __name__ == "__main__":
    main()
