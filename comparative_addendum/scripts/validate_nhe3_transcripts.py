#!/usr/bin/env python3
"""Validate NHE3 protein states against the assembled worker transcripts."""

from __future__ import annotations

import argparse
import csv
import gzip
import subprocess
import tempfile
from pathlib import Path


TARGET_POSITIONS = [44, 159, 232, 353]
UNIGENES = {
    "Apis_laboriosa": "references/transcriptome_2019/GSM3757258_AL.unigene.fasta.gz",
    "Apis_dorsata": "references/transcriptome_2019/GSM3757259_AD.unigene.fasta.gz",
}


def read_plain_fasta(path: Path) -> dict[str, str]:
    records = {}
    key = None
    chunks = []
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                if key is not None:
                    records[key] = "".join(chunks)
                key = line[1:].split(maxsplit=1)[0]
                chunks = []
            else:
                chunks.append(line.strip())
    if key is not None:
        records[key] = "".join(chunks)
    return records


def get_gzip_fasta_record(path: Path, wanted: str) -> str:
    key = None
    chunks = []
    with gzip.open(path, "rt") as handle:
        for line in handle:
            if line.startswith(">"):
                if key == wanted:
                    return "".join(chunks)
                key = line[1:].split(maxsplit=1)[0]
                chunks = []
            elif key == wanted:
                chunks.append(line.strip())
    if key == wanted:
        return "".join(chunks)
    raise RuntimeError(f"Unigene {wanted} not found in {path}")


def blast_alignment(tblastn: Path, protein: str, nucleotide: str) -> dict[str, object]:
    with tempfile.TemporaryDirectory(prefix="nhe3_transcript_") as tmp:
        tmpdir = Path(tmp)
        query = tmpdir / "query.faa"
        subject = tmpdir / "subject.fna"
        query.write_text(">query\n" + protein + "\n")
        subject.write_text(">subject\n" + nucleotide + "\n")
        command = [
            str(tblastn),
            "-query", str(query),
            "-subject", str(subject),
            "-seg", "no",
            "-evalue", "1e-20",
            "-max_hsps", "1",
            "-outfmt", "6 pident length qstart qend sstart send sframe qseq sseq",
        ]
        completed = subprocess.run(command, check=True, capture_output=True, text=True)
    lines = [line for line in completed.stdout.splitlines() if line]
    if len(lines) != 1:
        raise RuntimeError(f"Expected one TBLASTN HSP, found {len(lines)}")
    fields = lines[0].split("\t")
    if len(fields) != 9:
        raise RuntimeError(f"Unexpected TBLASTN output: {lines[0]}")
    return {
        "percent_identity": float(fields[0]),
        "alignment_length_aa": int(fields[1]),
        "query_start": int(fields[2]),
        "query_end": int(fields[3]),
        "subject_start": int(fields[4]),
        "subject_end": int(fields[5]),
        "subject_frame": int(fields[6]),
        "query_alignment": fields[7],
        "subject_alignment": fields[8],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--project", type=Path, required=True)
    parser.add_argument("--proteomes", type=Path, required=True)
    parser.add_argument("--members", type=Path, required=True)
    parser.add_argument("--expression-mappings", type=Path, required=True)
    parser.add_argument("--tblastn", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    nhe3_rows = {}
    with args.expression_mappings.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            if row["category"] == "NHE3":
                nhe3_rows[row["species"]] = row
    if set(nhe3_rows) != set(UNIGENES):
        raise RuntimeError(f"Unexpected NHE3 expression species: {sorted(nhe3_rows)}")

    primary_by_gene = {}
    with args.members.open(newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            primary_by_gene[(row["species"], f"gene-{row['gene']}")] = row["accession"]

    detail = []
    summaries = []
    for species, mapping in nhe3_rows.items():
        accession = primary_by_gene[(species, mapping["gene_parent"])]
        protein = read_plain_fasta(args.proteomes / f"{species}.faa")[accession]
        transcript = get_gzip_fasta_record(
            args.project / UNIGENES[species], mapping["best_unigene"]
        )
        alignment = blast_alignment(args.tblastn, protein, transcript)
        residues = {}
        query_position = int(alignment["query_start"]) - 1
        for query_aa, transcript_aa in zip(
            str(alignment["query_alignment"]), str(alignment["subject_alignment"])
        ):
            if query_aa != "-":
                query_position += 1
                residues[query_position] = (query_aa, transcript_aa)
        for position in TARGET_POSITIONS:
            if position not in residues:
                raise RuntimeError(f"{species} transcript does not cover NHE3 position {position}")
            reference_state, transcript_state = residues[position]
            detail.append({
                "species": species,
                "gene_parent": mapping["gene_parent"],
                "protein_accession": accession,
                "unigene": mapping["best_unigene"],
                "protein_position_1based": position,
                "reference_protein_state": reference_state,
                "transcript_translation_state": transcript_state,
                "state_matches": reference_state == transcript_state,
            })
        summaries.append({
            "species": species,
            "protein_accession": accession,
            "unigene": mapping["best_unigene"],
            "protein_length_aa": len(protein),
            "transcript_length_nt": len(transcript),
            "aligned_query_start": alignment["query_start"],
            "aligned_query_end": alignment["query_end"],
            "alignment_length_aa": alignment["alignment_length_aa"],
            "percent_identity": alignment["percent_identity"],
            "subject_frame": alignment["subject_frame"],
            "target_positions_covered": len(TARGET_POSITIONS),
            "target_states_matching_reference": sum(
                row["state_matches"] for row in detail if row["species"] == species
            ),
            "design_warning": "processed pooled whole-worker transcript assembly; no replicate or allele-frequency estimate",
        })

    with (args.output / "nhe3_transcript_residue_validation.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=list(detail[0]), delimiter="\t", lineterminator="\n"
        )
        writer.writeheader()
        writer.writerows(detail)
    with (args.output / "nhe3_transcript_alignment_summary.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=list(summaries[0]),
            delimiter="\t",
            lineterminator="\n",
        )
        writer.writeheader()
        writer.writerows(summaries)
    print(f"validated_states={sum(row['state_matches'] for row in detail)}/{len(detail)}")


if __name__ == "__main__":
    main()
