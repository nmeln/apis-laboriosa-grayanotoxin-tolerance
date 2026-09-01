#!/usr/bin/env python3
"""Verify input snapshots, generated results, and TSV structure."""

from __future__ import annotations

import argparse
import csv
import hashlib
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_inputs() -> list[str]:
    errors = []
    with (ROOT / "data_sources.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        path = ROOT / row["relative_path"]
        if not path.is_file():
            errors.append(f"missing input: {row['relative_path']}")
            continue
        observed_size = path.stat().st_size
        observed_hash = sha256(path)
        if observed_size != int(row["size_bytes"]):
            errors.append(
                f"input size mismatch: {row['relative_path']} "
                f"({observed_size} != {row['size_bytes']})"
            )
        if observed_hash != row["sha256"]:
            errors.append(
                f"input SHA-256 mismatch: {row['relative_path']} "
                f"({observed_hash} != {row['sha256']})"
            )
    if not errors:
        print(f"inputs: {len(rows)} exact matches")
    return errors


def read_checksum_file(path: Path) -> list[tuple[str, str]]:
    entries = []
    for line_number, line in enumerate(path.read_text().splitlines(), start=1):
        if not line.strip():
            continue
        try:
            digest, relative_path = line.split("  ", 1)
        except ValueError as error:
            raise RuntimeError(f"Invalid checksum line {line_number}: {line!r}") from error
        entries.append((digest, relative_path))
    return entries


def verify_results() -> list[str]:
    errors = []
    checksum_file = ROOT / "results.sha256"
    if not checksum_file.is_file():
        return ["missing results.sha256"]
    entries = read_checksum_file(checksum_file)
    expected_paths = set()
    for expected_hash, relative_path in entries:
        expected_paths.add(relative_path)
        path = ROOT / relative_path
        if not path.is_file():
            errors.append(f"missing result: {relative_path}")
            continue
        observed_hash = sha256(path)
        if observed_hash != expected_hash:
            errors.append(
                f"result SHA-256 mismatch: {relative_path} "
                f"({observed_hash} != {expected_hash})"
            )
    actual_paths = {
        str(path.relative_to(ROOT))
        for path in (ROOT / "results").iterdir()
        if path.is_file()
    }
    for unexpected in sorted(actual_paths - expected_paths):
        errors.append(f"untracked result file: {unexpected}")
    if not errors:
        print(f"results: {len(entries)} exact matches")
    return errors


def verify_tsv_shapes() -> list[str]:
    errors = []
    paths = sorted((ROOT / "results").glob("*.tsv")) + [
        ROOT / "data_sources.tsv",
        ROOT / "genome_manifest.tsv",
    ]
    checked = 0
    for path in paths:
        if not path.is_file():
            continue
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                errors.append(f"empty TSV: {path.relative_to(ROOT)}")
                continue
            width = len(header)
            for line_number, row in enumerate(reader, start=2):
                if len(row) != width:
                    errors.append(
                        f"TSV width mismatch: {path.relative_to(ROOT)}:{line_number} "
                        f"({len(row)} != {width})"
                    )
        checked += 1
    if not errors:
        print(f"TSV structure: {checked} files consistent")
    return errors


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", action="store_true")
    parser.add_argument("--results", action="store_true")
    args = parser.parse_args()
    if not args.inputs and not args.results:
        args.inputs = args.results = True

    errors = []
    if args.inputs:
        errors.extend(verify_inputs())
    if args.results:
        errors.extend(verify_results())
    errors.extend(verify_tsv_shapes())

    if errors:
        for error in errors:
            print(f"ERROR: {error}", file=sys.stderr)
        raise SystemExit(1)


if __name__ == "__main__":
    main()
