#!/usr/bin/env python3
"""Verify addendum inputs, committed outputs, scripts, and generated work files."""

from __future__ import annotations

import argparse
import csv
import hashlib
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
PROJECT = ROOT.parent


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def verify_sha_manifest(path: Path) -> set[str]:
    errors = []
    recorded: set[str] = set()
    for line in path.read_text().splitlines():
        expected, relative = line.split("  ", 1)
        if relative in recorded:
            errors.append(f"duplicate manifest path {relative}")
            continue
        recorded.add(relative)
        target = PROJECT / relative
        if not target.is_file():
            errors.append(f"missing {relative}")
        elif sha256(target) != expected:
            errors.append(f"SHA-256 mismatch {relative}")
    if errors:
        raise RuntimeError("; ".join(errors))
    print(f"verified {path.name}")
    return recorded


def require_exact_files(recorded: set[str], files: list[Path], label: str) -> None:
    actual = {str(path.relative_to(PROJECT)) for path in files if path.is_file()}
    missing = sorted(recorded - actual)
    extra = sorted(actual - recorded)
    if missing or extra:
        parts = []
        if missing:
            parts.append(f"unlisted missing {label} paths: {missing}")
        if extra:
            parts.append(f"unhashed extra {label} paths: {extra}")
        raise RuntimeError("; ".join(parts))


def verify_result_tables() -> None:
    errors = []
    for path in sorted((ROOT / "results").rglob("*.tsv")):
        if path.name == "Statistics_Overall.tsv":
            continue
        with path.open(newline="") as handle:
            reader = csv.reader(handle, delimiter="\t")
            try:
                header = next(reader)
            except StopIteration:
                errors.append(f"empty TSV {path.relative_to(PROJECT)}")
                continue
            for line_number, row in enumerate(reader, 2):
                if len(row) != len(header):
                    errors.append(
                        f"TSV width {path.relative_to(PROJECT)}:{line_number} "
                        f"expected {len(header)}, observed {len(row)}"
                    )
                    break
    if errors:
        raise RuntimeError("; ".join(errors))
    print("verified result TSV structure")


def verify_inputs() -> None:
    with (ROOT / "input_sources.tsv").open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    errors = []
    for row in rows:
        path = ROOT / row["relative_path"]
        if not path.is_file():
            errors.append(f"missing {row['relative_path']}")
        elif path.stat().st_size != int(row["size_bytes"]):
            errors.append(f"size mismatch {row['relative_path']}")
        elif sha256(path) != row["sha256"]:
            errors.append(f"SHA-256 mismatch {row['relative_path']}")
    if errors:
        raise RuntimeError("; ".join(errors))
    print(f"verified {len(rows)} addendum inputs")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--inputs", action="store_true")
    parser.add_argument("--results", action="store_true")
    parser.add_argument("--scripts", action="store_true")
    parser.add_argument("--work", action="store_true")
    args = parser.parse_args()
    if not any((args.inputs, args.results, args.scripts, args.work)):
        args.results = args.scripts = True
    if args.inputs:
        verify_inputs()
    if args.results:
        recorded = verify_sha_manifest(ROOT / "results.sha256")
        require_exact_files(
            recorded,
            [path for path in (ROOT / "results").rglob("*") if path.is_file()],
            "result",
        )
        verify_result_tables()
    if args.scripts:
        recorded = verify_sha_manifest(ROOT / "scripts.sha256")
        require_exact_files(
            recorded,
            sorted((ROOT / "scripts").glob("*.py")) + [ROOT / "run_analysis.sh"],
            "script",
        )
    if args.work:
        verify_sha_manifest(ROOT / "primary_proteomes.sha256")
        verify_sha_manifest(ROOT / "orthofinder_key_outputs.sha256")


if __name__ == "__main__":
    main()
