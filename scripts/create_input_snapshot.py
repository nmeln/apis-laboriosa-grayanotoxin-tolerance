#!/usr/bin/env python3
"""Create a deterministic tar archive of the pinned public inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import tarfile
import tempfile
from pathlib import Path


ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "data_sources.tsv"
DEFAULT_OUTPUT = ROOT / "apis-laboriosa-inputs-v1.tar"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def verify_inputs(rows: list[dict[str, str]]) -> None:
    errors = []
    for row in rows:
        path = ROOT / row["relative_path"]
        if not path.is_file():
            errors.append(f"missing: {row['relative_path']}")
            continue
        if path.stat().st_size != int(row["size_bytes"]):
            errors.append(f"size mismatch: {row['relative_path']}")
        if sha256(path) != row["sha256"]:
            errors.append(f"SHA-256 mismatch: {row['relative_path']}")
    if errors:
        raise RuntimeError("; ".join(errors))


def create_archive(rows: list[dict[str, str]], output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = tempfile.NamedTemporaryFile(
        prefix=f"{output.name}.",
        suffix=".part",
        dir=output.parent,
        delete=False,
    )
    temporary_path = Path(temporary.name)
    try:
        with temporary:
            with tarfile.open(
                fileobj=temporary,
                mode="w",
                format=tarfile.USTAR_FORMAT,
            ) as archive:
                for row in rows:
                    relative_path = row["relative_path"]
                    path = ROOT / relative_path
                    info = tarfile.TarInfo(name=relative_path)
                    info.size = path.stat().st_size
                    info.mode = 0o644
                    info.mtime = 0
                    info.uid = 0
                    info.gid = 0
                    info.uname = ""
                    info.gname = ""
                    with path.open("rb") as source:
                        archive.addfile(info, source)
        os.replace(temporary_path, output)
    except Exception:
        temporary_path.unlink(missing_ok=True)
        raise


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, default=DEFAULT_OUTPUT)
    args = parser.parse_args()
    output = args.output if args.output.is_absolute() else ROOT / args.output

    rows = load_manifest()
    verify_inputs(rows)
    create_archive(rows, output)
    print(f"snapshot files: {len(rows)}")
    print(f"snapshot bytes: {output.stat().st_size}")
    print(f"snapshot sha256: {sha256(output)}")
    print(f"snapshot path: {output}")


if __name__ == "__main__":
    main()
