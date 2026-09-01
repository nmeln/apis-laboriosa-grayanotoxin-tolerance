#!/usr/bin/env python3
"""Fetch and byte-verify the additional comparative-addendum inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import os
import shutil
import tarfile
import tempfile
import time
import urllib.request
import zipfile
from pathlib import Path


ADDENDUM_ROOT = Path(__file__).resolve().parents[1]
MANIFEST = ADDENDUM_ROOT / "input_sources.tsv"
SNAPSHOT_NAME = "comparative-addendum-inputs-v1.tar"
SNAPSHOT_HASH = ADDENDUM_ROOT / "input_snapshot.sha256"
SNAPSHOT_URL = (
    "https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/"
    "releases/download/comparative-addendum-inputs-v1/"
    f"{SNAPSHOT_NAME}"
)
USER_AGENT = "apis-laboriosa-comparative-addendum/1.0"
LIN_SUPPLEMENT_PATH = "inputs/lin2021_supplementary_data.zip"
LIN_SUPPLEMENT_MEMBER = "evab227_supplementary_data.zip"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_manifest() -> list[dict[str, str]]:
    with MANIFEST.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    expected = {"relative_path", "accession_or_dataset", "source_url", "size_bytes", "sha256"}
    if not rows or set(rows[0]) != expected:
        raise RuntimeError(f"Malformed {MANIFEST}")
    if len({row["relative_path"] for row in rows}) != len(rows):
        raise RuntimeError(f"Duplicate paths in {MANIFEST}")
    return rows


def exact_match(path: Path, row: dict[str, str]) -> bool:
    return (
        path.is_file()
        and path.stat().st_size == int(row["size_bytes"])
        and sha256(path) == row["sha256"]
    )


def download(url: str, directory: Path) -> Path:
    directory.mkdir(parents=True, exist_ok=True)
    last_error: Exception | None = None
    for attempt in range(1, 4):
        temporary = tempfile.NamedTemporaryFile(
            prefix="download-", suffix=".part", dir=directory, delete=False
        )
        path = Path(temporary.name)
        try:
            request = urllib.request.Request(url, headers={"User-Agent": USER_AGENT})
            with temporary, urllib.request.urlopen(request, timeout=180) as response:
                shutil.copyfileobj(response, temporary, length=1024 * 1024)
            return path
        except Exception as error:
            last_error = error
            path.unlink(missing_ok=True)
            if attempt < 3:
                time.sleep(attempt * 2)
    raise RuntimeError(f"Download failed after three attempts: {url}") from last_error


def install_verified(temporary: Path, row: dict[str, str]) -> None:
    expected_size = int(row["size_bytes"])
    observed_size = temporary.stat().st_size
    expected_hash = row["sha256"]
    observed_hash = sha256(temporary)
    if observed_size != expected_size or observed_hash != expected_hash:
        temporary.unlink(missing_ok=True)
        raise RuntimeError(
            f"Input mismatch for {row['relative_path']}: size {observed_size} "
            f"vs {expected_size}; sha256 {observed_hash} vs {expected_hash}"
        )
    target = ADDENDUM_ROOT / row["relative_path"]
    target.parent.mkdir(parents=True, exist_ok=True)
    os.replace(temporary, target)
    print(f"verified  {row['relative_path']}")


def expected_snapshot_hash() -> str:
    digest, filename = SNAPSHOT_HASH.read_text().strip().split("  ", 1)
    if filename != SNAPSHOT_NAME or len(digest) != 64:
        raise RuntimeError(f"Malformed {SNAPSHOT_HASH}")
    return digest


def fetch_snapshot(rows: list[dict[str, str]]) -> None:
    archive_path = download(SNAPSHOT_URL, ADDENDUM_ROOT / ".cache")
    try:
        observed = sha256(archive_path)
        expected = expected_snapshot_hash()
        if observed != expected:
            raise RuntimeError(f"Snapshot SHA-256 mismatch: {observed} vs {expected}")
        by_path = {row["relative_path"]: row for row in rows}
        with tarfile.open(archive_path, "r") as archive:
            members = {member.name: member for member in archive.getmembers() if member.isfile()}
            if set(members) != set(by_path):
                raise RuntimeError("Snapshot members do not match input_sources.tsv")
            for relative_path, row in by_path.items():
                source = archive.extractfile(members[relative_path])
                if source is None:
                    raise RuntimeError(f"Cannot read snapshot member: {relative_path}")
                target = ADDENDUM_ROOT / relative_path
                target.parent.mkdir(parents=True, exist_ok=True)
                with tempfile.NamedTemporaryFile(
                    prefix=target.name + ".", suffix=".part", dir=target.parent, delete=False
                ) as temporary:
                    shutil.copyfileobj(source, temporary, length=1024 * 1024)
                    temporary_path = Path(temporary.name)
                install_verified(temporary_path, row)
    finally:
        archive_path.unlink(missing_ok=True)


def fetch_official(rows: list[dict[str, str]], force: bool) -> None:
    for row in rows:
        target = ADDENDUM_ROOT / row["relative_path"]
        if not force and exact_match(target, row):
            print(f"cached    {row['relative_path']}")
            continue
        temporary = download(row["source_url"], target.parent)
        if row["relative_path"] == LIN_SUPPLEMENT_PATH:
            outer_archive = temporary
            try:
                with zipfile.ZipFile(outer_archive) as archive:
                    names = archive.namelist()
                    if names.count(LIN_SUPPLEMENT_MEMBER) != 1:
                        raise RuntimeError(
                            f"Expected one {LIN_SUPPLEMENT_MEMBER} in Europe PMC archive; "
                            f"observed {names}"
                        )
                    with archive.open(LIN_SUPPLEMENT_MEMBER) as source:
                        with tempfile.NamedTemporaryFile(
                            prefix=target.name + ".",
                            suffix=".part",
                            dir=target.parent,
                            delete=False,
                        ) as extracted:
                            shutil.copyfileobj(source, extracted, length=1024 * 1024)
                            temporary = Path(extracted.name)
            finally:
                outer_archive.unlink(missing_ok=True)
        install_verified(temporary, row)


def verify_all(rows: list[dict[str, str]]) -> None:
    bad = [
        row["relative_path"]
        for row in rows
        if not exact_match(ADDENDUM_ROOT / row["relative_path"], row)
    ]
    if bad:
        raise RuntimeError(f"Missing or mismatched inputs: {bad}")
    print(f"all {len(rows)} addendum inputs match")


def main() -> None:
    parser = argparse.ArgumentParser()
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--official-only", action="store_true")
    mode.add_argument("--snapshot-only", action="store_true")
    parser.add_argument("--force", action="store_true")
    args = parser.parse_args()
    rows = load_manifest()

    if not args.force and all(
        exact_match(ADDENDUM_ROOT / row["relative_path"], row) for row in rows
    ):
        verify_all(rows)
        return
    if args.official_only:
        fetch_official(rows, args.force)
    else:
        try:
            fetch_snapshot(rows)
        except Exception:
            if args.snapshot_only:
                raise
            fetch_official(rows, args.force)
    verify_all(rows)


if __name__ == "__main__":
    main()
