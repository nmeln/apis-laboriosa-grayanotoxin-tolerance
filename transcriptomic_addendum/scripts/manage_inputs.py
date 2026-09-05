#!/usr/bin/env python3
"""Fetch, verify, or archive the exact ten added source files."""
import argparse
import csv
import hashlib
import io
import subprocess
import sys
import tarfile
import urllib.request
import zipfile
from pathlib import Path

HERE = Path(__file__).resolve().parents[1]
NAME = 'transcriptomic-addendum-inputs-v1.tar'
URL = 'https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/download/transcriptomic-inputs-v1/' + NAME

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''):
            h.update(block)
    return h.hexdigest()

def rows():
    with (HERE/'input_sources.tsv').open() as f:
        return list(csv.DictReader(f, delimiter='\t'))

def verify(path, row):
    assert path.is_file(), f'Missing {path}'
    assert path.stat().st_size == int(row['size_bytes']), f'Size mismatch: {path}'
    assert digest(path) == row['sha256'], f'SHA-256 mismatch: {path}; preserve and investigate'

def download(url, path):
    request = urllib.request.Request(url, headers={'User-Agent':'bee-transcript-audit/1.0'})
    with urllib.request.urlopen(request, timeout=120) as response, path.open('wb') as out:
        for block in iter(lambda: response.read(1024*1024), b''):
            out.write(block)

def main():
    p = argparse.ArgumentParser()
    p.add_argument('action', choices=['fetch','verify','snapshot'])
    p.add_argument('--official', action='store_true', help='Skip release and restore from original repositories')
    args = p.parse_args()
    manifest = rows()
    assert len(manifest) == 10
    (HERE/'inputs').mkdir(exist_ok=True)
    (HERE/'work').mkdir(exist_ok=True)
    if args.action == 'fetch':
        missing = [r for r in manifest if not (HERE/r['relative_path']).exists()]
        # Existing mismatches are errors, never silently replaced.
        for row in manifest:
            if (HERE/row['relative_path']).exists(): verify(HERE/row['relative_path'], row)
        if missing and not args.official:
            archive = HERE/'work'/NAME
            try:
                download(URL, archive)
            except (OSError, urllib.error.URLError) as exc:
                print(f'Release unavailable ({exc}); trying official sources.', flush=True)
            else:
                expected = (HERE/'input_snapshot.sha256').read_text().split()[0]
                assert digest(archive) == expected, 'Snapshot hash mismatch; refusing fallback'
                with tarfile.open(archive) as tar:
                    assert {m.name for m in tar.getmembers()} == {r['relative_path'] for r in manifest}
                    for row in missing:
                        member = tar.getmember(row['relative_path'])
                        assert member.isfile()
                        path = HERE/row['relative_path']
                        path.write_bytes(tar.extractfile(member).read())
                        verify(path, row)
        for row in manifest:
            path = HERE/row['relative_path']
            if path.exists(): continue
            method = row['retrieval_method']
            if method.startswith('first_'):
                run = path.name.split('.')[0]
                subprocess.run([sys.executable, str(HERE/'scripts/fetch_raw_subset.py'), run], check=True)
            else:
                temporary = path.with_suffix(path.suffix+'.part')
                download(row['source_url'], temporary)
                if method.startswith('zip_member:'):
                    member = method.split(':', 1)[1]
                    with zipfile.ZipFile(temporary) as z:
                        matches = [n for n in z.namelist() if Path(n).name == member]
                        assert len(matches) == 1
                        contents = z.read(matches[0])
                    temporary.write_bytes(contents)
                verify(temporary, row)
                temporary.replace(path)
            verify(path, row)
    for row in manifest: verify(HERE/row['relative_path'], row)
    print(f'All {len(manifest)} added inputs match size and SHA-256.', flush=True)
    if args.action == 'snapshot':
        path = HERE/NAME
        with tarfile.open(path, 'w', format=tarfile.USTAR_FORMAT) as tar:
            for row in manifest:
                source = HERE/row['relative_path']
                info = tarfile.TarInfo(row['relative_path'])
                info.size = source.stat().st_size
                info.mode = 0o644
                info.uid = info.gid = info.mtime = 0
                with source.open('rb') as f: tar.addfile(info, f)
        line = digest(path)+'  '+NAME+'\n'
        checksum = HERE/'input_snapshot.sha256'
        if checksum.exists():
            assert checksum.read_text() == line, 'Archive changed; investigate before versioning a new snapshot'
        else:
            checksum.write_text(line)
        print(line, end='')

if __name__ == '__main__': main()
