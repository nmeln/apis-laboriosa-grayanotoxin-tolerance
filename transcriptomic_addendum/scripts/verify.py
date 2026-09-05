#!/usr/bin/env python3
"""Verify pinned source files, frozen orthology dependencies, and outputs."""
import argparse
import csv
import subprocess
import sys
from pathlib import Path
from manage_inputs import HERE, digest, rows, verify

def checksum_manifest(path, root):
    names=set()
    for line in path.read_text().splitlines():
        expected,name=line.split('  ',1)
        assert name not in names, 'Duplicate checksum path'
        names.add(name)
        assert digest(root/name)==expected, f'Hash mismatch: {name}'
    return names

def main():
    p=argparse.ArgumentParser();p.add_argument('--results',action='store_true');a=p.parse_args()
    if a.results:
        names=checksum_manifest(HERE/'results.sha256',HERE)
        actual={str(f.relative_to(HERE)) for f in (HERE/'results').iterdir() if f.is_file()}
        assert names==actual, f'Output file set changed: {names ^ actual}'
        for path in (HERE/'results').glob('*.tsv'):
            with path.open() as f:
                r=csv.reader(f,delimiter='\t'); width=len(next(r))
                assert all(len(row)==width for row in r), f'TSV shape: {path}'
        print(f'{len(names)} addendum results match byte for byte; TSV shapes valid.')
    else:
        subprocess.run([sys.executable,str(HERE.parent/'scripts/verify_project.py'),'--inputs'],check=True)
        for row in rows():verify(HERE/row['relative_path'],row)
        dependencies=checksum_manifest(HERE/'dependencies.sha256',HERE.parent)
        print(f'{len(rows())} added sources and {len(dependencies)} frozen orthology files verified.')

if __name__=='__main__':main()
