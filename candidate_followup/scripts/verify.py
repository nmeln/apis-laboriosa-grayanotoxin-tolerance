"""Verify exact outputs or frozen dependencies, and inspect TSV shapes."""
import argparse
import csv
from common import ROOT, HERE, OUT
from manage_inputs import digest


def main():
    p = argparse.ArgumentParser(); p.add_argument('--dependencies',action='store_true'); a = p.parse_args()
    manifest = HERE / ('dependencies.sha256' if a.dependencies else 'results.sha256')
    base = ROOT if a.dependencies else HERE
    names = set()
    for line in manifest.read_text().splitlines():
        expected, name = line.split('  ',1)
        assert name not in names; names.add(name)
        assert digest(base / name) == expected, 'Checksum mismatch: ' + name
    if not a.dependencies:
        assert names == {str(p.relative_to(HERE)) for p in OUT.iterdir() if p.is_file()}, 'Result file set changed'
        for path in OUT.glob('*.tsv'):
            with path.open() as f:
                rows = csv.reader(f, delimiter='\t'); width = len(next(rows))
                assert all(len(r) == width for r in rows), 'Malformed table: ' + str(path)
    print(str(len(names)) + (' frozen dependencies' if a.dependencies else ' outputs') + ' match SHA-256 exactly.')


if __name__ == '__main__': main()
