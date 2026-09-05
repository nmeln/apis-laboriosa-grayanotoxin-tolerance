"""Restore, verify or archive all added source files without changing bytes."""
import argparse
import hashlib
import tarfile
import urllib.request
import urllib.error
from pathlib import Path
from common import HERE, INPUT, WORK, table

NAME = 'candidate-followup-inputs-v1.tar'
URL = 'https://github.com/nmeln/apis-laboriosa-grayanotoxin-tolerance/releases/download/candidate-inputs-v1/' + NAME


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as f:
        for block in iter(lambda: f.read(1024*1024), b''): h.update(block)
    return h.hexdigest()


def verify(path, row):
    assert path.is_file(), 'Missing input: ' + str(path)
    assert path.stat().st_size == int(row['size_bytes']), 'Input size changed: ' + str(path)
    assert digest(path) == row['sha256'], 'Input SHA-256 changed; preserve and investigate: ' + str(path)


def download(url, destination):
    req = urllib.request.Request(url, headers={'User-Agent':'bee-candidate-followup/1.0'})
    with urllib.request.urlopen(req, timeout=120) as response, destination.open('wb') as out:
        for block in iter(lambda: response.read(1024*1024), b''): out.write(block)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('action', choices=['fetch','verify','snapshot'])
    parser.add_argument('--official', action='store_true', help='Use original source URLs; live query responses may have changed since the snapshot')
    args = parser.parse_args()
    rows = table(HERE / 'input_sources.tsv')
    assert len(rows) == 28
    assert len({r['relative_path'] for r in rows}) == len(rows)
    assert all(Path(r['relative_path']).parts[0] == 'inputs' and '..' not in Path(r['relative_path']).parts for r in rows)
    INPUT.mkdir(exist_ok=True); WORK.mkdir(exist_ok=True)
    for r in rows:
        if (HERE / r['relative_path']).exists(): verify(HERE / r['relative_path'], r)
    if args.action == 'fetch':
        missing = [r for r in rows if not (HERE / r['relative_path']).exists()]
        if missing and not args.official:
            archive = WORK / NAME
            download(URL, archive)
            expected = (HERE / 'input_snapshot.sha256').read_text().split()[0]
            assert digest(archive) == expected, 'Input archive mismatch; refusing extraction'
            with tarfile.open(archive) as tar:
                assert {m.name for m in tar.getmembers()} == {r['relative_path'] for r in rows}
                for r in missing:
                    member = tar.getmember(r['relative_path']); assert member.isfile()
                    path = HERE / r['relative_path']
                    path.write_bytes(tar.extractfile(member).read()); verify(path, r)
        elif missing:
            for r in missing:
                path = HERE / r['relative_path']; temporary = path.with_suffix(path.suffix + '.part')
                download(r['source_url'], temporary); verify(temporary, r); temporary.replace(path)
    for r in rows: verify(HERE / r['relative_path'], r)
    print(str(len(rows)) + ' candidate-followup inputs verified.')
    if args.action == 'snapshot':
        path = HERE / NAME
        with tarfile.open(path, 'w', format=tarfile.USTAR_FORMAT) as tar:
            for r in rows:
                source = HERE / r['relative_path']; info = tarfile.TarInfo(r['relative_path'])
                info.size = source.stat().st_size; info.mode = 0o644; info.uid = info.gid = info.mtime = 0
                with source.open('rb') as f: tar.addfile(info, f)
        line = digest(path) + '  ' + NAME + '\n'; checksum = HERE / 'input_snapshot.sha256'
        if checksum.exists(): assert checksum.read_text() == line, 'Snapshot changed; investigate before versioning'
        else: checksum.write_text(line)
        print(line, end='')


if __name__ == '__main__': main()
