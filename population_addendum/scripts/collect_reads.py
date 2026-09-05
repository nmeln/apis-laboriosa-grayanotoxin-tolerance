"""Collect candidate pairs from complete, MD5-verified public FASTQ files.

Run once per accession. Downloads can resume after transient network failures.
No sample is complete until both source digests and the BBDuk input totals
agree with ENA. The selected reads retain base qualities and both mates.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import gzip
import json
from pathlib import Path
import re
import subprocess
import tempfile
import time
import urllib.request
from common import HERE, INPUT, WORK, OUT, digest, json_write, table
from setup_tools import classpath, SHA256 as TOOL_SHA256


def fetch_verified(url, path, expected_bytes, expected_md5):
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        assert path.stat().st_size == expected_bytes and digest(path, 'md5') == expected_md5
        return
    part = path.with_suffix(path.suffix + '.partial')
    for attempt in range(4):
        try:
            offset = part.stat().st_size if part.exists() else 0
            assert offset <= expected_bytes
            h = hashlib.md5()
            if offset:
                with part.open('rb') as f:
                    while block := f.read(1024 * 1024):
                        h.update(block)
            if offset < expected_bytes:
                req = urllib.request.Request(url, headers={'Range': f'bytes={offset}-'} if offset else {})
                with urllib.request.urlopen(req, timeout=90) as src:
                    if offset:
                        assert src.status == 206 and src.headers['Content-Range'].startswith(f'bytes {offset}-')
                    with part.open('ab' if offset else 'wb') as out:
                        announced = offset // (512 * 1024 * 1024)
                        while block := src.read(1024 * 1024):
                            out.write(block)
                            h.update(block)
                            offset += len(block)
                            if offset // (512 * 1024 * 1024) > announced:
                                announced = offset // (512 * 1024 * 1024)
                                print(path.name, offset, '/', expected_bytes, 'bytes', flush=True)
            assert offset == expected_bytes, (path, offset, expected_bytes)
            assert h.hexdigest() == expected_md5, ('Source MD5 mismatch', path)
            part.rename(path)
            print('Verified complete source', path.name, flush=True)
            return
        except (OSError, TimeoutError) as error:
            if attempt == 3:
                raise
            print('Retrying interrupted source transfer:', path.name, type(error).__name__, flush=True)
            time.sleep(3)


def count_fastq(path):
    count = bases = 0
    with path.open() as f:
        while name := f.readline():
            sequence, plus, quality = f.readline().rstrip(), f.readline(), f.readline().rstrip()
            assert name.startswith('@') and plus.startswith('+') and len(sequence) == len(quality)
            assert sequence and quality
            count += 1
            bases += len(sequence)
    return count, bases


def restore_original_records(source, captured, output):
    """Use capture only for names; copy exact original FASTQ records.

    BBDuk 39.91 can normalize quality scores at N bases even without trimming.
    A second sequential pass preserves every selected source byte instead.
    """
    wanted = set()
    with captured.open('rb') as f:
        while header := f.readline():
            wanted.add(header[1:].rstrip(b'\r\n'))
            assert f.readline() and f.readline() and f.readline()
    seen, reads, bases = set(), 0, 0
    opener = gzip.open if source.suffix == '.gz' else open
    with opener(source, 'rb') as f, output.open('wb') as out:
        while header := f.readline():
            sequence, plus, quality = f.readline(), f.readline(), f.readline()
            assert header.startswith(b'@') and plus.startswith(b'+') and quality
            assert len(sequence.rstrip(b'\r\n')) == len(quality.rstrip(b'\r\n'))
            reads += 1
            bases += len(sequence.rstrip(b'\r\n'))
            key = header[1:].rstrip(b'\r\n')
            if key in wanted:
                assert key not in seen, 'Duplicate source identifier'
                seen.add(key)
                out.write(header + sequence + plus + quality)
    assert seen == wanted, ('Captured identifiers missing from original source', len(wanted-seen))
    return reads, bases, len(seen)


def collect(run, delete_raw=False):
    assert re.fullmatch(r'SRR[0-9]+', run)
    rows = table(INPUT / 'ena_population_runs.tsv') if (INPUT / 'ena_population_runs.tsv').exists() else table(OUT / 'sample_inventory.tsv')
    row = next(r for r in rows if r['run_accession'] == run)
    # Large temporary sources do not need to occupy the shared artifact tree.
    rawdir = Path(tempfile.gettempdir()) / 'bee_population_raw' / run
    destination = INPUT / 'selected_reads' / run
    destination.mkdir(parents=True, exist_ok=True)
    urls = ['https://' + u.removeprefix('https://') for u in row['fastq_ftp'].split(';')]
    sizes, md5s = list(map(int, row['fastq_bytes'].split(';'))), row['fastq_md5'].split(';')
    assert len(urls) == len(sizes) == len(md5s) == 2
    sources = [rawdir / f'{run}_{mate}.fastq.gz' for mate in [1, 2]]
    with ThreadPoolExecutor(2) as pool:
        list(pool.map(lambda args: fetch_verified(*args), zip(urls, sources, sizes, md5s)))
    output = [destination / f'{run}_{mate}.fastq' for mate in [1, 2]]
    captured = [rawdir / f'{run}_{mate}.captured.fastq' for mate in [1, 2]]
    args = [f'in={sources[0]}', f'in2={sources[1]}', f'outm={captured[0]}', f'outm2={captured[1]}',
            f'ref={OUT / "capture_baits_21.fna"}', 'k=21', 'hdist=0', 'maskmiddle=f',
            'rcomp=t', 'minkmerhits=1', 'ordered=t', 'threads=2', 'overwrite=t',
            'qtrim=f', 'minlength=1', 'maxns=-1', 'usejni=f', 'pigz=f', 'unpigz=f',
            'qin=33', 'qout=33', 'changequality=f']
    command = ['java', '-Xmx1g', '-cp', str(classpath()), 'jgi.BBDuk', *args]
    log = destination / 'bbduk.log'
    with log.open('w') as f:
        subprocess.run(command, stdout=f, stderr=subprocess.STDOUT, check=True)
    content = log.read_text()
    match = re.search(r'Input:\s+(\d+)\s+reads\s+(\d+)\s+bases', content)
    assert match, content
    input_reads, input_bases = map(int, match.groups())
    assert input_reads == 2 * int(row['read_count']), (run, input_reads, row['read_count'])
    assert input_bases == int(row['base_count']), (run, input_bases, row['base_count'])
    with ThreadPoolExecutor(2) as pool:
        restoration = list(pool.map(lambda args: restore_original_records(*args), zip(sources, captured, output)))
    assert all(r[0] == int(row['read_count']) for r in restoration)
    assert sum(r[1] for r in restoration) == input_bases
    counts = [count_fastq(p) for p in output]
    assert counts[0][0] == counts[1][0] and counts[0][0] > 0, (run, counts)
    evidence = dict(run_accession=run, complete_source_files_verified=True,
                    source_files=[dict(url=u, size_bytes=n, md5=m) for u, n, m in zip(urls, sizes, md5s)],
                    input_reads=input_reads, input_bases=input_bases, selected_pairs=counts[0][0],
                    selected_files=[dict(name=p.name, size_bytes=p.stat().st_size, sha256=digest(p),
                                         reads=c[0], bases=c[1]) for p, c in zip(output, counts)],
                    bbduk_version='39.91', tool_package_sha256=TOOL_SHA256,
                    selected_records_copied_from_original_fastqs=True,
                    capture_baits_sha256=digest(OUT / 'capture_baits_21.fna'),
                    options=[a for a in args if not a.startswith(('in=', 'in2=', 'outm=', 'outm2=', 'ref='))])
    json_write(destination / 'collection.json', evidence)
    print('Completed', run, 'selected pairs:', counts[0][0], flush=True)
    if delete_raw:
        # Only exact verified source downloads produced by this invocation.
        # Their public URLs, sizes and digests remain in collection.json.
        for path in sources + captured:
            path.unlink()


if __name__ == '__main__':
    p = argparse.ArgumentParser()
    p.add_argument('--run', required=True)
    p.add_argument('--delete-raw', action='store_true', help='Remove verified large downloads after successful capture')
    a = p.parse_args()
    collect(a.run, a.delete_raw)
