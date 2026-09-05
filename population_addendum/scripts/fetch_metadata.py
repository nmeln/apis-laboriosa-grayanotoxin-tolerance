"""Retrieve the public run/sample metadata and a pinned-tool package catalog."""
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor
import csv
import hashlib
import json
import time
import urllib.request

HERE = Path(__file__).resolve().parents[1]
ROOT = HERE.parent
DEST = HERE / 'inputs'


def download(url, path):
    path = Path(path)
    if path.exists():
        return
    for attempt in range(3):
        try:
            request = urllib.request.Request(url, headers={'User-Agent': 'bee-population-reanalysis/1.0'})
            with urllib.request.urlopen(request, timeout=60) as response:
                data = response.read()
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(data)
            print('Retrieved', path.name, len(data), flush=True)
            return
        except Exception:
            if attempt == 2:
                raise
            time.sleep(2)


def main():
    DEST.mkdir(exist_ok=True)
    fields = ('run_accession,study_accession,sample_accession,sample_alias,library_name,'
              'scientific_name,instrument_model,library_strategy,library_source,library_layout,'
              'read_count,base_count,fastq_ftp,fastq_md5,fastq_bytes,submitted_ftp')
    url = ('https://www.ebi.ac.uk/ena/portal/api/filereport?accession=PRJNA931733'
           '&result=read_run&fields=' + fields + '&format=tsv')
    records = [('ena_population_runs.tsv', url),
               ('bbmap_package_catalog.json', 'https://api.anaconda.org/package/bioconda/bbmap')]
    for name, source in records:
        download(source, DEST / name)
    with (DEST / records[0][0]).open() as f:
        runs = list(csv.DictReader(f, delimiter='\t'))
    assert len(runs) == len({r['run_accession'] for r in runs}) == 57
    assert len({r['sample_accession'] for r in runs}) == 57
    assert {sp: sum(r['scientific_name'] == sp for r in runs)
            for sp in ['Apis laboriosa', 'Apis dorsata']} == {'Apis laboriosa': 29, 'Apis dorsata': 28}
    samples = [('samples/' + r['sample_accession'] + '.xml',
                'https://www.ebi.ac.uk/ena/browser/api/xml/' + r['sample_accession']) for r in runs]
    with ThreadPoolExecutor(4) as pool:
        list(pool.map(lambda item: download(item[1], DEST / item[0]), samples))
    records += samples
    rows = []
    for name, source in sorted(records):
        data = (DEST / name).read_bytes()
        rows.append(dict(path='inputs/' + name, source_url=source,
                         size_bytes=len(data), sha256=hashlib.sha256(data).hexdigest()))
    out = HERE / 'metadata_sources.tsv'
    with out.open('w', newline='') as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t', lineterminator='\n')
        writer.writeheader()
        writer.writerows(rows)


if __name__ == '__main__':
    main()
