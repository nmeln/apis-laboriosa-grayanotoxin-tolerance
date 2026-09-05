"""Count identical candidate sequence without an aligner or species preference.

Use the central 100 nt (positions 26..125) of each 150-nt R1 read. A marker
must be identical in both candidate CDSs and absent, in either orientation,
from every other pooled primary CDS and both tested BQCV genomes. This tests
a shared subset of the candidate, not total transcript abundance.
"""
import gzip
import json
from collections import Counter
from Bio import SeqIO
from Bio.Seq import Seq
from common import ROOT, OUT, CANDIDATES, write


def canonical(s):
    return min(s, str(Seq(s).reverse_complement()))


def main():
    records = list(SeqIO.parse(ROOT / 'transcriptomic_addendum/work/pooled_primary_cds.fna', 'fasta'))
    targets = {sp + '|' + pid for sp, pid in CANDIDATES.items()}
    candidate = {r.id[:2]: str(r.seq) for r in records if r.id in targets}
    assert set(candidate) == {'AL', 'AD'}
    sets = [{canonical(s[i:i+100]) for i in range(len(s)-99)} for s in candidate.values()]
    shared = sets[0] & sets[1]
    before = len(shared)
    other = [str(r.seq) for r in records if r.id not in targets]
    other += [str(SeqIO.read(ROOT / 'transcriptomic_addendum/inputs' / n, 'genbank').seq) for n in ['OR496406.gb', 'KY741959.gb']]
    # Reverse complements of markers make scanning other sequences cheaper.
    oriented = shared | {str(Seq(s).reverse_complement()) for s in shared}
    collisions = set()
    for s in other:
        for i in range(len(s)-99):
            k = s[i:i+100]
            if k in oriented: collisions.add(canonical(k))
    shared -= collisions
    counts, per_marker, excluded = {}, {}, {}
    for sp, run in [('AL', 'SRR9034695'), ('AD', 'SRR9034696')]:
        observed = Counter(); total = 0; wrong_length = 0
        with gzip.open(ROOT / 'transcriptomic_addendum/inputs' / (run + '.first_2000000.R1.fastq.gz'), 'rt') as f:
            while True:
                header = f.readline()
                if not header: break
                seq = f.readline().strip(); plus = f.readline(); quality = f.readline().strip()
                assert header.startswith('@') and plus.startswith('+') and len(seq) == len(quality)
                total += 1
                if len(seq) != 150:
                    wrong_length += 1; continue
                k = seq[25:125]
                if k in shared: observed[k] += 1
                else:
                    k = str(Seq(k).reverse_complement())
                    if k in shared: observed[k] += 1
        assert total == 2_000_000
        counts[sp] = sum(observed.values()); per_marker[sp] = observed; excluded[sp] = wrong_length
    rows = []
    for k in sorted(shared):
        positions = {}
        for sp, s in candidate.items():
            positions[sp] = ';'.join(str(i+1) for i in range(len(s)-99) if s[i:i+100] in [k, str(Seq(k).reverse_complement())])
        rows.append(dict(marker=k, AL_positions_1based=positions['AL'], AD_positions_1based=positions['AD'],
                         AL_reads=per_marker['AL'][k], AD_reads=per_marker['AD'][k]))
    write('exact_shared_markers.tsv', rows)
    summary = dict(marker_length=100, read_interval_1based='26-125', input_reads_each=2_000_000,
                   shared_markers_before_specificity_filter=before, excluded_nonspecific_markers=len(collisions),
                   specific_shared_markers=len(shared), excluded_non150nt_reads=excluded,
                   counts=counts, raw_AL_AD_ratio=counts['AL']/counts['AD'] if counts['AD'] else None,
                   limitation='Exact-match subset; these are technical counts from the same two pools, without additional biological replication.')
    (OUT / 'exact_shared_marker_summary.json').write_text(json.dumps(summary, indent=2, sort_keys=True) + '\n')
    print(summary)


if __name__ == '__main__': main()
