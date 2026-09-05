"""Quantify the candidate with each host reference and inspect breadth.

Each reference contains all primary CDSs for that species. Gene assignment,
alignment filtering and normalization match the preceding raw-read screen.
Positions covered are aligned match/mismatch positions, excluding deletions.
No statistical significance is inferred from these technical counts.
"""
import itertools
import math
import re
import statistics
from collections import Counter
from common import ROOT, WORK, CANDIDATES, table, write


def analyze(path, metadata):
    counts = Counter(); candidate_hits = []; quality = Counter()
    with path.open() as f:
        for read, lines in itertools.groupby(f, key=lambda s: s.split('\t', 1)[0]):
            hits = {}
            for line in lines:
                p = line.rstrip().split('\t'); m = metadata[p[5]]
                if int(p[10]) < 100 or int(p[9])/int(p[10]) < .9 or (int(p[3])-int(p[2]))/int(p[1]) < .8: continue
                tags = {s.split(':', 2)[0]: s.split(':', 2)[2] for s in p[12:]}
                group = m['group'] or 'UNGROUPED:' + m['gene']
                hit = (int(tags['AS']), p[5], int(p[7]), int(p[8]), int(p[9]), int(p[10]))
                if group not in hits or hit > hits[group][0]: hits[group] = (hit, p, tags)
            if not hits: continue
            ranked = sorted(hits, key=lambda g: (-hits[g][0][0], g))
            if len(ranked) > 1 and hits[ranked[1]][0][0] >= .98 * hits[ranked[0]][0][0]:
                quality['ambiguous'] += 1; continue
            group = ranked[0]; counts[group] += 1
            if group == 'OG0000499':
                hit, p, tags = hits[group]
                candidate_hits.append((read, p, tags))
    return counts, candidate_hits, quality


def main():
    metadata = {r['reference_id']: r for r in table(ROOT / 'transcriptomic_addendum/results/raw_reference_metadata.tsv')}
    one_to_one = {r['orthogroup'] for r in table(ROOT / 'transcriptomic_addendum/results/orthogroup_annotations.tsv') if r['one_to_one_AL_AD'] == 'True'}
    summaries, coverage_rows, bins, assignments = [], [], [], []
    for reference in ['AL', 'AD']:
        results = {sp: analyze(WORK / (sp + '_vs_' + reference + '.paf'), metadata) for sp in ['AL', 'AD']}
        eligible = sorted(g for g in one_to_one if min(results[s][0][g] for s in ['AL', 'AD']) >= 10)
        normalization = 2 ** statistics.median(math.log2(results['AL'][0][g]/results['AD'][0][g]) for g in eligible)
        a, b = [results[s][0]['OG0000499'] for s in ['AL', 'AD']]
        summaries.append(dict(reference=reference, AL_reads=a, AD_reads=b, raw_ratio=a/b if b else None,
                              background_groups=len(eligible), background_median_AL_AD=normalization,
                              normalized_ratio=a/b/normalization if b else None))
        rid = reference + '|' + CANDIDATES[reference]
        length = int(metadata[rid]['cds_length'])
        for sp in ['AL', 'AD']:
            coverage = [0] * length; starts = Counter(); bin_starts = Counter()
            for read, p, tags in results[sp][1]:
                assert p[5] == rid
                target = int(p[7]); starts[(p[4], target, int(p[8]))] += 1
                bin_starts[min(9, int(10 * target / length))] += 1
                for n, op in re.findall(r'(\d+)([MIDNSHP=X])', tags['cg']):
                    n = int(n)
                    if op in 'M=X':
                        for i in range(target, target+n): coverage[i] += 1
                    if op in 'M=XDN': target += n
                assert target == int(p[8])
                assignments.append(dict(sample=sp, reference=reference, read_id=read, start_0based=p[7],
                                        end_exclusive=p[8], strand=p[4], matching_bases=p[9], alignment_nt=p[10], cigar=tags['cg']))
            coverage_rows.append(dict(sample=sp, reference=reference, candidate_length_nt=length,
                                      reads=len(results[sp][1]), distinct_strand_start_end=len(starts),
                                      covered_nt=sum(c > 0 for c in coverage), covered_fraction=sum(c > 0 for c in coverage)/length,
                                      covered_at_5=sum(c >= 5 for c in coverage), median_depth=statistics.median(coverage),
                                      max_depth=max(coverage)))
            for k in range(10):
                segment = coverage[k*length//10:(k+1)*length//10]
                bins.append(dict(sample=sp, reference=reference, decile=k+1, start_0based=k*length//10,
                                 end_exclusive=(k+1)*length//10, covered_fraction=sum(c > 0 for c in segment)/len(segment),
                                 mean_depth=sum(segment)/len(segment), read_starts=bin_starts[k]))
    write('reference_swap_counts.tsv', summaries)
    write('candidate_read_coverage.tsv', coverage_rows)
    write('candidate_coverage_deciles.tsv', bins)
    write('candidate_read_assignments.tsv', assignments)
    print(summaries)


if __name__ == '__main__': main()
