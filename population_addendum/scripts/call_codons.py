"""Transparent diploid codon likelihoods with fragment-level evidence."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations_with_replacement, product
import json
import math
from pathlib import Path
import pysam
from Bio.Seq import Seq
from common import OUT, WORK, json_write, table, write

CODONS = tuple(''.join(c) for c in product('ACGT', repeat=3))
GENOTYPES = tuple(combinations_with_replacement(range(64), 2))
COMPLEMENT = str.maketrans('ACGTN', 'TGCAN')
PROFILES = {'primary': (30, 10), 'sensitivity': (20, 8)}


def genotype(observations, min_depth):
    counts = Counter(o['codon'] for o in observations)
    depth = len(observations)
    if not depth:
        return dict(best_codons='', amino_acids='', gq='0.000', status='no_coverage', depth=0, codon_counts='', forward=0, reverse=0)
    likelihoods = []
    # Repeated codon/quality tuples have identical likelihood contributions.
    unique = Counter((o['codon'], tuple(o['qualities'])) for o in observations)
    for (observed, qualities), weight in unique.items():
        probabilities = []
        for allele in CODONS:
            p = 1.0
            for b, a, q in zip(observed, allele, qualities):
                e = max(0.001, 10 ** (-q/10))
                p *= 1-e if b == a else e/3
            probabilities.append(p)
        likelihoods.append((probabilities, weight))
    scores = [(sum(weight*math.log10((ps[a]+ps[b])/2) for ps,weight in likelihoods), a, b) for a,b in GENOTYPES]
    scores.sort(reverse=True)
    best, second = scores[:2]
    a, b = CODONS[best[1]], CODONS[best[2]]
    gq = min(99., 10*(best[0]-second[0]))
    forward = sum(not o['reverse'] for o in observations)
    reverse = depth-forward
    failures = []
    if depth < min_depth: failures.append('low_depth')
    if gq < 30: failures.append('low_genotype_quality')
    if a != b and (min(counts[a],counts[b]) < 3 or min(counts[a],counts[b])/depth < .2):
        failures.append('heterozygote_balance')
    if not forward or not reverse: failures.append('one_strand_only')
    return dict(best_codons='/'.join(sorted((a,b))), amino_acids='/'.join(sorted(str(Seq(c).translate()) for c in (a,b))),
                gq=f'{gq:.3f}', status=';'.join(failures) or 'PASS', depth=depth,
                codon_counts=';'.join(f'{c}:{counts[c]}' for c in sorted(counts)), forward=forward, reverse=reverse)


def read_observation(read, target, min_bq):
    if read.is_unmapped: return None, 'unmapped'
    if read.is_secondary or read.is_supplementary: return None, 'nonprimary_alignment'
    if read.is_qcfail: return None, 'qc_failed'
    if read.is_duplicate: return None, 'duplicate'
    if read.mapping_quality < 30: return None, 'low_mapping_quality'
    positions = [int(x)-1 for x in target['genomic_codon_positions_1based'].split(',')]
    aligned = {r:q for q,r in read.get_aligned_pairs() if r in positions}
    if any(aligned.get(p) is None for p in positions): return None, 'incomplete_codon'
    indices = [aligned[p] for p in positions]
    bases = ''.join(read.query_sequence[q] for q in indices)
    qualities = [read.query_qualities[q] for q in indices]
    if any(q < min_bq for q in qualities): return None, 'low_base_quality'
    if any(b not in 'ACGT' for b in bases): return None, 'ambiguous_base'
    if target['strand'] == '-': bases = bases.translate(COMPLEMENT)
    return dict(codon=bases, qualities=qualities, reverse=read.is_reverse, read1=read.is_read1,
                mapq=read.mapping_quality, start=read.reference_start+1, cigar=read.cigarstring), 'eligible'


def collapse_fragments(by_name):
    selected, conflicts = [], []
    for name, observations in sorted(by_name.items()):
        if len({o['codon'] for o in observations}) != 1:
            conflicts.append(name)
            continue
        best = max(observations, key=lambda o:(min(o['qualities']),sum(o['qualities']),o['read1']))
        selected.append(dict(best, fragment=name, eligible_mates=len(observations)))
    return selected, conflicts


def call_bam(bam_path, targets, runs, profile):
    min_bq, min_depth = PROFILES[profile]
    calls, evidence, audits = [], [], []
    with pysam.AlignmentFile(str(bam_path), 'rb') as bam:
        for target in targets:
            positions = [int(x)-1 for x in target['genomic_codon_positions_1based'].split(',')]
            observations, reasons = defaultdict(lambda:defaultdict(list)), defaultdict(Counter)
            for read in bam.fetch(target['scaffold'], min(positions), max(positions)+1):
                run = read.get_tag('RG')
                assert run in runs
                obs, reason = read_observation(read, target, min_bq)
                reasons[run][reason] += 1
                if obs: observations[run][read.query_name].append(obs)
            for run in runs:
                fragments, conflicts = collapse_fragments(observations[run])
                key = dict(run_accession=run, assembly=target['assembly'], focal_AL_position=target['focal_AL_position'],
                           role=target['role'], profile=profile)
                calls.append(dict(key, **genotype(fragments, min_depth)))
                for fragment in fragments:
                    evidence.append(dict(key, fragment=fragment['fragment'],codon=fragment['codon'],
                                         base_qualities=','.join(map(str,fragment['qualities'])),
                                         mapq=fragment['mapq'],strand='-' if fragment['reverse'] else '+',
                                         eligible_mates=fragment['eligible_mates']))
                for reason, count in sorted(reasons[run].items()): audits.append(dict(key,unit='alignment',reason=reason,count=count))
                audits.extend([dict(key,unit='fragment',reason='conflicting_mates',count=len(conflicts)),
                               dict(key,unit='fragment',reason='independent_eligible',count=len(fragments))])
    return calls,evidence,audits


def main():
    p=argparse.ArgumentParser();p.add_argument('--label',default='population');p.add_argument('--output-prefix',default='')
    a=p.parse_args(); destination=WORK/a.label
    runs=json.loads((destination/'mapping_inputs.json').read_text())['runs']
    targets=table(OUT/'target_coordinates.tsv')
    calls,evidence,audits=[],[],[]
    for assembly in ('AL_Shangrila','AD_Malaysia'):
        for profile in PROFILES:
            c,e,q=call_bam(destination/(assembly+'.bam'),[t for t in targets if t['assembly']==assembly],runs,profile)
            calls.extend(c);evidence.extend(e);audits.extend(q)
    for name,rows in [('sample_codon_calls.tsv',calls),('fragment_evidence.tsv',evidence),('genotype_filter_audit.tsv',audits)]:
        write(a.output_prefix+name,rows)
    print('Evaluated',len(calls),'sample/site/reference/threshold combinations.')

if __name__=='__main__': main()
