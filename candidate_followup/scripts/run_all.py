"""Run the complete candidate follow-up in a fixed order, with no downloads."""
import subprocess
import sys
from common import ROOT, HERE, WORK, OUT, MM, MP


def run(path):
    print('Running ' + str(path.relative_to(ROOT)), flush=True)
    subprocess.run([sys.executable, str(path)], check=True)


def main():
    WORK.mkdir(exist_ok=True); OUT.mkdir(exist_ok=True)
    assert subprocess.check_output([str(MM),'--version'],text=True).strip() == '2.31-r1302'
    assert subprocess.check_output([str(MP),'--version'],text=True).strip() == '0.18-r281'
    # Rebuild the complete pooled host CDS reference from the original inputs.
    # This reuses the original source code and checks its unchanged output.
    run(ROOT / 'transcriptomic_addendum/scripts/prepare_raw_reference.py')
    subprocess.run([sys.executable,str(HERE/'scripts/verify.py'),'--dependencies'],check=True)
    for name in ['audit_inventory.py','analyze_genomes.py','align_references.py',
                 'analyze_reference_swaps.py','exact_shared_markers.py','analyze_protein.py',
                 'design_variants.py','analyze_gut_proteomics.py','validate_claims.py']:
        run(HERE / 'scripts' / name)


if __name__ == '__main__': main()
