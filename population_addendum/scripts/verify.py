"""Compare generated results or prior-stage dependencies with frozen hashes."""
import argparse
from common import HERE,ROOT,OUT,digest


def main():
    p=argparse.ArgumentParser();p.add_argument('--dependencies',action='store_true');a=p.parse_args()
    manifest=HERE/('dependencies.sha256' if a.dependencies else 'results.sha256')
    root=ROOT if a.dependencies else HERE
    names=[]
    for line in manifest.read_text().splitlines():
        expected,name=line.split(maxsplit=1);names.append(name)
        assert digest(root/name)==expected,'Checksum mismatch: '+name
    if not a.dependencies:
        actual={str(p.relative_to(HERE)) for p in OUT.iterdir() if p.is_file()}
        assert actual==set(names),'Unexpected or missing result files: '+repr(actual^set(names))
    print('Verified',len(names),'dependency files.' if a.dependencies else 'result files.')

if __name__=='__main__':main()
