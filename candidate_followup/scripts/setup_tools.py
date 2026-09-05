"""Install the checksum-pinned official miniprot 0.18 Linux executable."""
import hashlib
import platform
import subprocess
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAME = 'miniprot-0.18_x64-linux.tar.bz2'
SHA = '794397d918c85cd55f42cd905350abe0528ceb09f6a1430c3a3d6650309ac28e'
URL = 'https://github.com/lh3/miniprot/releases/download/v0.18/' + NAME


def main():
    assert platform.system() == 'Linux' and platform.machine() in ['x86_64','AMD64']
    folder = ROOT / '.tools'; folder.mkdir(exist_ok=True)
    archive = folder / NAME
    if not archive.exists():
        temporary = archive.with_suffix('.part')
        with urllib.request.urlopen(URL, timeout=120) as response:
            temporary.write_bytes(response.read())
        assert hashlib.sha256(temporary.read_bytes()).hexdigest() == SHA
        temporary.replace(archive)
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == SHA
    binary = folder / 'miniprot-0.18_x64-linux/miniprot'
    with tarfile.open(archive) as tar:
        member = tar.getmember('miniprot-0.18_x64-linux/miniprot')
        assert member.isfile()
        binary.parent.mkdir(exist_ok=True)
        binary.write_bytes(tar.extractfile(member).read())
    binary.chmod(0o755)
    assert subprocess.check_output([str(binary), '--version'], text=True).strip() == '0.18-r281'
    print('Verified miniprot 0.18-r281')


if __name__ == '__main__': main()
