#!/usr/bin/env python3
"""Install the checksum-pinned upstream minimap2 Linux x86-64 binary."""
import hashlib
import platform
import subprocess
import tarfile
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
NAME = 'minimap2-2.31_x64-linux.tar.bz2'
SHA = '300bc287f05eb890c6211fa7db043ce98320a401621fadd1cfdbeabd1a6e4ab5'
URL = 'https://github.com/lh3/minimap2/releases/download/v2.31/'+NAME

def main():
    assert platform.system() == 'Linux' and platform.machine() in ['x86_64','AMD64'], 'Pinned binary requires Linux x86-64'
    folder = ROOT/'.tools'
    folder.mkdir(exist_ok=True)
    archive = folder/NAME
    if not archive.exists():
        temporary = archive.with_suffix('.part')
        urllib.request.urlretrieve(URL, temporary)
        assert hashlib.sha256(temporary.read_bytes()).hexdigest() == SHA
        temporary.replace(archive)
    assert hashlib.sha256(archive.read_bytes()).hexdigest() == SHA
    binary = folder/'minimap2-2.31_x64-linux/minimap2'
    # Extract only the executable from a checked archive. No ownership changes.
    with tarfile.open(archive) as tar:
        member = tar.getmember('minimap2-2.31_x64-linux/minimap2')
        binary.parent.mkdir(exist_ok=True)
        binary.write_bytes(tar.extractfile(member).read())
    binary.chmod(0o755)
    assert subprocess.check_output([str(binary),'--version'],text=True).strip() == '2.31-r1302'
    print('Verified minimap2 2.31-r1302')

if __name__ == '__main__': main()
