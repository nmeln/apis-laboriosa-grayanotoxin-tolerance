"""Install a checksum-pinned BBDuk Java distribution. No native JNI is used."""
import hashlib
import io
from pathlib import Path
import subprocess
import tarfile
import urllib.request
import zipfile
import zstandard

ROOT = Path(__file__).resolve().parents[2]
TOOLS = ROOT / '.tools'
PACKAGE = 'bbmap-39.91-h09cc210_0.conda'
URL = 'https://conda.anaconda.org/bioconda/linux-64/' + PACKAGE
SHA256 = '61a3354ececf3247bdcc77a67a9ba5e79e8a05b4b35ccf1455ab8bd7ce64d9de'
DEST = TOOLS / 'bbmap-39.91'


def classpath():
    matches = list(DEST.glob('opt/bbmap-*/current/jgi/BBDuk.class'))
    if len(matches) != 1:
        raise RuntimeError('Expected one installed BBDuk class: ' + repr(matches))
    return matches[0].parents[1]


def main():
    TOOLS.mkdir(exist_ok=True)
    archive = TOOLS / PACKAGE
    if not archive.exists():
        with urllib.request.urlopen(URL, timeout=60) as src, archive.open('wb') as out:
            while block := src.read(1024 * 1024):
                out.write(block)
    assert hashlib.file_digest(archive.open('rb'), 'sha256').hexdigest() == SHA256
    DEST.mkdir(exist_ok=True)
    if not list(DEST.glob('opt/bbmap-*/current/jgi/BBDuk.class')):
        with zipfile.ZipFile(archive) as z:
            name = next(n for n in z.namelist() if n.startswith('pkg-') and n.endswith('.tar.zst'))
            with z.open(name) as compressed:
                with zstandard.ZstdDecompressor().stream_reader(compressed) as raw:
                    with tarfile.open(fileobj=raw, mode='r|') as tar:
                        tar.extractall(DEST, filter='data')
    result = subprocess.run(['java', '-cp', str(classpath()), 'jgi.BBDuk', '--version'],
                            text=True, capture_output=True)
    output = result.stdout + result.stderr
    assert '39.91' in output, output
    print('Verified BBDuk 39.91 package SHA-256; Java class path:', classpath())


if __name__ == '__main__':
    main()
