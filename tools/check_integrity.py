"""Read-only SHA-256 verification of every file in the deposited supplement."""
from pathlib import Path
import hashlib
import json


def main():
    root = Path(__file__).resolve().parents[1]
    manifest = json.loads((root / 'deposit_manifest.json').read_text(encoding='utf-8'))
    errors = []
    for entry in manifest['files']:
        path = (root / entry['path']).resolve()
        if not path.is_relative_to(root) or not path.is_file():
            errors.append(entry['path'] + ': missing or invalid path')
            continue
        digest = hashlib.sha256()
        with path.open('rb') as stream:
            for block in iter(lambda: stream.read(4 * 1024 * 1024), b''):
                digest.update(block)
        if path.stat().st_size != entry['bytes'] or digest.hexdigest() != entry['sha256']:
            errors.append(entry['path'] + ': content mismatch')
    if errors:
        raise SystemExit('\n'.join(errors))
    print(f"PASS: all {len(manifest['files'])} deposited scientific files match their archived hashes.")


if __name__ == '__main__':
    main()
