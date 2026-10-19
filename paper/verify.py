"""Verify that the committed paper matches its recorded inputs and output hash."""
from pathlib import Path
import hashlib
import json

ROOT = Path(__file__).resolve().parents[1]
meta = json.loads((ROOT / 'output/pdf/build.json').read_text())
for group in ['input_sha256', 'output_sha256']:
    for relative, expected in meta[group].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f'{group} mismatch: {relative}; rebuild the paper')
print('Verified paper input/output hashes for BHRD-TN-001 v'+meta['version']+'.')
