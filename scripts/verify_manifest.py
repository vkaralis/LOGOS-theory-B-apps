"""Verify the SHA256 inventory for the extracted release directory."""
import hashlib,json
from pathlib import Path
root=Path(__file__).resolve().parents[1]
expected=json.loads((root/'MANIFEST_SHA256.json').read_text())
errors=[]
for rel,digest in expected.items():
    p=root/rel
    if not p.is_file() or hashlib.sha256(p.read_bytes()).hexdigest()!=digest:
        errors.append(rel)
if errors:raise SystemExit('Hash mismatch or missing file: '+', '.join(errors))
print(f'PASS: {len(expected)} release files verified')
