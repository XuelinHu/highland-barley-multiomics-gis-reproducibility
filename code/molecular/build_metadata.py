from pathlib import Path
import re, csv

src = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'geo_raw'))
out = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/metadata'))
out.mkdir(parents=True, exist_ok=True)
for page in src.glob('*_page.txt'):
    acc = page.stem.replace('_page', '')
    rows = {}
    for line in page.read_text(errors='ignore').splitlines():
        m = re.match(r'!Sample_(\S+)\s*=\s*(.*)', line)
        if m:
            key, val = m.groups()
            rows.setdefault(key, []).append(val.strip())
    keys = sorted(rows)
    n = max((len(v) for v in rows.values()), default=0)
    with (out / f'{acc}_samples.tsv').open('w', newline='') as f:
        w = csv.writer(f, delimiter='\t')
        w.writerow(['dataset'] + keys)
        for i in range(n):
            w.writerow([acc] + [rows[k][i] if i < len(rows[k]) else '' for k in keys])
    print(acc, n)
