#!/usr/bin/env python3
"""Independent disorder evidence: MobiDB-lite regions from the shared InterProScan
run (whole-protein, per-residue) for the 10 targets."""
import os, json, csv
ANN = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/annotation/iprscan_top100.tsv')
BASE = os.path.dirname(os.path.abspath(__file__))
targets = json.load(open(os.path.join(BASE, 'targets.json')))
want = {v['protein_id']: g for g, v in targets.items()}
out = {g: [] for g in targets}
with open(ANN) as fh:
    for row in csv.reader(fh, delimiter='\t'):
        if len(row) < 13 or row[0] not in want:
            continue
        if row[3] == 'MobiDBLite':
            out[want[row[0]]].append({'start': int(row[6]), 'end': int(row[7]),
                                      'desc': row[5], 'score': row[8]})
json.dump(out, open(os.path.join(BASE, 'disorder_annotation_mobidb.json'), 'w'), indent=1)
print('MobiDB-lite disorder regions (note: run on first 60 aa of each protein):')
for g, v in out.items():
    tot = sum(x['end'] - x['start'] + 1 for x in v)
    print(f"  {g:28s} n_regions={len(v):2d} covered={tot:3d} aa  {v}")
