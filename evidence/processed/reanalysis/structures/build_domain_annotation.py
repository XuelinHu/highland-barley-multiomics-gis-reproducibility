#!/usr/bin/env python3
"""Collect the domain annotation for the 10 targets from the shared annotation
tables (read-only) into domain_annotation.json."""
import os, json, csv

ANN_DIR = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/annotation')
IPR = os.path.join(ANN_DIR, 'iprscan_top100.tsv')
CDS = os.path.join(ANN_DIR, 'top100_cdsearch_hits.tsv')
BASE = os.path.dirname(os.path.abspath(__file__))

targets = json.load(open(os.path.join(BASE, 'targets.json')))
want = {}
for g, v in targets.items():
    want[v['protein_id']] = g
    want[g] = g

KEEP = {'Pfam', 'PANTHER', 'ProSiteProfiles', 'ProSitePatterns', 'SUPERFAMILY', 'Gene3D',
        'NCBIfam', 'PRINTS', 'SMART', 'CDD', 'Hamap', 'PIRSF'}
ann = {g: [] for g in targets}
seen = set()
with open(IPR) as fh:
    rd = csv.reader(fh, delimiter='\t')
    for row in rd:
        if len(row) < 13 or row[0] not in want:
            continue
        g = want[row[0]]
        analysis, sig, desc, start, stop = row[3], row[4], row[5], row[6], row[7]
        ipr, ipr_desc = row[11], row[12]
        if analysis not in KEEP:
            continue
        key = (g, analysis, sig, start, stop, ipr)
        if key in seen:
            continue
        seen.add(key)
        ann[g].append({'source': 'InterProScan', 'analysis': analysis, 'signature': sig,
                       'signature_desc': desc, 'start': int(start), 'end': int(stop),
                       'interpro': ipr, 'interpro_desc': ipr_desc})

with open(CDS) as fh:
    rd = csv.reader(fh, delimiter='\t')
    next(rd, None)
    for row in rd:
        if len(row) < 10:
            continue
        q = row[0]
        pid = q.split('>')[-1].split('[')[0] if '>' in q else q
        if pid not in want:
            continue
        g = want[pid]
        if '.CDS' in q and targets[g]['protein_id'].endswith('.CDS1') is False:
            continue
        try:
            f, t = int(row[3]), int(row[4])
        except ValueError:
            continue
        key = (g, 'CD-Search', row[7], f, t, '')
        if key in seen:
            continue
        seen.add(key)
        ann[g].append({'source': 'CD-Search', 'analysis': row[1], 'signature': row[7],
                       'signature_desc': row[8], 'start': f, 'end': t,
                       'interpro': '', 'interpro_desc': row[9].strip()[:100],
                       'evalue': row[5]})

json.dump(ann, open(os.path.join(BASE, 'domain_annotation.json'), 'w'), indent=1,
          ensure_ascii=False)
for g in targets:
    print(f"--- {g} ({len(targets[g]['seq'])} aa) ---")
    for d in ann[g]:
        print(f"   {d['source']:12s} {d['analysis']:16s} {d['signature']:22s} "
              f"{d['start']:>4}-{d['end']:<4} {d['interpro']:12s} {d['signature_desc'][:44]}")
