#!/usr/bin/env python3
"""Collect pTM / mean pLDDT of the rank_001 ColabFold model for each target."""
import os, json, glob
os.chdir(os.path.dirname(os.path.abspath(__file__)))
targets = json.load(open('targets.json'))
out = {}
for g in targets:
    f = sorted(glob.glob(f'cf_out/{g}_scores_rank_001_*.json'))
    if not f:
        continue
    d = json.load(open(f[0]))
    out[g] = {'model': f[0].split('_scores_rank_001_')[1].replace('.json',''),
              'plddt_mean': round(sum(d['plddt']) / len(d['plddt']), 2),
              'ptm': round(d.get('ptm', float('nan')), 3),
              'max_pae': round(d.get('max_pae', float('nan')), 2) if 'max_pae' in d else None,
              'n_res': len(d['plddt'])}
json.dump(out, open('colabfold_scores.json', 'w'), indent=1)
for g, v in out.items():
    print(f"{g:28s} {v['model']:38s} pLDDT={v['plddt_mean']:6.2f} pTM={v['ptm']:.3f}")
print('wrote colabfold_scores.json')
