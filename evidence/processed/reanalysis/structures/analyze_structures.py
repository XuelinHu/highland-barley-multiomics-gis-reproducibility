#!/usr/bin/env python3
"""Extract per-residue pLDDT from every target PDB, compute stats, write TSVs.

Source priority per target: local ColabFold model (pdb/<gene>_colabfold_rank001.pdb)
> AlphaFold-DB entry (afdb/AF-<ac>-F1.pdb).
"""
import os, json, glob, csv
from plddt_lib import read_plddt, plddt_stats, pdb_first_model_last_atom_resid

BASE = os.path.dirname(os.path.abspath(__file__)); os.chdir(BASE)
targets = json.load(open('targets.json'))
LABEL = json.load(open('figure_labels.json'))
afdb = json.load(open('afdb/afdb_meta.json'))
cf_scores = json.load(open('colabfold_scores.json')) if os.path.exists('colabfold_scores.json') else {}
AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M'}

def ca_seq(pdb):
    return ''.join(AA3.get(l[17:20].strip(), 'X') for l in open(pdb)
                   if l.startswith('ATOM') and l[12:16].strip() == 'CA')

# AlphaFold-DB entry whose sequence is IDENTICAL to the target -> use it directly
exact = {}
for a, v in afdb.items():
    try:
        if ca_seq(v['pdb_file']) == targets[v['gene']]['seq']:
            exact[v['gene']] = a
    except Exception:
        pass
print('AFDB entries identical to the target sequence:', exact)

rows = []
for g in targets:
    cf = sorted(glob.glob(f'pdb/{g}_colabfold*.pdb'))
    if g in exact:
        cand, src = [afdb[exact[g]]['pdb_file']], f'AFDB:{exact[g]} (exact match)'
    elif cf:
        cand, src = cf, 'ColabFold(local)'
    else:
        acs = [a for a, v in afdb.items() if v['gene'] == g]
        cand, src = [afdb[a]['pdb_file'] for a in acs], (f'AFDB:{acs[0]}' if acs else 'FAILED')
    if not cand:
        rows.append({'gene_id': g, 'label': LABEL[g], 'source': 'FAILED',
                     'seq_len': len(targets[g]['seq']), 'n_residues': None})
        continue
    p = cand[0]
    _, per = read_plddt(p)
    vals = [per[k][1] for k in per]
    st = plddt_stats(vals)
    st.update({'gene_id': g, 'label': LABEL[g], 'source': src, 'pdb': p,
               'seq_len': len(targets[g]['seq']), 'modeled_residues': len(vals)})
    if src == 'ColabFold(local)':
        st['ptm'] = cf_scores.get(g, {}).get('ptm')
        st['colabfold_model'] = cf_scores.get(g, {}).get('model')
    else:
        st['ptm'] = None
        st['colabfold_model'] = None
        st['afdb_global_plddt_reported'] = afdb[g in exact and exact[g] or '']['globalMetricValue'] if g in exact else None
    rows.append(st)
    with open(f'plddt_{g}.tsv', 'w') as fh:
        fh.write('residue_index\tresidue_aa\tplddt\n')
        for i, k in enumerate(sorted(per, key=lambda x: x[1]), 1):
            aa, b = per[k]
            fh.write(f'{i}\t{aa}\t{b:.2f}\n')

cols = ['gene_id', 'label', 'source', 'colabfold_model', 'seq_len', 'modeled_residues',
        'mean', 'median', 'min', 'max', 'std', 'ptm',
        'fraction_very_low_lt50', 'fraction_low_50_70',
        'fraction_confident_70_90', 'fraction_very_high_ge90', 'pdb']
with open('structure_plddt_summary.tsv', 'w', newline='') as fh:
    w = csv.DictWriter(fh, fieldnames=cols, delimiter='\t', extrasaction='ignore')
    w.writeheader()
    for r in rows:
        w.writerow(r)
json.dump(rows, open('structure_plddt_summary.json', 'w'), indent=1, ensure_ascii=False)
for r in rows:
    print(f"{r['gene_id']:28s} {r['source']:18s} n={str(r.get('modeled_residues')):>4} "
          f"mean={r.get('mean')} med={r.get('median')} pTM={r.get('ptm')} "
          f"<50={r.get('fraction_very_low_lt50')}")
print('wrote structure_plddt_summary.tsv/.json + plddt_<gene>.tsv')
