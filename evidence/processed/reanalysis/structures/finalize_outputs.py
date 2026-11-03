#!/usr/bin/env python3
"""Copy each target's best ColabFold model into pdb/ with clean names; write PDB_LIST.tsv."""
import os, json, glob
from plddt_lib import pdb_first_model_last_atom_resid

BASE = os.path.dirname(os.path.abspath(__file__)); os.chdir(BASE)
targets = json.load(open('targets.json'))
afdb = json.load(open('afdb/afdb_meta.json'))
os.makedirs('pdb', exist_ok=True)

rows = []
for g in targets:
    src = sorted(glob.glob(f'cf_out/{g}_unrelaxed_*rank_001*.pdb')) or \
          sorted(glob.glob(f'cf_out/{g}*rank_001*.pdb'))
    dst, model = '', ''
    if src:
        dst = f'pdb/{g}_colabfold_rank001.pdb'
        open(dst, 'w').write(pdb_first_model_last_atom_resid(src[0]))
        model = os.path.basename(src[0])
    acs = [a for a, v in afdb.items() if v['gene'] == g]
    rows.append({'gene_id': g, 'protein_id': targets[g]['protein_id'],
                 'length': len(targets[g]['seq']), 'colabfold_pdb': dst,
                 'colabfold_model': model, 'afdb_ac': acs[0] if acs else '',
                 'afdb_pdb': afdb[acs[0]]['pdb_file'] if acs else ''})

with open('PDB_LIST.tsv', 'w') as fh:
    fh.write('gene_id\tprotein_id\tlength\tcolabfold_pdb\tcolabfold_model\tafdb_ac\tafdb_pdb\n')
    for r in rows:
        fh.write('\t'.join(str(r[k]) for k in
                 ('gene_id', 'protein_id', 'length', 'colabfold_pdb', 'colabfold_model',
                  'afdb_ac', 'afdb_pdb')) + '\n')
print(open('PDB_LIST.tsv').read())
