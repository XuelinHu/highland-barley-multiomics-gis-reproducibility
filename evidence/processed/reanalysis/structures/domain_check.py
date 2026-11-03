#!/usr/bin/env python3
"""Structural consistency check vs CD-Search / InterPro domain annotation.

For every target with both a local ColabFold model and an AlphaFold-DB homolog
model: global sequence alignment (BLOSUM62) + Kabsch superposition of the
aligned CA atoms -> RMSD (does the predicted fold match the annotated family?).
Also reports mean pLDDT inside each annotated domain region on the *full-length*
local model, and the disorder evidence (low-pLDDT fraction, longest low-pLDDT run).
"""
import os, json, glob
import numpy as np
from Bio.Align import PairwiseAligner, substitution_matrices
from plddt_lib import read_plddt

BASE = os.path.dirname(os.path.abspath(__file__)); os.chdir(BASE)
targets = json.load(open('targets.json'))
afdb = json.load(open('afdb/afdb_meta.json'))
summary = {r['gene_id']: r for r in json.load(open('structure_plddt_summary.json'))}
ANN = json.load(open('domain_annotation.json'))

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G','HIS':'H',
       'ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S','THR':'T','TRP':'W',
       'TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def ca(pdb):
    d = {}
    for ln in open(pdb):
        if ln.startswith('ATOM') and ln[12:16].strip() == 'CA':
            d[int(ln[22:26])] = np.array([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
    return d

def ca_seq(pdb):
    s = []
    for ln in open(pdb):
        if ln.startswith('ATOM') and ln[12:16].strip() == 'CA':
            s.append(AA3.get(ln[17:20].strip(), 'X'))
    return ''.join(s)

al = PairwiseAligner()
al.substitution_matrix = substitution_matrices.load('BLOSUM62')
# local (Smith-Waterman): the question is 'does the target CONTAIN this fold',
# not 'are the two full chains co-linear' (barley LRR is 224 aa vs rice SERK2 628 aa)
al.open_gap_score, al.extend_gap_score, al.mode = -11, -1, 'local'

rows = []
for g in targets:
    rec = summary.get(g, {})
    pdb = rec.get('pdb')
    if not pdb or not os.path.exists(pdb):
        continue
    tca = ca(pdb)
    _, per = read_plddt(pdb)
    plddt = {k[1]: v[1] for k, v in per.items()}
    row = {'gene_id': g, 'seq_len': len(targets[g]['seq']), 'source': rec.get('source'),
           'pdb': pdb, 'mean_plddt': rec.get('mean')}
    # ---- disorder evidence ----
    vals = [plddt[k] for k in sorted(plddt)]
    low = [v < 50 for v in vals]
    best = cur = 0
    for b in low:
        cur = cur + 1 if b else 0
        best = max(best, cur)
    row['n_residues'] = len(vals)
    row['frac_plddt_lt50'] = round(sum(low) / len(vals), 3) if vals else None
    row['longest_low_plddt_run'] = best
    # ---- alignment vs the *other* structure source ----
    acs = [a for a, v in afdb.items() if v['gene'] == g]
    hp = None; ac = None; cmp_kind = ''
    if rec.get('source', '').startswith('AFDB'):
        cfp = sorted(glob.glob(f'pdb/{g}_colabfold*.pdb'))
        if cfp:                      # AFDB exact match -> cross-check vs our ColabFold run
            hp, cmp_kind, ac = cfp[0], 'ColabFold(independent AF2 run)', 'local'
    elif acs:                        # ColabFold model -> compare vs the AFDB homolog
        ac = acs[0]; hp = afdb[ac]['pdb_file']; cmp_kind = 'AFDB homolog'
    if hp:
        hca, hseq = ca(hp), ca_seq(hp)
        _, hper = read_plddt(hp)
        hplddt = {k[1]: v[1] for k, v in hper.items()}
        aln = al.align(targets[g]['seq'], hseq)[0]
        # aln.aligned[i] is an (n_blocks, 2) array of [start, end) per block, so
        # zip() pairs the two rows of the SAME block: (start1,end1) and (start2,end2).
        # The correct residue pairing walks both blocks from their STARTS.
        pairs = []
        for (ts, te), (hs, he) in zip(aln.aligned[0], aln.aligned[1]):
            assert te - ts == he - hs, 'aligned block lengths differ'
            for k in range(te - ts):
                pairs.append((ts + k + 1, hs + k + 1))
        common = [(a, b) for a, b in pairs if a in tca and b in hca]
        def kabsch(pairs):
            A = np.array([tca[a] for a, _ in pairs]); B = np.array([hca[b] for _, b in pairs])
            Ac, Bc = A - A.mean(0), B - B.mean(0)
            U, S, Vt = np.linalg.svd(Ac.T @ Bc)
            d = np.sign(np.linalg.det(Vt.T @ U.T))
            R = Vt.T @ np.diag([1, 1, d]) @ U.T
            return float(np.sqrt((((R @ Ac.T).T - Bc) ** 2).sum(1).mean()))
        if len(common) > 3:
            # "core" = residues both models call confident (pLDDT>=70); a whole-chain
            # Kabsch is dominated by disordered tails and domain motion.
            core = [(a, b) for a, b in common
                    if plddt.get(a, 0) >= 70 and hplddt.get(b, 0) >= 70]
            # best-superposing 60-residue window inside the core: for multi-domain
            # proteins and solenoids this asks whether the fold is reproduced locally
            # even when inter-domain packing differs.
            best_w = None
            for i in range(len(core)):
                w = [p for p in core if core[i][0] <= p[0] < core[i][0] + 60]
                if len(w) >= 25:
                    r = kabsch(w)
                    if best_w is None or r < best_w[0]:
                        best_w = (r, core[i][0], len(w))
            row.update({
                'ca_rmsd_core_A': round(kabsch(core), 2) if len(core) > 3 else None,
                'n_core_ca': len(core),
                'best_60aa_window_rmsd_A': round(best_w[0], 2) if best_w else None,
                'best_60aa_window_start': int(best_w[1]) if best_w else None,
                'best_60aa_window_n': best_w[2] if best_w else None,
                'cmp_kind': cmp_kind,
                'afdb_ac': ac or '',
                'afdb_organism': afdb[ac]['organismScientificName'] if ac in afdb else '',
                'afdb_desc': afdb[ac]['uniprotDescription'] if ac in afdb else 'same sequence (barley)',
                'afdb_global_plddt': afdb[ac]['globalMetricValue'] if ac in afdb else '',
                'n_aligned_ca': len(common),
                'seq_identity': round(sum(1 for a, b in common if targets[g]['seq'][a-1] == hseq[b-1]) / len(common), 3),
                'ca_rmsd_A': round(kabsch(common), 2),
                'coverage_of_target': round(len(common) / len(targets[g]['seq']), 3),
            })
    # ---- pLDDT inside annotated domain regions ----
    doms = []
    for d in ANN.get(g, []):
        s, e = max(d['start'], 1), d['end']
        v = [plddt[k] for k in range(s, e + 1) if k in plddt]
        if not v:
            continue
        doms.append({**{k: (int(d[k]) if k in ('start', 'end') else d[k])
                        for k in ('source', 'analysis', 'signature', 'signature_desc',
                                  'start', 'end', 'interpro', 'interpro_desc')},
                     'region_mean_plddt': round(float(np.mean(v)), 2),
                     'region_frac_ge70': round(sum(1 for x in v if x >= 70) / len(v), 3),
                     'region_n': int(len(v))})
    row['domain_regions'] = doms
    rows.append(row)

json.dump(rows, open('domain_consistency.json', 'w'), indent=1, ensure_ascii=False)

with open('domain_consistency.tsv', 'w') as fh:
    fh.write('gene_id\tseq_len\tsource\tmean_plddt\tn_residues\tfrac_plddt_lt50\t'
             'longest_low_plddt_run\tafdb_ac\tafdb_organism\tseq_identity\tcoverage_of_target\t'
             'n_aligned_ca\tca_rmsd_A\tn_core_ca\tca_rmsd_core_A\t'
             'best_60aa_window_rmsd_A\tbest_60aa_window_start\tafdb_global_plddt\n')
    for r in rows:
        fh.write('\t'.join(str(r.get(k, '')) for k in
                 ('gene_id', 'seq_len', 'source', 'mean_plddt', 'n_residues', 'frac_plddt_lt50',
                  'longest_low_plddt_run', 'afdb_ac', 'afdb_organism', 'seq_identity',
                  'coverage_of_target', 'n_aligned_ca', 'ca_rmsd_A', 'n_core_ca',
                  'ca_rmsd_core_A', 'best_60aa_window_rmsd_A', 'best_60aa_window_start',
                  'afdb_global_plddt')) + '\n')

for r in rows:
    out = f"{r['gene_id']:28s} mean={str(r.get('mean_plddt')):>6} <50={r['frac_plddt_lt50']:.3f} run={r['longest_low_plddt_run']:>4}"
    if 'ca_rmsd_A' in r:
        org = (r['afdb_organism'].split()[0] if r.get('afdb_organism') else 'barley-self')
        out += (f" | vs {r.get('afdb_ac')} ({org:11s}) ident={r['seq_identity']:.2f} "
                f"cov={r['coverage_of_target']:.2f} RMSD_all={r['ca_rmsd_A']:6.2f} "
                f"RMSD_core={r.get('ca_rmsd_core_A')} (n={r.get('n_core_ca')})")
    else:
        out += ' | no AFDB homolog'
    print(out)
print('\nwrote domain_consistency.json / domain_consistency.tsv')
