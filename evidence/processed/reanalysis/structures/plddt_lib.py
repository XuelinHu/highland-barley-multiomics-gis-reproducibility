#!/usr/bin/env python3
"""Shared helpers: read PDB, extract per-residue pLDDT from B-factor column."""
import os, json

AA3 = {'ALA':'A','ARG':'R','ASN':'N','ASP':'D','CYS':'C','GLN':'Q','GLU':'E','GLY':'G',
       'HIS':'H','ILE':'I','LEU':'L','LYS':'K','MET':'M','PHE':'F','PRO':'P','SER':'S',
       'THR':'T','TRP':'W','TYR':'Y','VAL':'V','MSE':'M','SEC':'U','PYL':'O'}

def read_plddt(pdb_path):
    """Return (order_list, {resseq: (aa1, plddt)}). Uses CA atoms when present."""
    per = {}
    order = []
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(('ATOM', 'HETATM')):
                continue
            atom = line[12:16].strip()
            resn = line[17:20].strip()
            chain = line[21]
            try:
                resseq = int(line[22:26])
                b = float(line[60:66])
            except ValueError:
                continue
            if atom != 'CA':
                continue
            key = (chain, resseq)
            if key in per:
                continue
            per[key] = (AA3.get(resn, 'X'), b)
            order.append(key)
    return order, per

def plddt_stats(vals):
    import statistics
    n = len(vals)
    if n == 0:
        return {}
    srt = sorted(vals)
    def frac(lo, hi):
        return sum(1 for v in vals if lo <= v < hi) / n
    return {
        'n_residues': n,
        'mean': round(statistics.fmean(vals), 2),
        'median': round(statistics.median(vals), 2),
        'min': round(min(vals), 2),
        'max': round(max(vals), 2),
        'std': round(statistics.pstdev(vals), 2),
        'fraction_very_low_lt50': round(frac(-1, 50), 4),
        'fraction_low_50_70': round(frac(50, 70), 4),
        'fraction_confident_70_90': round(frac(70, 90), 4),
        'fraction_very_high_ge90': round(sum(1 for v in vals if v >= 90) / n, 4),
    }

def pdb_first_model_last_atom_resid(pdb_path):
    """For ColabFold multi-model PDBs, keep MODEL 1 only; return new text."""
    lines = open(pdb_path).read().splitlines()
    if not any(l.startswith('MODEL') for l in lines):
        return "\n".join(lines) + "\n"
    out, in_m1, seen = [], False, False
    for l in lines:
        if l.startswith('MODEL'):
            if seen:
                break
            in_m1 = True
            seen = True
            continue
        if l.startswith('ENDMDL'):
            break
        if in_m1:
            out.append(l)
    return "\n".join(out) + "\n"
