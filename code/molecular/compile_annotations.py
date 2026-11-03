"""Compile top-100 annotation: best Swiss-Prot homolog, conserved domains
(CD-Search), and drought/stress keyword check (dehydrin, LEA, ABA, etc.)."""
import pandas as pd
from pathlib import Path

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
OUT = ROOT / 'processed/reanalysis/annotation'

top = pd.read_csv(ROOT / 'processed/reanalysis/top100_drought_candidates.tsv', sep='\t')
prot2gene = {}
for _, r in top.iterrows():
    if pd.notna(r.protein_ids):
        for p in str(r.protein_ids).split(';'):
            prot2gene.setdefault(p, []).append(r.gene_id)

# --- best Swiss-Prot hit per protein ---
sp = pd.read_csv(OUT / 'top100_swissprot.tsv', sep='\t', header=None,
                 names=['qseqid', 'sseqid', 'pident', 'length', 'mismatch',
                        'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue',
                        'bitscore', 'stitle'])
sp['uniprot_ac'] = sp.sseqid.str.split('|').str[1]
sp['prot_name'] = sp.sseqid.str.split('|').str[2]
best = sp.loc[sp.groupby('qseqid').evalue.idxmin()]

# --- CD-Search domains ---
cd = pd.read_csv(OUT / 'top100_cdsearch_hits.tsv', sep='\t')
cd = cd.rename(columns={'Query': 'query', 'Hit type': 'hit_type',
                        'Short name': 'short_name', 'E-Value': 'evalue',
                        'Definition': 'definition'})
cd['query'] = cd['query'].str.extract(r'>([^\[\]]+)')[0]  # strip Q#n - >[id] wrapper
cd = cd.astype({'evalue': float})
cd = cd[cd.hit_type == 'specific']  # keep specific hits, drop superfamilies

drought_kw = ('dehydrin|LEA|late embryogenesis|ABA|abscisic|water stress|drought|'
              'cold|heat shock|HSP|aquaporin|TIP|PIP|trehalose|proline|osmotic|'
              'desiccation|wounding|jasmonate|stress')
rows = []
for _, r in top.iterrows():
    prots = str(r.protein_ids).split(';') if pd.notna(r.protein_ids) else []
    b = best[best.qseqid.isin(prots)]
    sp_name = b.stitle.min() if len(b) else ''
    sp_ac = ';'.join(b.uniprot_ac.tolist())
    hits = cd[cd['query'].isin(prots)].sort_values('evalue')
    doms = ';'.join(f"{h.short_name}({h.evalue:.1e})" for _, h in hits.iterrows())
    stress_hit = bool(pd.concat([
        b.stitle.str.contains(drought_kw, case=False, regex=True),
        hits.short_name.str.contains('dehydrin|LEA|DSP|K-seg', case=False, regex=True)
    ]).any()) if (len(b) or len(hits)) else False
    rows.append(dict(gene_id=r.gene_id, direction=r.direction,
                     n_probes=r.n_probes, probe=r.probe,
                     mean_abs_logFC=round(r.mean_abs_logFC, 3),
                     best_swissprot=sp_name, uniprot_ac=sp_ac,
                     domains=doms, stress_related=stress_hit))
ann = pd.DataFrame(rows)
ann.to_csv(OUT / 'top100_annotation.tsv', sep='\t', index=False)
print('annotated rows:', len(ann))
print('stress-related (keyword):', ann.stress_related.sum())
print('with any Swiss-Prot hit:', ann.best_swissprot.str.len().gt(0).sum())
print('with any domain:', ann.domains.str.len().gt(0).sum())
print()
print(ann[ann.stress_related][['gene_id', 'direction', 'best_swissprot', 'domains']]
      .head(15).to_string(max_colwidth=90))
