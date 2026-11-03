"""Drought consensus: probe-level intersection of the two reanalysed datasets.

FDR < 0.05 and |logFC| >= 1 in at least one contrast per dataset; direction
consistency required between datasets; probes mapped to MorexV3 via the
existing blastn mapping. Outputs land in processed/reanalysis (v2 results).
"""
import pandas as pd
from pathlib import Path

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
OUT = ROOT / 'processed/reanalysis'
map_ = pd.read_csv(ROOT / 'processed/mapping/GPL14927_morex_mapping.tsv', sep='\t')
map_ = map_.set_index('probe_id')
if 'protein_ids' not in map_.columns:  # old mapping file lacked the protein column
    import gzip, re
    proteins = {}
    with gzip.open(ROOT / 'reference/Hordeum_vulgare.MorexV3_pseudomolecules_assembly.pep.all.fa.gz', 'rt') as f:
        for line in f:
            if line.startswith('>'):
                gene = re.search(r'gene:(\S+)', line).group(1)
                proteins.setdefault(gene, []).append(line[1:].split()[0])
    def to_proteins(s):
        if not isinstance(s, str):
            return ''
        return ';'.join(sorted({p for g in s.split(';') for p in proteins.get(g, [])}))
    map_['protein_ids'] = map_.gene_ids.map(to_proteins)
    map_.reset_index().to_csv(ROOT / 'processed/mapping/GPL14927_morex_mapping.tsv', sep='\t', index=False)
    map_ = map_.reset_index().set_index('probe_id')

contrasts = ['CamB_leaf', 'CamB_root', 'Maresi_leaf', 'Maresi_root']


def load_dataset(acc):
    frames = {}
    for i, c in enumerate(contrasts, start=1):
        t = pd.read_csv(OUT / f'{acc}_{i}_limma.tsv', sep='\t')
        t['contrast'] = c
        frames[c] = t
    all_t = pd.concat(frames, ignore_index=True)
    probe = all_t.groupby('ProbeName').agg(
        mean_logFC=('logFC', 'mean'),
        min_padj=('adj.P.Val', 'min'),
        max_abs_logFC=('logFC', lambda s: s.abs().max()),
        n_sig=('logFC', lambda s: ((s.abs() >= 1) & (all_t.loc[s.index, 'adj.P.Val'] < 0.05)).sum()),
        seq=('Sequence', 'first'),
    )
    sig = (all_t['adj.P.Val'] < 0.05) & (all_t['logFC'].abs() >= 1)
    probe['significant'] = probe.index.isin(all_t.loc[sig, 'ProbeName'].unique())
    return probe


gse103 = load_dataset('GSE103278')
gse128 = load_dataset('GSE128048')
print('significant probes:', gse103.significant.sum(), gse128.significant.sum())

both = pd.DataFrame({'gse103_logFC': gse103.mean_logFC, 'gse128_logFC': gse128.mean_logFC})
both = both[gse103.significant & gse128.significant]
both['same_direction'] = (both.gse103_logFC * both.gse128_logFC) > 0
consensus = both[both.same_direction].copy()
consensus['direction'] = 'up' if consensus.gse103_logFC.mean() > 0 else 'down'
consensus['direction'] = consensus.gse103_logFC.map(lambda v: 'up' if v > 0 else 'down')
consensus['mean_abs_logFC'] = (consensus.gse103_logFC.abs() + consensus.gse128_logFC.abs()) / 2
consensus = consensus.sort_values('mean_abs_logFC', ascending=False)
print('intersection (significant in both):', len(both),
      '| direction-consistent:', len(consensus))

ann = consensus.join(map_[['status', 'gene_ids', 'transcript_ids', 'protein_ids']], how='left')
ann.to_csv(OUT / 'drought_consensus_probes.tsv', sep='\t')
print(ann.status.value_counts(dropna=False).to_string())

# gene-level view: one row per unique mapped gene, best probe as representative
uniq = ann[ann.status == 'unique_gene'].copy()
gene_rows = []
for gene, grp in uniq.groupby('gene_ids'):
    best = grp.loc[grp.mean_abs_logFC.idxmax()]
    gene_rows.append(dict(
        gene_id=gene, n_probes=len(grp),
        probe=best.name, direction=best.direction,
        gse103_logFC=best.gse103_logFC, gse128_logFC=best.gse128_logFC,
        mean_abs_logFC=best.mean_abs_logFC, protein_ids=best.protein_ids))
genes = pd.DataFrame(gene_rows).sort_values('mean_abs_logFC', ascending=False)
genes.to_csv(OUT / 'drought_consensus_genes.tsv', sep='\t', index=False)
print('unique consensus genes:', len(genes), '| up:', (genes.direction == 'up').sum(),
      '| down:', (genes.direction == 'down').sum())

top100 = genes.head(100)
top100.to_csv(OUT / 'top100_drought_candidates.tsv', sep='\t', index=False)
print('top100 saved; direction counts:', top100.direction.value_counts().to_dict())
