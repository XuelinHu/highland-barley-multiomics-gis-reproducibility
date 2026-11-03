"""Validate drought-consensus genes against the cold datasets (within platform).

GSE41513 (GPL16130) is the only cold platform with a probe->MorexV3 gene
mapping (60-mer blastn), so direction-consistency statistics are computed on
it. GSE27821/GSE27822 (GPL1340) and GSE47882 (GPL17286) lack gene-level
mappings and are reported at probe level only.
"""
import pandas as pd
from pathlib import Path
from scipy.stats import binomtest

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
COLD = ROOT / 'processed/reanalysis/cold'
OUT = COLD / 'validation'
OUT.mkdir(exist_ok=True)

genes = pd.read_csv(ROOT / 'processed/reanalysis/drought_consensus_genes.tsv', sep='\t')
print(f'drought consensus genes: {len(genes)} (up {sum(genes.direction == "up")}, '
      f'down {sum(genes.direction == "down")})')

mapping = pd.read_csv(ROOT / 'processed/mapping/GPL16130_morex_mapping.tsv', sep='\t')
mapping = mapping[mapping.status == 'unique_gene'].set_index('probe_id')
mapping.index = mapping.index.astype('int64')  # GSE41513 matrix IDs are numeric
print(f'GPL16130 unique-gene probes: {len(mapping)}')

rows = []
for grp, label in [('shoot_3H', 'shoot 3H vs 0H'), ('shoot_6H', 'shoot 6H vs 0H'),
                   ('root_3H', 'root 3H vs 0H'), ('root_6H', 'root 6H vs 0H')]:
    de = pd.read_csv(COLD / f'GSE41513_{grp}_limma.tsv', sep='\t', index_col=0)
    de = de.rename_axis('probe_id').reset_index()
    de['probe_id'] = de.probe_id.astype('int64')
    de['gene_id'] = de.probe_id.map(mapping.gene_ids)
    de = de.dropna(subset=['gene_id'])
    sig = de[(de['adj.P.Val'] < 0.05) & (de.logFC.abs() >= 1)]  # same cutoffs as drought consensus
    overlap = sig.merge(genes, on='gene_id', how='inner')
    same = sum(overlap.direction == (overlap.logFC > 0).map({True: 'up', False: 'down'}))
    n = len(overlap)
    p = binomtest(same, n, 0.5) if n else None
    rows.append(dict(contrast=label, cold_sig_mapped_genes=len(sig.gene_id.unique()),
                     drought_genes_tested=len(sig.gene_id.unique()),
                     overlap=n, same_direction=same,
                     binom_p=round(p.pvalue, 4) if p else None))
    print(f'{label}: cold-DE genes {len(sig.gene_id.unique())}, overlap with drought '
          f'consensus {n}, same direction {same}'
          + (f', binomial p={p.pvalue:.3g}' if p else ''))

pd.DataFrame(rows).to_csv(OUT / 'GSE41513_direction_validation.tsv', sep='\t', index=False)

# probe-level reports for the platforms without gene mapping
for acc, cs in {
    'GSE27821': ['Shorttermcold_vs_Nonvernalized', 'Vernalized_vs_Nonvernalized',
                 'X1daypostvernalized_vs_Vernalized'],
    'GSE27822': ['Postvernalized_vs_Nonvernalized'],
    'GSE47882': ['d1_vs_d0', 'd7_vs_d0', 'd21_vs_d0']}.items():
    for c in cs:
        de = pd.read_csv(COLD / f'{acc}_{c}_limma.tsv', sep='\t', index_col=0)
        n = ((de["adj.P.Val"] < 0.05) & (de.logFC.abs() >= 1)).sum()
        print(f'{acc} {c}: probe-level DE (FDR<0.05 & |logFC|>=1) = {n}')
print('validation summary written to', OUT)
