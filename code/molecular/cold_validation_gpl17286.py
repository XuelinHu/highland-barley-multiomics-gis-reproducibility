"""Gene-level direction-consistency validation of GSE47882 (cold acclimation,
GPL17286) against the drought consensus genes, using the new probe->MorexV3
gene mapping.

Probes are collapsed to genes by taking the probe with the largest |logFC| per
gene.  GPL17286 probe IDs are integers in the limma tables, so both sides are
forced to str before joining.
"""
from pathlib import Path
import pandas as pd
from scipy.stats import binomtest

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
COLD = ROOT / 'processed/reanalysis/cold'
OUT = COLD / 'validation'
OUT.mkdir(exist_ok=True)

genes = pd.read_csv(ROOT / 'processed/reanalysis/drought_consensus_genes.tsv', sep='\t')
print(f'drought consensus genes: {len(genes)} (up {sum(genes.direction == "up")}, '
      f'down {sum(genes.direction == "down")})', flush=True)

mapping = pd.read_csv(ROOT / 'processed/mapping/GPL17286_morex_mapping.tsv',
                      sep='\t', dtype=str)
print('GPL17286 mapping status:\n' + mapping.status.value_counts().to_string(), flush=True)
mapping = mapping[mapping.status == 'unique_gene'].set_index('probe_id')
mapping.index = mapping.index.astype(str)  # GPL17286 IDs are ints -> force str
print(f'GPL17286 unique-gene probes: {len(mapping)}', flush=True)

rows = []
for contrast in ['d1_vs_d0', 'd7_vs_d0', 'd21_vs_d0']:
    de = pd.read_csv(COLD / f'GSE47882_{contrast}_limma.tsv', sep='\t', index_col=0)
    de = de.rename_axis('probe_id').reset_index()
    de['probe_id'] = de.probe_id.astype(str)
    de['gene_id'] = de.probe_id.map(mapping.gene_ids)
    # significant probes on the full platform, before the mapping restriction
    sig_all = int(((de['adj.P.Val'] < 0.05) & (de.logFC.abs() >= 1)).sum())
    de = de.dropna(subset=['gene_id'])

    sig = de[(de['adj.P.Val'] < 0.05) & (de.logFC.abs() >= 1)]
    sig_mapped_probes = int(len(sig))
    # one representative probe per gene: the one with the largest |logFC|
    sig = (sig.reindex(sig.logFC.abs().sort_values(ascending=False).index)
              .drop_duplicates('gene_id'))
    sig = sig.assign(direction=(sig.logFC > 0).map({True: 'up', False: 'down'}))

    overlap = sig.merge(genes[['gene_id', 'direction']], on='gene_id', how='inner',
                        suffixes=('_cold', '_drought'))
    n = len(overlap)
    same = int((overlap.direction_cold == overlap.direction_drought).sum())
    p = binomtest(same, n, 0.5) if n else None
    rows.append(dict(
        contrast=contrast,
        cold_sig_probes=sig_all,
        cold_sig_mapped_probes=sig_mapped_probes,
        cold_sig_genes=int(sig.gene_id.nunique()),
        overlap=n,
        same_direction=same,
        opposite_direction=n - same,
        same_up=int(((overlap.direction_cold == 'up') &
                     (overlap.direction_drought == 'up')).sum()),
        same_down=int(((overlap.direction_cold == 'down') &
                       (overlap.direction_drought == 'down')).sum()),
        binom_p='%.3g' % p.pvalue if p else None))
    print(f'{contrast}: cold sig genes {rows[-1]["cold_sig_genes"]}, overlap {n}, '
          f'same direction {same}' + (f', binomial p={p.pvalue:.3g}' if p else ''),
          flush=True)

pd.DataFrame(rows).to_csv(OUT / 'GSE47882_direction_validation.tsv', sep='\t', index=False)
print('written', OUT / 'GSE47882_direction_validation.tsv', flush=True)
