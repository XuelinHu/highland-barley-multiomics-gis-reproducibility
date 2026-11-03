"""Recover deposited platform tables and sample expression, preserving native scales."""
from pathlib import Path
import GEOparse

root = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed'))
for platform in ['GPL1340', 'GPL16130', 'GPL17286']:
    g = GEOparse.get_GEO(geo=platform, destdir=str(root / 'metadata'), silent=True)
    g.table.to_csv(root / 'metadata' / f'{platform}.tsv', sep='\t', index=False)
    print(platform, g.table.shape, list(g.table.columns), flush=True)
for acc in ['GSE27821', 'GSE27822', 'GSE41513', 'GSE47882']:
    g = GEOparse.get_GEO(filepath=str(root / 'metadata' / f'{acc}_family.soft.gz'), silent=True)
    matrix = g.pivot_samples('VALUE')
    matrix.to_csv(root / 'expression' / f'{acc}_deposited.tsv.gz', sep='\t')
    print(acc, matrix.shape, flush=True)
