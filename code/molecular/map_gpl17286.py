"""Map GPL17286 (IPK Hordeum vulgare 12k SEED2 cDNA array) spotted EST clones
to MorexV3 genes.

The platform annotation only carries EST clone names, so clone sequences were
recovered from NCBI nuccore (fetch_gpl17286_clones.py) and are blastn'ed
against MorexV3 cDNA here, mirroring the GPL14927/GPL16130 mapping logic.
Mapping is at probe level because the same clone is spotted more than once.

Note: GPL17286 probe IDs are integers, so probe_id is written as str to avoid
int/str join failures downstream.
"""
from pathlib import Path
import re, subprocess
import pandas as pd

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
REF = ROOT / 'reference'
OUT = ROOT / 'processed/mapping'
BIN = ROOT / 'tools/ncbi-blast-2.17.0+/bin'

genes = {}
for line in (REF / 'MorexV3.cdna.fa').open():
    if line.startswith('>'):
        genes[line[1:].split()[0]] = re.search(r'gene:(\S+)', line).group(1)
print('MorexV3 transcripts:', len(genes), flush=True)

platform = pd.read_csv(ROOT / 'processed/metadata/GPL17286.tsv', sep='\t',
                       dtype={'ID': str, 'POSITION': str, 'CLONE_ID': str}).fillna('')
print('probes on platform:', len(platform), flush=True)

cols = 'qseqid sseqid pident length qlen qcovhsp evalue bitscore'
subprocess.run([str(BIN / 'blastn'),
                '-query', str(OUT / 'GPL17286_clones.fa'),
                '-db', str(REF / 'MorexV3.cdna.fa'),
                '-num_threads', '8', '-evalue', '1e-10',
                '-max_target_seqs', '10', '-outfmt', '6 ' + cols,
                '-out', str(OUT / 'GPL17286_clone_hits.tsv')], check=True)
hits = pd.read_csv(OUT / 'GPL17286_clone_hits.tsv', sep='\t', names=cols.split(),
                   dtype={'qseqid': str})
hits['gene_id'] = hits.sseqid.map(genes)
print('raw blast hits:', len(hits), flush=True)

good = hits[(hits.pident >= 95) & (hits.qcovhsp >= 95)]
print('qualified hits (pident>=95, qcovhsp>=95):', len(good), flush=True)
grouped = {key: group for key, group in good.groupby('qseqid')}
none = good.iloc[:0]

rows = []
for row in platform.itertuples():
    h = grouped.get(row.CLONE_ID, none)
    gs = sorted(h.gene_id.dropna().unique())
    ts = sorted(h.sseqid.dropna().unique())
    rows.append(dict(probe_id=str(row.ID),
                     clone_id=row.CLONE_ID,
                     status=('unique_gene' if len(gs) == 1 else
                             'ambiguous' if gs else 'no_qualified_hit'),
                     gene_ids=';'.join(gs),
                     transcript_ids=';'.join(ts)))
result = pd.DataFrame(rows)
result.to_csv(OUT / 'GPL17286_morex_mapping.tsv', sep='\t', index=False)

print(result.status.value_counts().to_string(), flush=True)
clones_total = platform.CLONE_ID.nunique()
clones_mapped = platform.loc[result.status != 'no_qualified_hit', 'CLONE_ID'].nunique()
print(f'clones mapped: {clones_mapped}/{clones_total}', flush=True)
