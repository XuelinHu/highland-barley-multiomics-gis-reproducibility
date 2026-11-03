"""Map GPL16130 (Agilent barley) 60-mer probes to MorexV3 transcripts, same
pipeline as GPL14927 mapping. Only ControlType == 0 probes are queried."""
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

platform = pd.read_csv(ROOT / 'processed/metadata/GPL16130.tsv', sep='\t').fillna('')
platform = platform[platform.ControlType == 0]
valid = platform.SEQUENCE.str.fullmatch('[ACGTacgt]+')
probes = platform.loc[valid]
with (OUT / 'GPL16130.probes.fa').open('w') as f:
    for row in probes.itertuples():
        f.write(f'>{row.ID}\n{row.SEQUENCE}\n')
print('probes queried:', len(probes), flush=True)

cols = 'qseqid sseqid pident length qlen qcovhsp evalue bitscore'
subprocess.run([str(BIN / 'blastn'), '-task', 'blastn-short',
                '-query', str(OUT / 'GPL16130.probes.fa'),
                '-db', str(REF / 'MorexV3.cdna.fa'),
                '-num_threads', '8', '-dust', 'no', '-evalue', '1e-5',
                '-max_target_seqs', '100', '-outfmt', '6 ' + cols,
                '-out', str(OUT / 'GPL16130_probe_hits.tsv')], check=True)
hits = pd.read_csv(OUT / 'GPL16130_probe_hits.tsv', sep='\t', names=cols.split(),
                   dtype={'qseqid': str})
hits['gene_id'] = hits.sseqid.map(genes)
good = hits[(hits.pident >= 95) & (hits.qcovhsp >= 95)]
grouped = {key: group for key, group in good.groupby('qseqid')}
rows = []
for probe in platform.ID:
    h = grouped.get(str(probe), good.iloc[:0])
    gs = sorted(h.gene_id.dropna().unique())
    rows.append(dict(probe_id=str(probe),
                     status='unique_gene' if len(gs) == 1 else 'ambiguous' if gs else 'no_qualified_hit',
                     gene_ids=';'.join(gs)))
result = pd.DataFrame(rows)
result.to_csv(OUT / 'GPL16130_morex_mapping.tsv', sep='\t', index=False)
print(result.status.value_counts().to_string(), flush=True)
