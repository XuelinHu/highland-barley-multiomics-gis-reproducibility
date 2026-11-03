"""Map all GPL14927 probe sequences to MorexV3 transcripts, retaining ambiguity."""
from pathlib import Path
import gzip, re, shutil, subprocess
import pandas as pd

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
REF = ROOT / 'reference'
OUT = ROOT / 'processed/mapping'
OUT.mkdir(exist_ok=True)
BIN = ROOT / 'tools/ncbi-blast-2.17.0+/bin'
name = 'Hordeum_vulgare.MorexV3_pseudomolecules_assembly.cdna.all.fa.gz'
url = 'https://ftp.ensemblgenomes.ebi.ac.uk/pub/plants/release-60/fasta/hordeum_vulgare/cdna/' + name
subprocess.run(['wget', '-c', '-q', '--timeout=60', '--tries=4', '-P', str(REF), url], check=True)
with gzip.open(REF / name, 'rb') as src, (REF / 'MorexV3.cdna.fa').open('wb') as dst:
    shutil.copyfileobj(src, dst)
genes = {}
for line in (REF / 'MorexV3.cdna.fa').open():
    if line.startswith('>'):
        genes[line[1:].split()[0]] = re.search(r'gene:(\S+)', line).group(1)
platform = pd.read_csv(ROOT / 'processed/metadata/GPL14927.tsv', sep='\t').fillna('')
valid = platform.SEQUENCE.str.fullmatch('[ACGTacgt]+')
with (OUT / 'GPL14927.probes.fa').open('w') as f:
    for row in platform.loc[valid].itertuples():
        f.write(f'>{row.ID}\n{row.SEQUENCE}\n')
subprocess.run([str(BIN / 'makeblastdb'), '-in', str(REF / 'MorexV3.cdna.fa'), '-dbtype', 'nucl'], check=True, stdout=subprocess.DEVNULL)
cols = 'qseqid sseqid pident length qlen qcovhsp evalue bitscore'
subprocess.run([str(BIN / 'blastn'), '-task', 'blastn-short', '-query', str(OUT / 'GPL14927.probes.fa'), '-db', str(REF / 'MorexV3.cdna.fa'), '-num_threads', '8', '-dust', 'no', '-evalue', '1e-5', '-max_target_seqs', '100', '-outfmt', '6 ' + cols, '-out', str(OUT / 'probe_hits.tsv')], check=True)
hits = pd.read_csv(OUT / 'probe_hits.tsv', sep='\t', names=cols.split())
hits['gene_id'] = hits.sseqid.map(genes)
good = hits[(hits.pident >= 95) & (hits.qcovhsp >= 95)]
proteins = {}
with gzip.open(REF / 'Hordeum_vulgare.MorexV3_pseudomolecules_assembly.pep.all.fa.gz', 'rt') as f:
    for line in f:
        if line.startswith('>'):
            gene = re.search(r'gene:(\S+)', line).group(1)
            proteins.setdefault(gene, []).append(line[1:].split()[0])
grouped = {key: group for key, group in good.groupby('qseqid')}
rows = []
for probe in platform.ID:
    h = grouped.get(probe, good.iloc[:0])
    gs = sorted(h.gene_id.dropna().unique())
    rows.append(dict(probe_id=probe, status='unique_gene' if len(gs)==1 else 'ambiguous' if gs else 'no_qualified_hit', gene_ids=';'.join(gs), transcript_ids=';'.join(sorted(h.sseqid.unique()))))
result = pd.DataFrame(rows)
result['protein_ids'] = result.gene_ids.map(lambda s: ';'.join(sorted({p for g in s.split(';') for p in proteins.get(g, [])})))
result.to_csv(OUT / 'GPL14927_morex_mapping.tsv', sep='\t', index=False)
print(result.status.value_counts().to_string(), flush=True)
