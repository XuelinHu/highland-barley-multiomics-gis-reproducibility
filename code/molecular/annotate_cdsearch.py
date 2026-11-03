"""Batch NCBI CD-Search (conserved domain) annotation of the top-100 proteins."""
import time
import urllib.parse
import urllib.request
from pathlib import Path
import pandas as pd

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
OUT = ROOT / 'processed/reanalysis/annotation'
OUT.mkdir(exist_ok=True)
URL = 'https://www.ncbi.nlm.nih.gov/Structure/bwrpsb/bwrpsb.cgi'


def post(params):
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(URL, data=data)
    return urllib.request.urlopen(req, timeout=120).read().decode()


queries = (ROOT / 'processed/reanalysis/top100_proteins.fasta').read_text()
resp = post({'queries': queries, 'cddefl': 'true', 'tdata': 'hits',
             'qdefl': 'true', 'db': 'cdd', 'smode': 'auto'})
cdsid = None
for line in resp.splitlines():
    if line.startswith('#cdsid'):
        cdsid = line.split('\t')[1]
print('cdsid:', cdsid, flush=True)


def parse_status(resp):
    for line in resp.splitlines():
        if line.startswith('#status'):
            parts = line.split('\t')
            out = {}
            i = 0
            while i + 1 < len(parts):
                out[parts[i].lstrip('#')] = parts[i + 1]
                i += 2
            return out
    return {}


result = None
for attempt in range(40):                       # up to ~10 min
    time.sleep(15)
    resp = post({'cdsid': cdsid, 'qdefl': 'true', 'tdata': 'hits'})
    status = parse_status(resp)
    print('status:', status.get('status'), status.get('msg'), flush=True)
    if status.get('status') in ('4', '5'):
        result = resp
        break
if result is None:
    raise SystemExit('CD-search did not finish')

rows = []
for line in result.splitlines():
    if line.startswith('#') or not line.strip():
        continue
    parts = line.split('\t')
    if len(parts) < 10:
        continue
    # query, hits-type, pssm-id, from, to, e-value, bit-score, accession, short-name, ...
    rows.append(dict(query=parts[0], hit_type=parts[1], pssm_id=parts[2],
                     q_from=parts[3], q_to=parts[4], evalue=parts[5],
                     bitscore=parts[6], accession=parts[7], short_name=parts[8],
                     incomplete=parts[9] if len(parts) > 9 else ''))
pd.DataFrame(rows).to_csv(OUT / 'top100_cdsearch_hits.tsv', sep='\t', index=False)
print('total hit rows:', len(rows))
