"""Fetch EST clone sequences for the GPL17286 (IPK Hordeum vulgare 12k SEED2)
cDNA array from NCBI nuccore.

The platform annotation only carries EST clone names (HU04A01, HA01A05, ...),
so sequences are recovered by searching nuccore titles for each clone name in
OR-batches, then efetching the fasta records and parsing the clone name back
out of the defline.  Politely rate limited (<=3 requests/s, no api_key).
"""
from pathlib import Path
import re, sys, time, json, hashlib, urllib.parse, urllib.request

ROOT = Path(__import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), ''))
META = ROOT / 'processed/metadata/GPL17286.tsv'
OUT = ROOT / 'processed/mapping'
OUT.mkdir(parents=True, exist_ok=True)
CACHE = Path('/tmp/gpl17286_cache')
CACHE.mkdir(exist_ok=True)

EUTILS = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils'
BATCH = 100          # clone names per esearch OR-query
FETCH_BATCH = 500    # accessions per efetch
SLEEP = 0.4          # >= 3 req/s limit without an api_key


def eutils_post(endpoint, params):
    """POST to an eutils endpoint so long OR-queries stay out of the URL."""
    data = urllib.parse.urlencode(params).encode()
    req = urllib.request.Request(f'{EUTILS}/{endpoint}', data=data)
    for attempt in range(4):
        try:
            with urllib.request.urlopen(req, timeout=120) as fh:
                return fh.read().decode('utf-8', 'replace')
        except Exception as exc:                      # transient NCBI hiccups
            if attempt == 3:
                raise
            wait = 2 ** attempt
            print(f'  retry {endpoint} after {wait}s ({exc})', flush=True)
            time.sleep(wait)
    raise RuntimeError('unreachable')


def load_clone_ids():
    rows = [line.split('\t') for line in META.read_text().splitlines()[1:] if line.strip()]
    clones = sorted({r[2] for r in rows if len(r) > 2 and r[2]})
    return clones


def esearch_batches(clones):
    """Return {uid} found for the clone names, via cached OR-batch esearches."""
    all_ids = set()
    n_batches = (len(clones) + BATCH - 1) // BATCH
    for i in range(0, len(clones), BATCH):
        batch = clones[i:i + BATCH]
        tag = f'search_{i // BATCH:04d}'
        cache_file = CACHE / f'{tag}.json'
        if cache_file.exists():
            ids = json.loads(cache_file.read_text())
        else:
            term = ' OR '.join(f'"{c}"[Title]' for c in batch)
            raw = eutils_post('esearch.fcgi', dict(db='nuccore', term=term,
                                                   retmax=100000, retmode='json'))
            ids = json.loads(raw)['esearchresult']['idlist']
            cache_file.write_text(json.dumps(ids))
            time.sleep(SLEEP)
        all_ids.update(ids)
        if (i // BATCH) % 20 == 0:
            print(f'  esearch batch {i // BATCH + 1}/{n_batches}: '
                  f'{len(all_ids)} uids so far', flush=True)
    return sorted(all_ids, key=int)


CLONE_RE = re.compile(r'cDNA clone (\S+)')


def efetch_sequences(uids):
    """Fetch fasta for the uids, returning {(accession, clone_name): sequence}.

    Cache files are keyed by a hash of the uid batch, not by its position: a
    position-keyed cache silently served a pilot run's records to a later run
    whose uid list differed.
    """
    records = {}
    n_batches = (len(uids) + FETCH_BATCH - 1) // FETCH_BATCH
    for i in range(0, len(uids), FETCH_BATCH):
        batch = uids[i:i + FETCH_BATCH]
        key = hashlib.md5(','.join(batch).encode()).hexdigest()[:12]
        cache_file = CACHE / f'fasta_{i // FETCH_BATCH:04d}_{key}.txt'
        if cache_file.exists():
            text = cache_file.read_text()
        else:
            text = eutils_post('efetch.fcgi', dict(db='nuccore', id=','.join(batch),
                                                  rettype='fasta', retmode='text'))
            cache_file.write_text(text)
            time.sleep(SLEEP)
        header, chunks = None, []
        def flush(header, chunks):
            if header is None:
                return
            acc = header[1:].split()[0]
            match = CLONE_RE.search(header)
            if match and 'Hordeum vulgare' in header:
                records[(acc, match.group(1))] = ''.join(chunks)
        for line in text.splitlines():
            if line.startswith('>'):
                flush(header, chunks)
                header, chunks = line, []
            elif line.strip():
                chunks.append(line.strip())
        flush(header, chunks)
        if (i // FETCH_BATCH) % 10 == 0:
            print(f'  efetch batch {i // FETCH_BATCH + 1}/{n_batches}: '
                  f'{len(records)} clone records so far', flush=True)
    return records


def main():
    clones = load_clone_ids()
    print(f'unique CLONE_IDs on GPL17286: {len(clones)}', flush=True)

    uids = esearch_batches(clones)
    print(f'total nuccore uids retrieved: {len(uids)}', flush=True)

    records = efetch_sequences(uids)
    print(f'parsed cDNA clone records: {len(records)}', flush=True)

    # keep the longest sequence per clone (5' and 3' ESTs exist for many clones)
    best = {}
    for (acc, clone), seq in records.items():
        if clone not in best or len(seq) > len(best[clone][1]):
            best[clone] = (acc, seq)

    clone_set = set(clones)
    matched = {c: v for c, v in best.items() if c in clone_set}
    print(f'clones matched to >=1 sequence: {len(matched)}/{len(clones)}', flush=True)

    # tolerate a trailing suffix such as HU04A01u / HA01A05f seen in some titles
    unmatched = clone_set - set(matched)
    extra = {}
    for clone, (acc, seq) in best.items():
        base = re.sub(r'[a-z]$', '', clone)
        if base in unmatched and base not in extra:
            extra[base] = (acc, seq)
    print(f'extra matches after stripping 1-char suffix: {len(extra)}', flush=True)
    matched.update(extra)
    print(f'total clones with sequence: {len(matched)}/{len(clones)}', flush=True)

    with (OUT / 'GPL17286_clones.fa').open('w') as fh:
        for clone in sorted(matched):
            acc, seq = matched[clone]
            fh.write(f'>{clone} {acc} len={len(seq)}\n')
            for j in range(0, len(seq), 60):
                fh.write(seq[j:j + 60] + '\n')

    (OUT / 'GPL17286_clone_accessions.tsv').write_text(
        'clone_id\taccession\tlength\n'
        + ''.join(f'{c}\t{matched[c][0]}\t{len(matched[c][1])}\n' for c in sorted(matched)))
    print('written', OUT / 'GPL17286_clones.fa', flush=True)


if __name__ == '__main__':
    main()
