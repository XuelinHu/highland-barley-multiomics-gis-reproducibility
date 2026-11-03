#!/usr/bin/env python3
"""Fetch functional annotation for the top100 barley drought proteins from UniProt REST.

Why this exists: the eggNOG-mapper online service (now https://eggnog-mapper.cgmlab.org/)
requires an emailed activation link before a job will start, so no job could be run
without a working mailbox (see ANNOTATION_API_REPORT.md).  UniProtKB independently
imports eggNOG ortholog assignments, KEGG, Pfam, InterPro and GO, so the same
categories of information are retrieved here as a documented substitute.

Outputs (annotation/):
  top100_uniprot_eggnog_fallback.tsv   raw per-accession records
"""
import csv
import os
import re
import time
import urllib.error
import urllib.parse
import urllib.request

OUTDIR = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/annotation')
ANNOT = os.path.join(OUTDIR, "top100_annotation.tsv")
SWISSPROT = os.path.join(OUTDIR, "top100_swissprot.tsv")
OUT_TSV = os.path.join(OUTDIR, "top100_uniprot_eggnog_fallback.tsv")

FIELDS = "accession,protein_name,xref_eggnog,xref_kegg,xref_pfam,xref_interpro,go_id,ec,id"
REQ_TIMEOUT = 60


def log(m):
    print(f"[{time.strftime('%H:%M:%S')}] {m}", flush=True)


def collect_accessions():
    """Unique UniProt accessions + protein->accession map from existing local annotation."""
    accs = []
    prot2acc = {}
    with open(ANNOT) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            for a in re.split(r"[;,]", r.get("uniprot_ac", "") or ""):
                a = a.strip()
                if a and a not in accs:
                    accs.append(a)
    # per-protein accessions from the swissprot blast table (headerless, col0=protein id, col1=sp|acc|name)
    if os.path.exists(SWISSPROT):
        with open(SWISSPROT) as fh:
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) < 2:
                    continue
                m = re.match(r"sp\|([A-Z0-9]+)\|", f[1])
                if m:
                    prot2acc.setdefault(f[0], m.group(1))
                    if m.group(1) not in accs:
                        accs.append(m.group(1))
    return accs, prot2acc


def fetch_batch(accs):
    """Query UniProt for a batch of accessions (TSV). Retries with backoff."""
    query = " OR ".join(f"accession:{a}" for a in accs)
    url = "https://rest.uniprot.org/uniprotkb/search?" + urllib.parse.urlencode(
        {"query": f"({query})", "format": "tsv", "fields": FIELDS, "size": "500"}
    )
    last = None
    for attempt in range(5):
        try:
            with urllib.request.urlopen(url, timeout=REQ_TIMEOUT) as r:
                return r.read().decode("utf-8", "replace")
        except urllib.error.HTTPError as e:
            body = ""
            try:
                body = e.read()[:200].decode("utf-8", "replace")
            except Exception:
                pass
            last = f"HTTP {e.code}: {body}"
            if 400 <= e.code < 500 and e.code != 429:
                raise RuntimeError(f"{url[:120]} -> {last}")
        except Exception as e:  # noqa: BLE001
            last = f"{type(e).__name__}: {e}"
        wait = min(60, 5 * (2**attempt))
        log(f"  retry {attempt + 1}/5 in {wait}s ({last})")
        time.sleep(wait)
    raise RuntimeError(f"uniprot fetch failed: {last}")


def main():
    accs, prot2acc = collect_accessions()
    log(f"{len(accs)} unique accessions from local annotation tables")

    rows = []
    hdr = None
    for i in range(0, len(accs), 40):
        batch = accs[i : i + 40]
        log(f"fetching batch {i // 40 + 1} ({len(batch)} accessions) ...")
        try:
            tsv = fetch_batch(batch)
        except Exception as e:  # noqa: BLE001
            log(f"  batch failed: {e}")
            continue
        lines = [ln for ln in tsv.strip().split("\n") if ln.strip()]
        if not lines:
            continue
        if hdr is None:
            hdr = lines[0]
        for ln in lines[1:]:
            f = ln.split("\t")
            if f and f[0]:
                rows.append(f)
        time.sleep(2)

    if hdr is None:
        log("no data retrieved")
        return 1

    # de-duplicate by accession
    seen, uniq = set(), []
    for f in rows:
        if f[0] not in seen:
            seen.add(f[0])
            uniq.append(f)
    log(f"got {len(uniq)} unique UniProt entries")

    with open(OUT_TSV, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(hdr.split("\t"))
        w.writerows(uniq)
    log(f"wrote {OUT_TSV}")

    # also persist the protein -> accession map for the merge step
    map_path = os.path.join(OUTDIR, "protein_to_uniprot_ac.tsv")
    with open(map_path, "w", newline="") as out:
        w = csv.writer(out, delimiter="\t")
        w.writerow(["protein_id", "uniprot_ac"])
        for p, a in sorted(prot2acc.items()):
            w.writerow([p, a])
    log(f"wrote {map_path} ({len(prot2acc)} proteins)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
