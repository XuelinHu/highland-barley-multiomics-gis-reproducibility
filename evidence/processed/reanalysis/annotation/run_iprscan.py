#!/usr/bin/env python3
"""Submit top100 barley proteins to InterProScan 5 (EBI Job Dispatcher REST API).

NOTE: the EBI endpoint layout is the newer Job Dispatcher v2 form:
    POST /run                  (form fields: email, title, stype, sequence)
    GET  /status/{jobId}
    GET  /result/{jobId}/tsv
(The older /submit + /run/{id}/status paths now return 404/405.)

Splits 107 sequences into 4 jobs (<=30 seqs each), polls, downloads TSV,
merges to iprscan_top100.tsv.  Retries with backoff, per-request timeout 60s.
"""
import json
import os
import sys
import time
import urllib.error
import urllib.parse
import urllib.request

BASE = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
FASTA = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/top100_proteins.fasta')
OUTDIR = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/annotation')
STATE = os.path.join(OUTDIR, "iprscan_jobs_state.json")
# reserved placeholder domain; a@b.com is rejected by EBI's validator
EMAIL = "placeholder@example.org"
CHUNK = 30
JOB_TIMEOUT = 2 * 3600  # 2 hours max per job
REQ_TIMEOUT = 60

os.makedirs(OUTDIR, exist_ok=True)


def log(msg):
    print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


def read_fasta(path):
    seqs, name, buf = [], None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs.append((name, "".join(buf)))
                name, buf = line[1:].strip(), []
            elif line.strip():
                buf.append(line.strip())
    if name is not None:
        seqs.append((name, "".join(buf)))
    return seqs


def http(url, data=None, method=None):
    """HTTP request with retries + exponential backoff. Returns bytes or raises."""
    last = None
    for attempt in range(5):
        try:
            req = urllib.request.Request(url, data=data, method=method)
            if data is not None:
                req.add_header(
                    "Content-Type", "application/x-www-form-urlencoded; charset=utf-8"
                )
            with urllib.request.urlopen(req, timeout=REQ_TIMEOUT) as resp:
                return resp.read()
        except urllib.error.HTTPError as e:
            body = b""
            try:
                body = e.read()
            except Exception:
                pass
            last = f"HTTP {e.code}: {body[:300].decode('utf-8', 'replace')}"
            # non-retryable client errors (keep 429 retryable)
            if 400 <= e.code < 500 and e.code != 429:
                raise RuntimeError(f"{url} -> {last}")
        except Exception as e:  # noqa: BLE001 - timeouts, dns, connection reset
            last = f"{type(e).__name__}: {e}"
        wait = min(60, 5 * (2**attempt))
        log(f"  retry {attempt + 1}/5 in {wait}s ({last})")
        time.sleep(wait)
    raise RuntimeError(f"request failed after retries: {url} -> {last}")


def submit(title, fasta_text):
    body = urllib.parse.urlencode(
        {
            "email": EMAIL,
            "title": title,
            "stype": "p",
            "goterms": "true",
            "pathways": "true",
            "sequence": fasta_text,
        }
    ).encode()
    raw = http(f"{BASE}/run", data=body, method="POST")
    text = raw.decode("utf-8", "replace").strip()
    if text.startswith("<?xml") or "<error" in text:
        raise RuntimeError(f"submit rejected: {text[:300]}")
    job_id = text
    if not job_id.startswith("iprscan5-"):
        raise RuntimeError(f"unexpected submit response: {text[:300]}")
    return job_id


def wait_for(job_id, deadline):
    interval = 90
    last = None
    while time.time() < deadline:
        try:
            st = http(f"{BASE}/status/{job_id}").decode().strip()
        except Exception as e:  # noqa: BLE001
            log(f"  status check failed ({e}); retrying")
            st = last or "RUNNING"
        if st != last:
            log(f"  job {job_id}: {st}")
            last = st
        if st in ("FINISHED", "FAILURE", "ERROR", "NOT_FOUND"):
            return st
        time.sleep(interval)
    log(f"  job {job_id}: TIMEOUT after {JOB_TIMEOUT}s")
    return "TIMEOUT"


def main():
    seqs = read_fasta(FASTA)
    log(f"loaded {len(seqs)} sequences")

    state = {}
    if os.path.exists(STATE):
        state = json.load(open(STATE))
        log(f"resumed state with {len(state)} jobs")

    # ---- submit phase -------------------------------------------------
    for i in range(0, len(seqs), CHUNK):
        chunk = seqs[i : i + CHUNK]
        idx = i // CHUNK + 1
        key = f"job{idx}"
        if state.get(key, {}).get("job_id"):
            log(f"{key}: already submitted ({state[key]['job_id']}), skipping")
            continue
        fasta_text = "".join(f">{n}\n{s}\n" for n, s in chunk)
        log(f"{key}: submitting {len(chunk)} sequences ...")
        try:
            jid = submit(f"barley_drought_top100_part{idx}", fasta_text)
            state[key] = {
                "job_id": jid,
                "n_seqs": len(chunk),
                "ids": [n for n, _ in chunk],
                "submitted_at": time.time(),
                "status": "SUBMITTED",
                "attempt": 1,
            }
            log(f"{key}: submitted -> {jid}")
        except Exception as e:  # noqa: BLE001
            log(f"{key}: SUBMIT FAILED: {e}")
            state[key] = {
                "job_id": None,
                "n_seqs": len(chunk),
                "ids": [n for n, _ in chunk],
                "error": str(e),
            }
        json.dump(state, open(STATE, "w"), indent=2)
        time.sleep(3)

    # ---- poll phase ---------------------------------------------------
    for key, info in sorted(state.items()):
        jid = info.get("job_id")
        if not jid or info.get("status") == "FINISHED":
            continue
        st = wait_for(jid, time.time() + JOB_TIMEOUT)
        # one resubmission attempt on failure
        if st in ("FAILURE", "ERROR", "NOT_FOUND") and info.get("attempt", 1) < 2:
            log(f"{key}: got {st}, resubmitting once ...")
            wanted = set(info["ids"])
            sub = [s for s in seqs if s[0] in wanted]
            fasta_text = "".join(f">{n}\n{s}\n" for n, s in sub)
            try:
                newid = submit(f"barley_drought_top100_part{key[-1]}_retry", fasta_text)
                info["job_id"] = newid
                info["attempt"] = 2
                json.dump(state, open(STATE, "w"), indent=2)
                log(f"{key}: resubmitted -> {newid}")
                st = wait_for(newid, time.time() + JOB_TIMEOUT)
            except Exception as e:  # noqa: BLE001
                log(f"{key}: resubmit failed: {e}")
        info["status"] = st
        info["finished_at"] = time.time()
        json.dump(state, open(STATE, "w"), indent=2)

    # ---- download phase -----------------------------------------------
    for key, info in sorted(state.items()):
        jid = info.get("job_id")
        if not jid or info.get("status") != "FINISHED":
            log(f"{key}: skip download (status={info.get('status')})")
            continue
        dest = os.path.join(OUTDIR, f"iprscan_{key}.tsv")
        if os.path.exists(dest) and os.path.getsize(dest) > 0:
            log(f"{key}: {dest} exists, skipping")
            continue
        try:
            raw = http(f"{BASE}/result/{jid}/tsv")
            open(dest, "wb").write(raw)
            log(f"{key}: saved {dest} ({len(raw)} bytes)")
        except Exception as e:  # noqa: BLE001
            log(f"{key}: TSV download failed ({e}), trying XML")
            try:
                raw = http(f"{BASE}/result/{jid}/xml")
                open(os.path.join(OUTDIR, f"iprscan_{key}.xml"), "wb").write(raw)
                log(f"{key}: saved XML fallback")
            except Exception as e2:  # noqa: BLE001
                log(f"{key}: XML fallback also failed: {e2}")

    # ---- merge phase --------------------------------------------------
    # InterProScan's REST TSV ships without a header row; add the documented one.
    IPR_HEADER = [
        "protein_accession",
        "sequence_md5",
        "sequence_length",
        "analysis",
        "signature_accession",
        "signature_description",
        "start",
        "stop",
        "score",
        "status",
        "date",
        "interpro_accession",
        "interpro_description",
        "go_annotations",
        "pathways",
    ]
    merged = os.path.join(OUTDIR, "iprscan_top100.tsv")
    header = None
    n_rows = 0
    with open(merged, "w") as out:
        out.write("\t".join(IPR_HEADER) + "\n")
        for key in sorted(state):
            part = os.path.join(OUTDIR, f"iprscan_{key}.tsv")
            if not os.path.exists(part) or os.path.getsize(part) == 0:
                continue
            with open(part) as fh:
                lines = [ln.rstrip("\n") for ln in fh if ln.strip()]
            if not lines:
                continue
            if lines[0].lower().startswith("protein accession"):
                if header is None:
                    header = lines[0]
                    out.write(header + "\n")
                lines = lines[1:]
            for ln in lines:
                out.write(ln + "\n")
                n_rows += 1
    log(f"merged -> {merged} ({n_rows} data rows)")
    json.dump(state, open(STATE, "w"), indent=2)
    log("DONE")


if __name__ == "__main__":
    sys.exit(main())
