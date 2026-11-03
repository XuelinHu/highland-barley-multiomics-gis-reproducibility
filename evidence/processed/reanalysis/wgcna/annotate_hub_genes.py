#!/usr/bin/env python3
"""Swiss-Prot annotation of WGCNA hub genes (and optionally all network genes).

Reuses the BLAST database already built for the top100 annotation:
  reference/MorexV3.pep.fa   (MorexV3 proteins)
  reference/sprot_db         (makeblastdb of reference/uniprot_sprot.fasta)
  tools/ncbi-blast-2.17.0+/bin/blastp

Usage:
  annotate_hub_genes.py <genes.txt> <out.tsv>
where <genes.txt> is one Morex gene id per line
  (HORVU.MOREX.r3.<chr>HG<nnnnnnn>).
"""
import os
import re
import subprocess
import sys

ROOT = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), '')
REF = os.path.join(ROOT, "reference")
PEP = os.path.join(REF, "MorexV3.pep.fa")
SPROT_DB = os.path.join(REF, "sprot_db")
BLASTP = os.path.join(ROOT, "tools/ncbi-blast-2.17.0+/bin/blastp")
MAPPING = os.path.join(ROOT, "processed/mapping/GPL14927_morex_mapping.tsv")
TMP = "/tmp/wgcna_hub_annot"


def read_fasta(path):
    """Return {seq_id: sequence} for a (optionally gzipped) FASTA."""
    import gzip
    op = gzip.open if path.endswith(".gz") else open
    seqs, cur, buf = {}, None, []
    with op(path, "rt") as fh:
        for line in fh:
            if line.startswith(">"):
                if cur is not None:
                    seqs[cur] = "".join(buf)
                cur = line[1:].split()[0]
                buf = []
            else:
                buf.append(line.strip())
    if cur is not None:
        seqs[cur] = "".join(buf)
    return seqs


def main():
    genes_file, out_tsv = sys.argv[1], sys.argv[2]
    genes = [l.strip() for l in open(genes_file) if l.strip()]
    genes = list(dict.fromkeys(genes))
    print(f"[annot] {len(genes)} unique genes requested", flush=True)

    # ---- gene -> protein id (from the probe->Morex mapping) ----
    gene2prot = {}
    with open(MAPPING) as fh:
        header = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[idx["status"]] != "unique_gene":
                continue
            g = f[idx["gene_ids"]]
            p = f[idx["protein_ids"]]
            if g and p and g not in gene2prot:
                gene2prot[g] = p.split(";")[0]

    os.makedirs(TMP, exist_ok=True)
    seqs = read_fasta(PEP)
    print(f"[annot] MorexV3 protein fasta: {len(seqs)} sequences", flush=True)

    # gene id -> first available protein isoform, built once (avoids an O(n^2) scan)
    gene2seqid = {}
    for sid in seqs:
        gid = re.sub(r"\.\d+(\.CDS\d+)?$", "", sid)
        gene2seqid.setdefault(gid, sid)

    # ---- collect the protein sequences for the requested genes ----
    query_fa = os.path.join(TMP, "query.fa")
    picked = 0
    with open(query_fa, "w") as out:
        for g in genes:
            pid = gene2prot.get(g)
            hit = None
            if pid is not None:
                cand = [pid, re.sub(r"\.CDS\d+$", "", pid), pid + ".CDS1"]
                hit = next((c for c in cand if c in seqs), None)
            if hit is None:
                # last resort: any protein isoform of this gene
                hit = gene2seqid.get(g)
            if hit is None:
                continue
            out.write(f">{g}\n{seqs[hit]}\n")
            picked += 1
    print(f"[annot] {picked}/{len(genes)} genes had a recoverable protein sequence",
          flush=True)

    if picked == 0:
        with open(out_tsv, "w") as out:
            out.write("gene_id\tswissprot_hit\tswissprot_evalue\tidentity_pct\tcoverage_pct\n")
        return

    # ---- blastp vs Swiss-Prot ----
    raw = os.path.join(TMP, "blast_raw.tsv")
    cmd = [BLASTP, "-query", query_fa, "-db", SPROT_DB,
           "-outfmt", "6 qseqid sseqid pident length evalue bitscore stitle",
           "-evalue", "1e-5", "-max_target_seqs", "1",
           "-num_threads", os.environ.get("WGCNA_BLAST_THREADS", "14"),
           "-out", raw]
    print("[annot] running blastp ...", flush=True)
    subprocess.run(cmd, check=True)

    # ---- best hit per gene ----
    best = {}
    qlen = {}
    with open(query_fa) as fh:
        cur = None
        for line in fh:
            if line.startswith(">"):
                cur = line[1:].strip()
                qlen[cur] = 0
            elif cur:
                qlen[cur] += len(line.strip())

    for line in open(raw):
        f = line.rstrip("\n").split("\t")
        if len(f) < 7:
            continue
        q, s, pid, ln, ev, bs, st = f[0], f[1], float(f[2]), int(f[3]), float(f[4]), float(f[5]), f[6]
        if q not in best or ev < best[q][0]:
            best[q] = (ev, s, pid, ln, bs, st)

    with open(out_tsv, "w") as out:
        out.write("gene_id\tswissprot_hit\tswissprot_evalue\tidentity_pct\tcoverage_pct\n")
        for g in genes:
            if g in best:
                ev, s, pid, ln, bs, st = best[g]
                cov = 100.0 * ln / qlen.get(g, ln)
                out.write(f"{g}\t{st}\t{ev:.3g}\t{pid:.1f}\t{cov:.1f}\n")
            else:
                out.write(f"{g}\t\t\t\t\n")
    print(f"[annot] wrote {out_tsv}: {len(best)} genes with a Swiss-Prot hit", flush=True)


if __name__ == "__main__":
    main()
