#!/usr/bin/env python3
"""Merge InterProScan + eggNOG-substitute annotation into the top100 drought table.

Skeleton: top100_annotation.tsv  (left join, original rows and column order preserved)
Added columns: interpro_entries, go_terms, cog_categories, kegg, eggnog_desc, pfam_from_eggnog

Protein-ID handling: MorexV3 protein IDs may carry a `.CDS1`/`.CDS3` suffix
(e.g. HORVU.MOREX.r3.7HG0637810.1.CDS1).  We try an exact match first and fall back
to the ID with the trailing `.CDS<n>` removed.

Output: top100_annotation_full.tsv
"""
import collections
import csv
import os
import re
import sys

OUTDIR = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/annotation')
SKELETON = os.path.join(OUTDIR, "top100_annotation.tsv")
FASTA = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/top100_proteins.fasta')
IPR = os.path.join(OUTDIR, "iprscan_top100.tsv")
UP_FALLBACK = os.path.join(OUTDIR, "top100_uniprot_eggnog_fallback.tsv")
OUT = os.path.join(OUTDIR, "top100_annotation_full.tsv")

NEW_COLS = [
    "interpro_entries",
    "go_terms",
    "cog_categories",
    "kegg",
    "eggnog_desc",
    "pfam_from_eggnog",
]

# COG functional-category names -> one-letter code (NCBI COG standard)
COG_NAME2CODE = {
    "translation, ribosomal structure and biogenesis": "J",
    "rna processing and modification": "A",
    "transcription": "K",
    "replication, recombination and repair": "L",
    "chromatin structure and dynamics": "B",
    "cell cycle control, cell division, chromosome partitioning": "D",
    "nuclear structure": "Y",
    "defense mechanisms": "V",
    "signal transduction mechanisms": "T",
    "cell wall/membrane/envelope biogenesis": "M",
    "cell motility": "N",
    "cytoskeleton": "Z",
    "extracellular structures": "W",
    "intracellular trafficking, secretion, and vesicular transport": "U",
    "posttranslational modification, protein turnover, chaperones": "O",
    "energy production and conversion": "C",
    "carbohydrate transport and metabolism": "G",
    "amino acid transport and metabolism": "E",
    "nucleotide transport and metabolism": "F",
    "coenzyme transport and metabolism": "H",
    "lipid transport and metabolism": "I",
    "inorganic ion transport and metabolism": "P",
    "secondary metabolites biosynthesis, transport and catabolism": "Q",
    "general function prediction only": "R",
    "function unknown": "S",
}


# CD-Search short name -> strings that appear in the InterPro entry name for the
# same domain family (different naming vocabularies, same biology).
CDS_ALIASES = {
    "lea_1": ["lea_1", "lea"],
    "nam": ["nac domain", "nam"],
    "lrrnt_2": ["leucine-rich repeat-containing n-terminal", "lrrnt"],
    "a_thal_3526": ["chp01589", "a_thal_3526"],
}


def norm_prot(pid):
    """Strip a trailing .CDS<n> suffix."""
    return re.sub(r"\.CDS\d+$", "", pid.strip())


def gene_of(protein_id):
    """HORVU.MOREX.r3.1HG0002450.1[.CDS3] -> HORVU.MOREX.r3.1HG0002450"""
    return ".".join(protein_id.split(".")[:4])


def read_fasta_ids(path):
    ids = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                ids.append(line[1:].strip().split()[0])
    return ids


def read_tsv(path):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def split_multi(v):
    return [x.strip() for x in re.split(r"[;,]", v or "") if x.strip()]


def clean(v):
    """InterProScan writes '-' for 'no value'; normalise it away."""
    v = (v or "").strip()
    return "" if v in ("-", "NULL", "None") else v


def main():
    skeleton = read_tsv(SKELETON)
    fieldnames = list(skeleton[0].keys())
    print(f"skeleton: {len(skeleton)} rows, {len(fieldnames)} columns")

    # ---- gene -> ordered protein ids (fasta order) --------------------
    prot_ids = read_fasta_ids(FASTA)
    gene2prots = collections.OrderedDict()
    for p in prot_ids:
        gene2prots.setdefault(gene_of(p), []).append(p)
    print(f"fasta: {len(prot_ids)} proteins in {len(gene2prots)} genes")

    def lookup(container, pid):
        """Exact match first, then match after stripping .CDS<n>."""
        if pid in container:
            return container[pid]
        n = norm_prot(pid)
        for key in (n,):
            if key in container:
                return container[key]
        # last resort: any key that normalises to the same id
        for k, v in container.items():
            if norm_prot(k) == n:
                return v
        return None

    # ---- InterProScan -------------------------------------------------
    ipr_by_prot = collections.defaultdict(list)
    ipr_rows = read_tsv(IPR)
    print(f"interproscan: {len(ipr_rows)} rows for {len(set(r['protein_accession'] for r in ipr_rows))} proteins")
    for r in ipr_rows:
        ipr_by_prot[r["protein_accession"]].append(r)

    # global IPR accession -> entry name, so UniProt-supplied accessions can be named too
    ipr_acc2name = {}
    for r in ipr_rows:
        acc, name = clean(r.get("interpro_accession")), clean(r.get("interpro_description"))
        if acc.startswith("IPR") and name:
            ipr_acc2name.setdefault(acc, name)

    # ---- UniProt / eggNOG substitute ----------------------------------
    up_rows = {r["Entry"]: r for r in read_tsv(UP_FALLBACK)}
    print(f"uniprot fallback: {len(up_rows)} entries")

    # protein -> accession map produced by the fetcher (blast-level mapping)
    p2a = {}
    map_path = os.path.join(OUTDIR, "protein_to_uniprot_ac.tsv")
    if os.path.exists(map_path):
        for r in read_tsv(map_path):
            p2a[r["protein_id"]] = r["uniprot_ac"]

    # ---- CD-Search COG categories (from the local hit table) ----------
    cog_by_gene = collections.defaultdict(set)
    cds_hits = os.path.join(OUTDIR, "top100_cdsearch_hits.tsv")
    if os.path.exists(cds_hits):
        for r in read_tsv(cds_hits):
            m = re.search(r">([^\[\]]+)", r.get("Query", ""))
            if not m:
                continue
            g = gene_of(m.group(1).strip())
            for name in re.findall(r"\[([^\]]+)\]", r.get("Definition", "") or ""):
                code = COG_NAME2CODE.get(name.strip().lower())
                if code:
                    cog_by_gene[g].add(code)

    # ---- build merged rows --------------------------------------------
    stats = collections.Counter()
    go_counter = collections.Counter()
    consistency = []  # (gene, cds_domain, interpro_entries, matched)
    out_rows = []

    for row in skeleton:
        gene = row["gene_id"]
        prots = gene2prots.get(gene, [])
        # primary protein for this gene = first in fasta order
        primary = prots[0] if prots else None

        # ------ InterProScan annotation for the primary protein --------
        ipr_hits = []
        if primary:
            ipr_hits = lookup(ipr_by_prot, primary) or []
        if not ipr_hits:
            # fall back to any protein of this gene
            for p in prots:
                h = lookup(ipr_by_prot, p)
                if h:
                    ipr_hits = h
                    break

        ipr_names, ipr_accs, go_terms, pfam_is = [], [], [], []
        for h in ipr_hits:
            acc = clean(h.get("interpro_accession"))
            name = clean(h.get("interpro_description"))
            if acc.startswith("IPR") and acc not in ipr_accs:
                ipr_accs.append(acc)
                ipr_names.append(name or acc)
            if h.get("analysis") == "Pfam":
                sig = clean(h.get("signature_accession"))
                if sig and sig not in pfam_is:
                    pfam_is.append(sig)
            for g in re.findall(r"GO:\d+", clean(h.get("go_annotations"))):
                if g not in go_terms:
                    go_terms.append(g)

        # ------ UniProt-derived (eggNOG / KEGG / Pfam / GO) ------------
        up = None
        for cand in ([primary] if primary else []) + prots:
            if not cand:
                continue
            ac = p2a.get(cand) or p2a.get(norm_prot(cand))
            if ac and ac in up_rows:
                up = up_rows[ac]
                break
        # also honour the accession already present in the skeleton
        if up is None:
            for ac in split_multi(row.get("uniprot_ac", "")):
                if ac in up_rows:
                    up = up_rows[ac]
                    break

        up_pfam, up_go, up_ipr, up_eggnog, up_kegg, up_desc = [], [], [], [], [], ""
        if up:
            up_pfam = split_multi(up.get("Pfam", ""))
            up_go = re.findall(r"GO:\d+", up.get("Gene Ontology IDs", "") or "")
            up_ipr = split_multi(up.get("InterPro", ""))
            up_eggnog = split_multi(up.get("eggNOG", ""))
            up_kegg = split_multi(up.get("KEGG", ""))
            up_desc = (up.get("Protein names") or "").strip()

        # suppress IPR ids that UniProt lists as non-InterPro placeholders
        up_ipr = [i for i in up_ipr if i.startswith("IPR")]

        # InterPro entry names: InterProScan descriptions, else UniProt accessions
        # (UniProt accessions get a global name from the InterProScan set when known)
        if not ipr_names and up_ipr:
            for a in up_ipr:
                nm = ipr_acc2name.get(a)
                ipr_names.append(f"{a}:{nm}" if nm else a)

        # GO terms: keep InterProScan terms, add UniProt ones not already present
        if not go_terms:
            go_terms = list(up_go)
        else:
            for g in up_go:
                if g not in go_terms:
                    go_terms.append(g)

        # Pfam: InterProScan PF accessions + UniProt Pfam accessions
        pfams = list(pfam_is)
        for p in up_pfam:
            if re.match(r"PF\d+", p) and p not in pfams:
                pfams.append(p)

        # COG categories + eggNOG ortholog group ids
        parts = []
        letters = sorted(cog_by_gene.get(gene, set()))
        if letters:
            parts.append("COG_cat:" + ";".join(letters))
        if up_eggnog:
            parts.append("eggNOG_OGs:" + ";".join(up_eggnog))
        cog_categories = "|".join(parts)

        new = {
            "interpro_entries": ";".join(ipr_names),
            "go_terms": ";".join(go_terms),
            "cog_categories": cog_categories,
            "kegg": ";".join(up_kegg),
            "eggnog_desc": up_desc,
            "pfam_from_eggnog": ";".join(pfams),
        }
        out_rows.append({**row, **new})

        # ---- stats ----------------------------------------------------
        if new["interpro_entries"]:
            stats["genes_with_interpro"] += 1
        if new["pfam_from_eggnog"]:
            stats["genes_with_pfam"] += 1
        if new["go_terms"]:
            stats["genes_with_go"] += 1
            for g in go_terms:
                go_counter[g] += 1
        if new["kegg"]:
            stats["genes_with_kegg"] += 1
        if new["cog_categories"]:
            stats["genes_with_cog_or_eggnog"] += 1
        if new["eggnog_desc"]:
            stats["genes_with_desc"] += 1
        if ipr_hits:
            stats["genes_interproscan_hit"] += 1

        # ---- CD-Search vs InterPro consistency (keyed on the CD-Search domain)
        dom = clean(row.get("domains", ""))
        if dom:
            dom_key = dom.split("(")[0].strip()
            blob = " ".join(ipr_names).lower()
            # CD-Search short names vs InterPro entry names are different vocabularies;
            # accept the documented synonyms so the check reflects biology, not naming.
            needles = CDS_ALIASES.get(dom_key.lower(), [dom_key.lower()])
            matched = any(n in blob for n in needles)
            stats["cdsearch_domains"] += 1
            if matched:
                stats["cdsearch_confirmed_by_interpro"] += 1
            consistency.append((gene, dom_key, new["interpro_entries"], matched))

    # ---- write output --------------------------------------------------
    with open(OUT, "w", newline="") as out:
        w = csv.DictWriter(out, fieldnames=fieldnames + NEW_COLS, delimiter="\t")
        w.writeheader()
        for r in out_rows:
            w.writerow(r)
    print(f"\nwrote {OUT} ({len(out_rows)} rows, {len(fieldnames) + len(NEW_COLS)} columns)")

    # ---- report statistics --------------------------------------------
    n = len(out_rows)
    print("\n===== coverage (n=%d genes) =====" % n)
    for k in [
        "genes_interproscan_hit",
        "genes_with_interpro",
        "genes_with_pfam",
        "genes_with_go",
        "genes_with_kegg",
        "genes_with_cog_or_eggnog",
        "genes_with_desc",
    ]:
        print(f"  {k:32s} {stats[k]:3d}  ({100.0 * stats[k] / n:.0f}%)")

    print("\n===== top 10 GO terms by gene count =====")
    for g, c in go_counter.most_common(10):
        print(f"  {g}  {c}")

    print("\n===== CD-Search vs InterProScan consistency =====")
    print(
        f"  CD-Search domain calls: {stats['cdsearch_domains']}; "
        f"also found by InterPro: {stats['cdsearch_confirmed_by_interpro']}"
    )
    for gene, dom, ipr, matched in consistency:
        flag = "OK " if matched else "MISS"
        print(f"  [{flag}] {gene:32s} CD-Search={dom:14s} InterPro={ipr or '-'}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
