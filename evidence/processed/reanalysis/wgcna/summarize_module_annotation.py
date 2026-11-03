#!/usr/bin/env python3
"""Summarise Swiss-Prot annotation per WGCNA module.

Inputs
  <module_genes.tsv>      columns: gene, module, ... (one row per gene)
  <annotation.tsv>        output of annotate_hub_genes.py
                          (gene_id, swissprot_hit, swissprot_evalue, ...)
  <hub_genes.tsv>         optional; hub gene table for a focused view

Output
  module_annotation_summary.tsv  one row per module:
    module, n_genes, n_annotated, pct_annotated,
    top_descriptions (semicolon separated, most frequent first),
    n_stress_related, stress_related_genes
  Also prints a per-module block to stdout for the report.
"""
import os
import re
import sys
from collections import Counter, defaultdict

OUT = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'processed/reanalysis/wgcna')

STRESS_KW = [
    "dehydrin", "lea ", "late embryogenesis", "heat shock", "hsp",
    "abscisic", "aba", "drought", "cold", "freezing", "osmotic",
    "salt", "water stress", "desiccation", "aquaporin", "peroxidase",
    "glutathione", "catalase", "superoxide", "thioredoxin", "chaperone",
    "transcription factor", "myb", "nac", "bzip", "dreb", "wrky", "ap2",
    "erf", "cbf", "aquaporin", "pip", "tip", "sucrose", "proline",
    "betaine", "trehalose", "universal stress", "usp", "cold shock",
    "damage", "senescence", "cyp", "cytochrome p450",
]


def short_desc(hit):
    """Reduce a Swiss-Prot hit to a readable protein name.

    Handles both forms:
      'sp|Q39138|LEA6_ARATH Late embryogenesis abundant protein 6 OS=...'
      'Late embryogenesis abundant protein 6 OS=Arabidopsis thaliana OX=...'
    """
    if not hit:
        return ""
    hit = hit.strip()
    acc = ""
    m = re.match(r"sp\|([A-Z0-9_]+)\|([A-Z0-9_]+)\s+(.*)$", hit)
    if m:
        acc, hit = m.group(2), m.group(3)
    # drop the OS=/OX=/GN=/PE=/SV= annotation tail
    name = re.split(r"\s+OS=", hit, maxsplit=1)[0].strip()
    name = re.sub(r"\s*\(Fragment\)$", "", name).strip()
    return f"{acc}: {name}" if acc else name


def is_stress(desc):
    d = desc.lower()
    return any(k in d for k in STRESS_KW)


def main():
    module_genes_tsv, annotation_tsv = sys.argv[1], sys.argv[2]
    hub_tsv = sys.argv[3] if len(sys.argv) > 3 else None

    ann = {}
    with open(annotation_tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(head)}
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if len(f) > idx["swissprot_hit"] and f[idx["swissprot_hit"]]:
                ann[f[idx["gene_id"]]] = f[idx["swissprot_hit"]]

    by_mod = defaultdict(list)
    with open(module_genes_tsv) as fh:
        head = fh.readline().rstrip("\n").split("\t")
        gi, mi = head.index("gene"), head.index("module")
        for line in fh:
            f = line.rstrip("\n").split("\t")
            if f[mi] and f[mi] != "grey":
                by_mod[f[mi]].append(f[gi])

    hub_genes = set()
    if hub_tsv:
        with open(hub_tsv) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            gi = head.index("gene")
            for line in fh:
                hub_genes.add(line.split("\t")[gi])

    rows = []
    for mod in sorted(by_mod):
        genes = by_mod[mod]
        descs = [short_desc(ann[g]) for g in genes if g in ann]
        cnt = Counter(d for d in descs if d)
        stress = sorted({g for g in genes if g in ann and is_stress(ann[g])})
        top = "; ".join(f"{d} ({c})" for d, c in cnt.most_common(8))
        rows.append((mod, len(genes), len(descs), 100.0 * len(descs) / max(1, len(genes)),
                     top, len(stress)))

        print(f"\n### module {mod}  (n={len(genes)}, annotated {len(descs)}, "
              f"stress-related {len(stress)})")
        for d, c in cnt.most_common(10):
            print(f"   {c:5d}  {d}")
        hubs = [g for g in genes if g in hub_genes and g in ann]
        if hubs:
            print("   -- hub genes with annotation:")
            for g in hubs:
                mark = " *" if g in hub_genes else ""
                print(f"      {g}  {short_desc(ann[g])}{mark}")

    out = os.path.join(OUT, "module_annotation_summary.tsv")
    with open(out, "w") as fh:
        fh.write("module\tn_genes\tn_annotated\tpct_annotated\ttop_descriptions\tn_stress_related\n")
        for r in rows:
            fh.write(f"{r[0]}\t{r[1]}\t{r[2]}\t{r[3]:.1f}\t{r[4]}\t{r[5]}\n")
    print(f"\n[summary] wrote {out}")


if __name__ == "__main__":
    main()
