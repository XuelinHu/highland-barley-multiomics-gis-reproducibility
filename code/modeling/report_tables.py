#!/usr/bin/env python3
"""Print the report tables as markdown, straight from the result TSVs (no manual transcription)."""
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modeling_lib import OUTDIR

TAG = sys.argv[1] if len(sys.argv) > 1 else ""


def md(df, floatfmt=4):
    df = df.copy()
    for c in df.columns:
        if df[c].dtype.kind == "f":
            df[c] = df[c].map(lambda v: "" if pd.isna(v) else f"{v:.{floatfmt}f}")
    head = "| " + " | ".join(map(str, df.columns)) + " |"
    sep = "|" + "|".join(["---"] * len(df.columns)) + "|"
    body = ["| " + " | ".join(map(str, r)) + " |" for r in df.itertuples(index=False)]
    return "\n".join([head, sep] + body)


def main():
    cv = pd.read_csv(os.path.join(OUTDIR, f"cv_metrics{TAG}.tsv"), sep="\t")
    summ = pd.read_csv(os.path.join(OUTDIR, f"metrics_summary{TAG}.tsv"), sep="\t")
    dist = pd.read_csv(os.path.join(OUTDIR, f"block_positive_distribution{TAG}.tsv"), sep="\t")

    print("### TABLE block distribution")
    d = dist[dist["n_pos"] > 0].copy()
    print(md(d))
    print()
    print("### TABLE per-block positive summary")
    g = (dist.groupby("target")
         .agg(n_blocks=("block", "size"), blocks_with_pos=("n_pos", lambda s: (s > 0).sum()),
              blocks_ge5=("n_pos", lambda s: (s >= 5).sum()), n_pos=("n_pos", "sum"),
              max_block_share=("pos_share_pct", "max"))
         .reset_index())
    print(md(g))
    print()
    print("### TABLE metrics summary (full feature sets)")
    s = summ[summ["feature_set"].isin(["A", "B"])].copy()
    print(md(s[["target", "feature_set", "model", "n_folds", "n_folds_skipped", "auc_mean",
                "auc_sd", "pr_auc_mean", "pr_auc_sd", "f1_mean", "f1_sd", "bacc_mean",
                "bacc_sd"]]))
    print()
    print("### TABLE metrics summary (no-elev ablations)")
    s = summ[~summ["feature_set"].isin(["A", "B"])].copy()
    print(md(s[["target", "feature_set", "model", "n_folds", "auc_mean", "auc_sd",
                "pr_auc_mean", "f1_mean", "bacc_mean"]]))
    print()
    print("### TABLE per-fold detail")
    print(md(cv[["target", "feature_set", "model", "fold", "valid", "n_pos_train",
                 "n_pos_test", "auc", "pr_auc", "f1", "bacc", "threshold",
                 "threshold_source"]]))

    fi = os.path.join(OUTDIR, f"feature_importance{TAG}.tsv")
    if os.path.exists(fi):
        imp = pd.read_csv(fi, sep="\t")
        print()
        for (t, f, m), sub in imp.groupby(["target", "feature_set", "model"]):
            print(f"### TABLE top-12 SHAP {t} | set {f} | {m}")
            print(md(sub.sort_values("shap_rank").head(12)[
                ["shap_rank", "feature", "mean_abs_shap", "rf_gini_importance"]]))
            print()

    st = os.path.join(OUTDIR, f"stability_summary{TAG}.tsv")
    if os.path.exists(st):
        print("### TABLE seed stability")
        print(md(pd.read_csv(st, sep="\t")))


if __name__ == "__main__":
    main()
