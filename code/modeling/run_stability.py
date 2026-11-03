#!/usr/bin/env python3
"""Stage 5d: seed-stability check.

Re-runs the LeaveOneBlockOut protocol with three different random seeds, applied both to the
k-means block partition and to the models, keeping the hyper-parameters fixed at the values the
main nested run selected most often (so the only thing changing is the seed).
"""
import argparse
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modeling_lib import (FEATURE_SETS, OUTDIR, best_f1_threshold, load_points,
                          make_estimator, make_targets, metrics_at, prune_collinear,
                          spatial_blocks)

SEEDS = [42, 7, 2024]
SPECS = [("T1_qingke_global", "A"), ("T2_highalt_global", "A"), ("T3_qingke_plateau", "B")]
MIN_POS = 5


def modal_params(cv, tname, fs, model, fallback):
    sub = cv[(cv.target == tname) & (cv.feature_set == fs) & (cv.model == model) & cv.valid]
    vals = [v for v in sub["best_params"].tolist() if isinstance(v, str) and v.startswith("{")]
    if not vals:
        return fallback
    return eval(Counter(vals).most_common(1)[0][0])


def default_params(model, pos_weight):
    return {"LR": {"C": 0.3},
            "RF": {"max_depth": None, "min_samples_leaf": 1},
            "XGB": {"max_depth": 4, "scale_pos_weight": pos_weight}}[model]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", default=",".join(map(str, SEEDS)))
    ap.add_argument("--tag", default="")
    args = ap.parse_args()
    seeds = [int(s) for s in args.seeds.split(",")]

    df = load_points()
    targets = make_targets(df)
    cv = pd.read_csv(os.path.join(OUTDIR, f"cv_metrics{args.tag}.tsv"), sep="\t")

    rows = []
    for tname, fs in SPECS:
        tinfo = targets[tname]
        sub = df[tinfo["mask"]].reset_index(drop=True)
        y = tinfo["y"]
        X = sub[FEATURE_SETS[fs]].copy()
        for name in ("LR", "RF", "XGB"):
            try:
                import xgboost  # noqa: F401
            except ImportError:
                break
            for seed in seeds:
                blocks = spatial_blocks(sub[["lon", "lat"]].values, k=10, seed=seed)
                aucs = []
                for b in np.unique(blocks):
                    te = np.where(blocks == b)[0]
                    tr = np.where(blocks != b)[0]
                    if y[tr].sum() < MIN_POS or y[te].sum() < MIN_POS:
                        continue
                    if y[te].sum() in (0, len(te)):
                        continue
                    pos_weight = float((y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1))
                    params = modal_params(cv, tname, fs, name, default_params(name, pos_weight))
                    if name == "XGB":
                        params = dict(params, scale_pos_weight=pos_weight)
                    Xtr, Xte = X.iloc[tr], X.iloc[te]
                    if name == "LR":
                        kept, _ = prune_collinear(Xtr)
                        Xtr, Xte = Xtr[kept], Xte[kept]
                    est = make_estimator(name, params, seed=seed)
                    est.fit(Xtr, y[tr])
                    p_tr = est.predict_proba(Xtr)[:, 1]
                    thr, _ = best_f1_threshold(y[tr], p_tr)
                    p = est.predict_proba(Xte)[:, 1]
                    m = metrics_at(y[te], p, thr)
                    aucs.append(m["auc"])
                    rows.append(dict(target=tname, feature_set=fs, model=name, seed=seed,
                                     fold=int(b), n_pos_test=int(y[te].sum()),
                                     auc=m["auc"], pr_auc=m["pr_auc"], f1=m["f1"],
                                     bacc=m["bacc"], params=str(params)))
                print(f"[seed {seed}] {tname}/{fs}/{name}: mean AUC="
                      f"{np.nanmean(aucs) if aucs else float('nan'):.4f} over {len(aucs)} folds",
                      flush=True)

    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(OUTDIR, f"stability_seeds{args.tag}.tsv"), sep="\t", index=False)
    summ = (out.groupby(["target", "feature_set", "model", "seed"])
            .agg(n_folds=("auc", "size"), auc_mean=("auc", "mean"), auc_sd=("auc", "std"),
                 pr_auc_mean=("pr_auc", "mean"), f1_mean=("f1", "mean"))
            .reset_index())
    across = (summ.groupby(["target", "feature_set", "model"])
              .agg(seed_mean_auc_min=("auc_mean", "min"), seed_mean_auc_max=("auc_mean", "max"),
                   seed_mean_auc_sd=("auc_mean", "std"),
                   mean_within_seed_sd=("auc_sd", "mean")).reset_index())
    across.to_csv(os.path.join(OUTDIR, f"stability_summary{args.tag}.tsv"), sep="\t", index=False)
    print("\n", across.round(4).to_string(index=False))
    print("\n[done] stability stage")


if __name__ == "__main__":
    main()
