#!/usr/bin/env python3
"""Stage 5a: spatial block cross-validation for the barley environmental models.

Usage
  python run_cv.py                # full run (T1/T2/T3, LR/RF/XGB)
  python run_cv.py --targets T3_qingke_plateau --models LR   # smoke test
"""
import argparse
import os
import sys
import time
import traceback

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modeling_lib import (FEATURE_SETS, MODEL_NAMES, OUTDIR, best_f1_threshold,
                          block_distribution, load_points, make_estimator,
                          make_targets, metrics_at, param_grid, prune_collinear,
                          spatial_blocks)

from sklearn.model_selection import GroupKFold

MIN_POS_TRAIN = 5     # a fold is only scored if the training part can support a model
MIN_POS_TEST = 5      # and if the held-out block contains enough positives to score it
INNER_SPLITS = 3
SCALE_K = {"G": 10, "P": 10}


def inner_tune(name, X, y, groups, seed, pos_weight, log):
    """Search the grid with an inner GroupKFold over training blocks.

    Returns (best_params, inner_oof_proba, per_config_auc).
    """
    grid = param_grid(name, pos_weight)
    uniq = np.unique(groups)
    n_splits = min(INNER_SPLITS, len(uniq))
    results, oof_store = [], {}
    if n_splits >= 2:
        gkf = GroupKFold(n_splits=n_splits)
        splits = list(gkf.split(X, y, groups))
    else:
        splits = []
    for params in grid:
        oof = np.full(len(y), np.nan)
        aucs = []
        for tr, te in splits:
            if len(np.unique(y[tr])) < 2 or len(np.unique(y[te])) < 2:
                continue
            est = make_estimator(name, params, seed=seed)
            est.fit(X.iloc[tr], y[tr])
            p = est.predict_proba(X.iloc[te])[:, 1]
            oof[te] = p
            aucs.append(metrics_at(y[te], p, 0.5)["auc"])
        mean_auc = float(np.nanmean(aucs)) if aucs else np.nan
        results.append((mean_auc, params))
        oof_store[id(params)] = oof
    valid = [r for r in results if np.isfinite(r[0])]
    if not valid:
        best = grid[0]
        return best, np.full(len(y), np.nan), results
    best_auc, best = max(valid, key=lambda r: (r[0], -len(str(r[1]))))
    log(f"      inner: " + ", ".join(
        f"{p}->{a:.4f}" if np.isfinite(a) else f"{p}->NA" for a, p in results) +
        f" | best {best} ({best_auc:.4f})")
    return best, oof_store[id(best)], results


def run_fold(name, fs_name, X, y, tr, te, groups_tr, seed, log):
    """One LeaveOneBlockOut fold: prune (LR) -> tune -> threshold -> score."""
    pos_weight = float((y[tr] == 0).sum() / max((y[tr] == 1).sum(), 1))
    Xtr, Xte = X.iloc[tr], X.iloc[te]
    dropped = []
    if name == "LR":
        kept, dropped = prune_collinear(Xtr)
        if not kept:
            kept = list(Xtr.columns)
            dropped = []
        Xtr, Xte = Xtr[kept], Xte[kept]
    best, oof, grid_res = inner_tune(name, Xtr, y[tr], groups_tr, seed, pos_weight, log)
    thr_inner = np.nan
    if np.isfinite(oof).any():
        thr_inner, _ = best_f1_threshold(y[tr][np.isfinite(oof)], oof[np.isfinite(oof)])
    est = make_estimator(name, best, seed=seed)
    est.fit(Xtr, y[tr])
    # primary rule (as designed): F1-optimal threshold on the training fold itself.  When the
    # inner CV could produce out-of-fold probabilities those are used instead (less optimistic).
    thr_train, _ = best_f1_threshold(y[tr], est.predict_proba(Xtr)[:, 1])
    if np.isfinite(thr_inner):
        thr, source = float(thr_inner), "inner_oof"
    else:
        thr, source = float(thr_train), "train_insample"
    proba = est.predict_proba(Xte)[:, 1]
    m = metrics_at(y[te], proba, thr)
    return dict(best_params=str(best), threshold=float(thr), threshold_source=source,
                threshold_train_insample=float(thr_train), threshold_inner_oof=float(thr_inner),
                n_features=int(Xtr.shape[1]),
                n_dropped=len(dropped), dropped=";".join(dropped), **m)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--targets", default="all")
    ap.add_argument("--models", default=",".join(MODEL_NAMES))
    ap.add_argument("--seed", type=int, default=42)
    ap.add_argument("--blocks", type=int, default=SCALE_K["G"])
    ap.add_argument("--tag", default="")
    ap.add_argument("--no-final", action="store_true", help="skip saving full-data models")
    args = ap.parse_args()

    os.makedirs(OUTDIR, exist_ok=True)
    logf = open(os.path.join(OUTDIR, f"run_cv{args.tag}.log"), "a")

    def log(msg):
        print(msg, flush=True)
        logf.write(str(msg) + "\n")
        logf.flush()

    t_start = time.time()
    log(f"\n===== run_cv {time.strftime('%Y-%m-%d %H:%M:%S')} args={vars(args)} =====")
    df = load_points()
    targets = make_targets(df)
    models = [m for m in args.models.split(",") if m in MODEL_NAMES]
    if "XGB" in models:
        try:
            import xgboost  # noqa: F401
        except ImportError:
            log("[warn] xgboost not installed -> XGB skipped")
            models = [m for m in models if m != "XGB"]
    want = list(targets) if args.targets == "all" else args.targets.split(",")

    rows, dist_rows, drop_rows, final_rows = [], [], [], []
    for tname in want:
        tinfo = targets[tname]
        sub = df[tinfo["mask"]].reset_index(drop=True)
        y = tinfo["y"]
        k = args.blocks if tinfo["scale"] == "G" else SCALE_K["P"]
        blocks = spatial_blocks(sub[["lon", "lat"]].values, k=k, seed=args.seed)
        dist, dsum = block_distribution(y, blocks)
        log(f"\n[{tname}] n={tinfo['n']} pos={tinfo['n_pos']} scale={tinfo['scale']} k={k}")
        log(f"  blocks: {dsum}")
        dist.insert(0, "target", tname)
        dist_rows.append(dist)
        np.save(os.path.join(OUTDIR, f"blocks_{tname}{args.tag}.npy"), blocks)

        for fs_name in tinfo["feature_sets"]:
            X = sub[FEATURE_SETS[fs_name]].copy()
            for name in models:
                for b in np.unique(blocks):
                    te = np.where(blocks == b)[0]
                    tr = np.where(blocks != b)[0]
                    n_pos_tr, n_pos_te = int(y[tr].sum()), int(y[te].sum())
                    row = dict(target=tname, feature_set=fs_name, model=name, fold=int(b),
                               n_train=len(tr), n_test=len(te),
                               n_pos_train=n_pos_tr, n_pos_test=n_pos_te, seed=args.seed)
                    if n_pos_tr < MIN_POS_TRAIN or n_pos_te < MIN_POS_TEST:
                        row.update(valid=False, skip_reason="too_few_positives",
                                   auc=np.nan, pr_auc=np.nan, f1=np.nan,
                                   bacc=np.nan, threshold=np.nan, threshold_source="",
                                   threshold_train_insample=np.nan, threshold_inner_oof=np.nan,
                                   best_params="", n_features=X.shape[1], n_dropped=0, dropped="")
                        rows.append(row)
                        continue
                    if n_pos_te in (0, len(te)):
                        # held-out block carries a single class -> AUC undefined
                        row.update(valid=False, skip_reason="single_class_block",
                                   auc=np.nan, pr_auc=np.nan, f1=np.nan, bacc=np.nan,
                                   threshold=np.nan, threshold_source="",
                                   threshold_train_insample=np.nan, threshold_inner_oof=np.nan,
                                   best_params="", n_features=X.shape[1],
                                   n_dropped=0, dropped="")
                        rows.append(row)
                        continue
                    try:
                        res = run_fold(name, fs_name, X, y, tr, te, blocks[tr], args.seed,
                                       log=lambda s, n=tname, f=fs_name, m=name, b=b:
                                       log(f"    {n}/{f}/{m}/block{b}: {s}"))
                        row.update(valid=True, skip_reason="", **res)
                        log(f"  {tname}/{fs_name}/{name} block{b}: "
                            f"AUC={res['auc']:.4f} PR={res['pr_auc']:.4f} "
                            f"F1={res['f1']:.4f} BACC={res['bacc']:.4f} "
                            f"thr={res['threshold']:.3f} nfeat={res['n_features']} "
                            f"({time.time()-t_start:.0f}s)")
                        if res["dropped"]:
                            drop_rows.append(dict(target=tname, feature_set=fs_name,
                                                  fold=int(b), model=name,
                                                  dropped=res["dropped"]))
                    except Exception:
                        log(f"  !! fold failed {tname}/{fs_name}/{name}/block{b}\n"
                            + traceback.format_exc())
                        row.update(valid=False, skip_reason="error", auc=np.nan, pr_auc=np.nan,
                                   f1=np.nan, bacc=np.nan, threshold=np.nan,
                                   threshold_source="", threshold_train_insample=np.nan,
                                   threshold_inner_oof=np.nan, best_params="ERROR",
                                   n_features=X.shape[1], n_dropped=0, dropped="")
                    rows.append(row)
                # end folds
                # ---- final model on the whole analysis set (for SHAP / importance)
                if not args.no_final:
                    folds = [r for r in rows if r["target"] == tname and r["feature_set"] == fs_name
                             and r["model"] == name and r.get("valid") and r.get("best_params")]
                    if folds:
                        params = eval(max(set(r["best_params"] for r in folds),
                                          key=lambda s: sum(1 for r in folds
                                                            if r["best_params"] == s)))
                        Xf = X
                        dropf = []
                        if name == "LR":
                            keptf, dropf = prune_collinear(Xf)
                            Xf = Xf[keptf]
                        est = make_estimator(name, params, seed=args.seed)
                        est.fit(Xf, y)
                        import joblib
                        mdir = os.path.join(OUTDIR, "models")
                        os.makedirs(mdir, exist_ok=True)
                        fn = os.path.join(mdir, f"{tname}__{fs_name}__{name}{args.tag}.joblib")
                        joblib.dump({"model": est, "features": list(Xf.columns),
                                     "params": params, "target": tname,
                                     "feature_set": fs_name, "name": name,
                                     "n_pos": int(y.sum()), "n": int(len(y)),
                                     "dropped": dropf}, fn)
                        final_rows.append(dict(target=tname, feature_set=fs_name, model=name,
                                               params=str(params), n_features=Xf.shape[1],
                                               dropped=";".join(dropf), file=fn))
                        log(f"  [final] {tname}/{fs_name}/{name} -> {os.path.basename(fn)} "
                            f"({len(folds)} fold votes)")

    cv = pd.DataFrame(rows)
    cvf = os.path.join(OUTDIR, f"cv_metrics{args.tag}.tsv")
    cv.to_csv(cvf, sep="\t", index=False)
    pd.concat(dist_rows, ignore_index=True).to_csv(
        os.path.join(OUTDIR, f"block_positive_distribution{args.tag}.tsv"), sep="\t", index=False)
    if drop_rows:
        pd.DataFrame(drop_rows).to_csv(
            os.path.join(OUTDIR, f"collinearity_dropped{args.tag}.tsv"), sep="\t", index=False)
    if final_rows:
        pd.DataFrame(final_rows).to_csv(
            os.path.join(OUTDIR, f"final_models{args.tag}.tsv"), sep="\t", index=False)

    valid = cv[cv["valid"]].copy()
    if len(valid):
        agg = (valid.groupby(["target", "feature_set", "model"])
               .agg(n_folds=("auc", "size"), n_pos_test_sum=("n_pos_test", "sum"),
                    auc_mean=("auc", "mean"), auc_sd=("auc", "std"),
                    pr_auc_mean=("pr_auc", "mean"), pr_auc_sd=("pr_auc", "std"),
                    f1_mean=("f1", "mean"), f1_sd=("f1", "std"),
                    bacc_mean=("bacc", "mean"), bacc_sd=("bacc", "std"))
               .reset_index())
        # invalid folds are informative in their own right
        inv = (cv[~cv["valid"]].groupby(["target", "feature_set", "model"])
               .agg(n_folds_skipped=("auc", "size"),
                    skip_reasons=("skip_reason", lambda s: ";".join(sorted(set(s)))))
               .reset_index())
        agg = agg.merge(inv, on=["target", "feature_set", "model"], how="outer")
        agg.to_csv(os.path.join(OUTDIR, f"metrics_summary{args.tag}.tsv"), sep="\t", index=False)
        log("\n===== summary =====")
        log(agg.round(4).to_string(index=False))
    log(f"===== run_cv done in {time.time()-t_start:.0f}s =====")
    logf.close()


if __name__ == "__main__":
    main()
