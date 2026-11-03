#!/usr/bin/env python3
"""Stage 5c: TreeSHAP explanation of the fitted RF / XGB models.

Reads the models saved by run_cv.py (trained on the whole analysis set with the
hyper-parameters that won the nested block CV) and writes
  shap_summary.png            beeswarm of the top features, RF and XGB side by side
  shap_top10_dependence.png   dependence scatter of the 10 most important features
  feature_importance.tsv      RF impurity importance + XGB gain + mean |SHAP|
"""
import argparse
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import shap

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modeling_lib import FEATURE_SETS, OUTDIR, load_points, make_targets

plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.titlesize": 10,
                     "axes.grid": True, "grid.alpha": 0.25, "savefig.bbox": "tight"})

MAX_EXPLAIN = 10000        # SHAP is O(n * trees * depth); subsample the global targets
BG = 2000                  # background sample handed to the explainer


def load_model(tname, fs, name, tag=""):
    import joblib
    fn = os.path.join(OUTDIR, "models", f"{tname}__{fs}__{name}{tag}.joblib")
    if not os.path.exists(fn):
        return None
    return joblib.load(fn), fn


def tree_importance(pipeline):
    clf = pipeline.named_steps["clf"]
    if hasattr(clf, "feature_importances_"):
        return np.asarray(clf.feature_importances_, dtype=float)
    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tag", default="")
    ap.add_argument("--seed", type=int, default=42)
    args = ap.parse_args()

    df = load_points()
    targets = make_targets(df)
    rows = []
    pair_specs = [("T1_qingke_global", "A"), ("T3_qingke_plateau", "B")]

    for tname, fs in pair_specs:
        tinfo = targets[tname]
        X_all = df.loc[tinfo["mask"], FEATURE_SETS[fs]].reset_index(drop=True)
        y_all = pd.Series(tinfo["y"])
        n = len(X_all)
        if n > MAX_EXPLAIN:
            rng = np.random.default_rng(args.seed)
            idx = rng.choice(n, MAX_EXPLAIN, replace=False)
        else:
            idx = np.arange(n)
        Xs = X_all.iloc[idx].reset_index(drop=True)
        ys = y_all.iloc[idx].reset_index(drop=True)

        for name in ("RF", "XGB"):
            loaded = load_model(tname, fs, name, args.tag)
            if loaded is None:
                print(f"[skip] no model for {tname}/{fs}/{name}")
                continue
            bundle, fn = loaded
            pipe = bundle["model"]
            feats = bundle["features"]
            imp = pipe.named_steps["imp"]
            clf = pipe.named_steps["clf"]
            Ximp = pd.DataFrame(imp.transform(Xs[feats]), columns=feats)

            expl = shap.TreeExplainer(clf, feature_perturbation="tree_path_dependent")
            sv = expl.shap_values(Ximp, check_additivity=False)
            if isinstance(sv, list):                       # RF in some shap versions
                sv = sv[1] if len(sv) > 1 else sv[0]
            sv = np.asarray(sv)
            if sv.ndim == 3:                               # (n, features, classes)
                sv = sv[:, :, 1]
            print(f"[shap] {tname}/{fs}/{name}: values {sv.shape} on {len(Ximp)} rows")

            mean_abs = np.abs(sv).mean(axis=0)
            order = np.argsort(mean_abs)[::-1]
            np.save(os.path.join(OUTDIR, f"shap_values_{tname}_{fs}_{name}{args.tag}.npy"), sv)

            # ---- global importance table
            gini = tree_importance(pipe)
            for j, f in enumerate(feats):
                rows.append(dict(target=tname, feature_set=fs, model=name, feature=f,
                                 mean_abs_shap=mean_abs[j],
                                 shap_rank=int(np.where(order == j)[0][0]) + 1,
                                 rf_gini_importance=(gini[j] if gini is not None else np.nan),
                                 xgb_gain=(float(clf.feature_importances_[j])
                                           if name == "XGB" else np.nan),
                                 n_explained=len(Ximp)))

            # ---- beeswarm
            plt.figure(figsize=(7.0, 6.2))
            shap.summary_plot(sv, Ximp, max_display=15, show=False, plot_size=None)
            plt.title(f"{name}: SHAP value (impact on qingke-candidate probability)\n"
                      f"{tname} | feature set {fs} | n={len(Ximp)}")
            plt.xlabel("SHAP value (log-odds)")
            plt.tight_layout()
            plt.savefig(os.path.join(OUTDIR, f"shap_beeswarm_{tname}_{fs}_{name}{args.tag}.png"))
            plt.close()

            # ---- top-10 dependence
            top = [feats[j] for j in order[:10]]
            fig, axes = plt.subplots(2, 5, figsize=(19, 7.2))
            for ax, f in zip(axes.ravel(), top):
                j = feats.index(f)
                color_feat = feats[order[0]] if f != feats[order[0]] else feats[order[1]]
                cj = feats.index(color_feat)
                sc = ax.scatter(Ximp[f], sv[:, j], c=Ximp[color_feat], cmap="coolwarm",
                                s=4, alpha=0.5, linewidths=0)
                ax.axhline(0, color="k", lw=0.6, ls="--")
                ax.set_xlabel(f)
                ax.set_ylabel(f"SHAP({f})")
                ax.set_title(f, fontsize=9)
                cb = plt.colorbar(sc, ax=ax, pad=0.02)
                cb.set_label(color_feat, fontsize=7)
                cb.ax.tick_params(labelsize=6)
            fig.suptitle(f"{name}: SHAP dependence, top-10 features — {tname} | feature set {fs} "
                         f"(colour = second feature)", fontsize=12)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            fig.savefig(os.path.join(OUTDIR, f"shap_top10_dependence_{tname}_{fs}_{name}"
                                             f"{args.tag}.png"))
            plt.close(fig)

    if rows:
        imp_df = pd.DataFrame(rows)
        fn = os.path.join(OUTDIR, f"feature_importance{args.tag}.tsv")
        imp_df.sort_values(["target", "feature_set", "model", "shap_rank"]).to_csv(
            fn, sep="\t", index=False)
        print(f"[write] {fn} ({len(imp_df)} rows)")

    # canonical deliverable names -> the primary model (XGB) on the main target (T1, set A)
    import shutil
    for primary in ("XGB", "RF"):
        src = os.path.join(OUTDIR, f"shap_beeswarm_T1_qingke_global_A_{primary}{args.tag}.png")
        if os.path.exists(src):
            shutil.copyfile(src, os.path.join(OUTDIR, "shap_summary.png"))
            shutil.copyfile(
                os.path.join(OUTDIR,
                             f"shap_top10_dependence_T1_qingke_global_A_{primary}{args.tag}.png"),
                os.path.join(OUTDIR, "shap_top10_dependence.png"))
            print(f"[write] shap_summary.png / shap_top10_dependence.png (from {primary}, T1 set A)")
            break
    print("[done] shap stage")


if __name__ == "__main__":
    main()
