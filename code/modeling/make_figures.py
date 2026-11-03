#!/usr/bin/env python3
"""Stage 5e: CV result figures (per-fold AUC distributions + mean AUC by model)."""
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modeling_lib import OUTDIR

plt.rcParams.update({"figure.dpi": 130, "font.size": 9, "axes.grid": True,
                     "grid.alpha": 0.25, "savefig.bbox": "tight"})
COLORS = {"LR": "#4C72B0", "RF": "#55A868", "XGB": "#C44E52"}


def main(tag=""):
    cv = pd.read_csv(os.path.join(OUTDIR, f"cv_metrics{tag}.tsv"), sep="\t")
    v = cv[cv["valid"]].copy()
    short = {"T1_qingke_global": "T1", "T2_highalt_global": "T2", "T3_qingke_plateau": "T3"}
    fs_short = {"A": "A", "A_noelev": "A'", "B": "B", "B_noelev": "B'"}
    v["combo"] = v["target"].map(short) + "\u00b7" + v["feature_set"].map(fs_short)
    combos = sorted(v["combo"].unique())
    models = ["LR", "RF", "XGB"]

    fig, axes = plt.subplots(1, 3, figsize=(18, 5.4),
                             gridspec_kw={"width_ratios": [2.4, 1.6, 1.2]})
    # --- per-fold AUC
    ax = axes[0]
    rng = np.random.default_rng(0)
    for i, c in enumerate(combos):
        for j, m in enumerate(models):
            d = v[(v["combo"] == c) & (v["model"] == m)]["auc"].values
            if not len(d):
                continue
            x = i * 3.4 + j + rng.uniform(-0.12, 0.12, len(d))
            ax.scatter(x, d, s=26, color=COLORS[m], alpha=0.85, zorder=3,
                       label=m if i == 0 else None)
            ax.hlines(np.mean(d), i * 3.4 + j - 0.32, i * 3.4 + j + 0.32, color="k", lw=1.6,
                      zorder=4)
    ax.set_xticks([i * 3.4 + j for i in range(len(combos)) for j in range(3)])
    ax.set_xticklabels([m for _ in combos for m in models], fontsize=7)
    for i, c in enumerate(combos):
        ax.text(i * 3.4 + 1, 0.035, c, ha="center", va="bottom", fontsize=8.5,
                transform=ax.get_xaxis_transform(), fontweight="bold",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.85, pad=1.5))
    ax.set_ylabel("ROC-AUC (held-out spatial block)")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title("Per-fold spatial block CV AUC (bar = fold mean; A/B keep elev, A'/B' drop it)")
    ax.legend(loc="center left", fontsize=8, framealpha=0.95)

    # --- mean AUC by model / target
    ax = axes[1]
    tgts = sorted(v["target"].unique())
    width = 0.26
    for j, m in enumerate(models):
        means, sds, xs = [], [], []
        for i, t in enumerate(tgts):
            d = v[(v["target"] == t) & (v["model"] == m) & (v["feature_set"].str.endswith("noelev") == False)]  # noqa: E712
            d = d[d["feature_set"].isin(["A", "B"])]
            if not len(d):
                continue
            means.append(d["auc"].mean())
            sds.append(d["auc"].std())
            xs.append(i + (j - 1) * width)
        ax.bar(xs, means, width=width, yerr=sds, capsize=3, color=COLORS[m], label=m)
    ax.set_xticks(range(len(tgts)))
    ax.set_xticklabels([t.replace("_", "\n") for t in tgts], fontsize=8)
    ax.set_ylim(0.72, 1.03)
    ax.set_ylabel("mean ROC-AUC (+/- sd over folds)")
    ax.set_title("Mean AUC, feature sets A / B")
    ax.legend(fontsize=8, loc="lower right")

    # --- valid folds per combo (one row per target x feature set, models collapsed)
    ax = axes[2]
    cnt = (cv.drop_duplicates(["target", "feature_set", "fold"])
           .groupby(["target", "feature_set"])["valid"].agg(["sum", "size"]).reset_index())
    lbl = cnt["target"].str.replace("_global", "").str.replace("_plateau", "") + " " + \
          cnt["feature_set"]
    ax.barh(range(len(cnt)), cnt["size"], color="#DDDDDD", label="blocks of the partition")
    ax.barh(range(len(cnt)), cnt["sum"], color="#4C72B0", label="folds with >=5 positives")
    ax.set_yticks(range(len(cnt)))
    ax.set_yticklabels(lbl, fontsize=7)
    ax.invert_yaxis()
    ax.set_xlabel("number of folds")
    ax.set_title("Spatial CV fold coverage")
    ax.legend(fontsize=8, loc="lower right")

    fig.suptitle("Environmental prediction models for qingke/highland barley — spatial block CV",
                 fontsize=13)
    fig.tight_layout(rect=(0, 0.055, 1, 0.95))
    out = os.path.join(OUTDIR, f"cv_auc_comparison{tag}.png")
    fig.savefig(out)
    print("[write]", out)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else "")
