"""Figure 5 - Qingke-candidate suitability surface over the plateau box.

Model
-----
L2 logistic regression (sklearn) on the 1,002 plateau-box records: 558 Qingke
candidates (positive) vs 444 plateau-box records that are not candidates (negative).
Features are standardised.

Two feature sets are fitted and mapped, because the two answer different questions:

  * "with elevation"  - wc_bio1, wc_bio5, wc_bio6, wc_bio12, wc_bio15 + elev
  * "climate only"    - the same five bioclimatic variables, no elevation

`is_qingke_candidate` is DEFINED as `on_plateau & elev >= 2000` m, so elevation is a
tautological predictor and the "with elevation" model can only recover the label
rule (it is reported for completeness, not as evidence). The climate-only model is
the interpretable surface.

Validation: (i) repeated stratified 5-fold CV and (ii) spatially blocked CV, where
folds are contiguous 5-degree longitude blocks, so a fold's test points cannot be
spatially adjacent to its training points. The block estimate is the honest one here
because the records are strongly spatially autocorrelated.

Outputs: maps/fig5_suitability.png, maps/fig5_model_metrics.tsv
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib import colors as mcolors
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import GroupKFold, RepeatedStratifiedKFold
from sklearn.pipeline import make_pipeline
from sklearn.preprocessing import StandardScaler

import fig_common as fc

fc.use_style()

CLIM = ["wc_bio1", "wc_bio5", "wc_bio6", "wc_bio12", "wc_bio15"]
FEATSETS = {"climate_only": CLIM, "with_elevation": CLIM + ["elev"]}

NROWS, NCOLS = 300, 700
SHAPE = (NROWS, NCOLS)
BOX = fc.PLATEAU
KEY = "wc_bio1"        # geography key: any layer on the common WorldClim grid

df = fc.load_matrix(usecols=["decimalLongitude", "decimalLatitude", "occurrenceStatus", "elev",
                             "on_plateau", "is_qingke_candidate"] + CLIM)
box = df[df.on_plateau].dropna(subset=CLIM + ["elev"]).copy()
box["y"] = box.is_qingke_candidate.astype(int)
print(f"training records: {len(box)}  positive={int(box.y.sum())}  "
      f"negative={int((1 - box.y).sum())}  "
      f"prevalence={box.y.mean():.3f}")

# ---- optional: reuse an externally trained model if the modelling agent left one
EXT = os.path.join(fc.PROC, "modeling")
ext_models = []
if os.path.isdir(EXT):
    for root, _, files in os.walk(EXT):
        ext_models += [os.path.join(root, f) for f in files
                       if f.endswith((".joblib", ".pkl"))]
print(f"external model artefacts found: {len(ext_models)} "
      f"{[os.path.basename(m) for m in ext_models]}")
cv_metrics = os.path.join(EXT, "cv_metrics.tsv")
print(f"modeling/cv_metrics.tsv present: {os.path.exists(cv_metrics)}")
MODEL_PROVENANCE = ("logistic regression fitted by fig5_suitability.py "
                    "(no cv_metrics.tsv and no non-smoke trained model in "
                    "gis/processed/modeling/ at build time)")


def make_model(C):
    return make_pipeline(StandardScaler(),
                         LogisticRegression(C=C, max_iter=5000,
                                            class_weight="balanced"))


def lon_blocks(lon):
    """Contiguous 5-degree longitude blocks - the spatial CV unit."""
    return np.floor((lon - BOX[0]) / 5.0).astype(int)


def repeated_cv(X, y, C, n_splits=5, n_repeats=20, seed=0):
    rskf = RepeatedStratifiedKFold(n_splits=n_splits, n_repeats=n_repeats,
                                   random_state=seed)
    aucs = []
    for tr, te in rskf.split(X, y):
        m = make_model(C).fit(X[tr], y[tr])
        aucs.append(roc_auc_score(y[te], m.predict_proba(X[te])[:, 1]))
    return np.array(aucs)


def blocked_oof(X, y, blk, C, n_splits=5):
    """Out-of-fold predictions with whole 5-degree longitude blocks held out."""
    oof = np.full(len(y), np.nan)
    folds = list(GroupKFold(n_splits=n_splits).split(X, y, blk))
    for tr, te in folds:
        oof[te] = make_model(C).fit(X[tr], y[tr]).predict_proba(X[te])[:, 1]
    detail = [(sorted(set(blk[te])), int(len(te)),
               round(roc_auc_score(y[te], oof[te]), 4)) for _, te in folds]
    return oof, detail


# --- pick C by spatially blocked log-loss (a proper scoring rule, so the chosen
# --- model is the best-calibrated one rather than the most overconfident)
C_GRID = [0.001, 0.003, 0.01, 0.03, 0.1, 0.3, 1.0]
metrics, fitted = [], {}
for name, feats in FEATSETS.items():
    X = box[feats].to_numpy(float)
    y = box.y.to_numpy()
    blk = lon_blocks(box.decimalLongitude.to_numpy())

    scan = []
    for C in C_GRID:
        oof, _ = blocked_oof(X, y, blk, C)
        scan.append(dict(C=C, logloss=log_loss(y, oof), brier=brier_score_loss(y, oof),
                         auc=roc_auc_score(y, oof)))
    scan = pd.DataFrame(scan)
    best = scan.loc[scan.logloss.idxmin()]
    C = float(best.C)
    print(f"\n[{name}] features={feats}")
    print(scan.round(4).to_string(index=False))
    print(f"  -> C={C} selected (minimum blocked log-loss)")

    oof, blk_detail = blocked_oof(X, y, blk, C)
    rep = repeated_cv(X, y, C)
    full = make_model(C).fit(X, y)
    fitted[name] = (full, feats)
    coefs = full[-1].coef_[0]
    print(f"  blocked 5-fold (5-deg lon blocks) AUC = {best.auc:.4f} "
          f"log-loss = {best.logloss:.4f} Brier = {best.brier:.4f}")
    print(f"  repeated stratified 5-fold CV AUC = {rep.mean():.4f} +/- {rep.std():.4f}")
    print(f"  apparent (in-sample) AUC = "
          f"{roc_auc_score(y, full.predict_proba(X)[:, 1]):.4f}")
    print(f"  block detail (lon blocks, n_test, AUC): {blk_detail}")
    print("  standardised coefficients: "
          + ", ".join(f"{f}={c:+.2f}" for f, c in zip(feats, coefs)))
    metrics.append(dict(model=name, C=C, n_features=len(feats), n=len(box),
                        n_positive=int(y.sum()), n_negative=int((1 - y).sum()),
                        auc_blocked_cv=round(float(best.auc), 4),
                        logloss_blocked_cv=round(float(best.logloss), 4),
                        brier_blocked_cv=round(float(best.brier), 4),
                        auc_repeated_cv=round(rep.mean(), 4),
                        auc_repeated_cv_sd=round(rep.std(), 4),
                        auc_apparent=round(float(roc_auc_score(
                            y, full.predict_proba(X)[:, 1])), 4),
                        coefs="; ".join(f"{f}={c:+.3f}" for f, c in zip(feats, coefs))))

mt = pd.DataFrame(metrics)
mt.to_csv(os.path.join(fc.OUT, "fig5_model_metrics.tsv"), sep="\t", index=False)

# ---- predict every 0.05 deg cell in the box ---------------------------------
grid = {}
for var in CLIM + ["elev"]:
    fname = {"elev": "wc2.1_30s_elev.tif"}.get(var, f"{var.replace('wc_bio', 'wc2.1_30s_bio_')}.tif")
    arr, extent = fc.read_window(os.path.join(fc.WC, fname), BOX, SHAPE)
    arr = np.where(np.abs(arr) > 1e30, np.nan, arr)
    if var == "elev":
        arr = np.where(arr < -100, np.nan, arr)
    grid[var] = arr
    print(f"grid {var}: {fname} valid={int(np.isfinite(arr).sum())}/{arr.size}")

valid = np.all([np.isfinite(grid[v]) for v in CLIM + ["elev"]], axis=0)
Xg = np.column_stack([grid[v][valid] for v in CLIM + ["elev"]])
print(f"predicting {Xg.shape[0]} cells")

surfaces = {}
for name, (model, feats) in fitted.items():
    Xs = np.column_stack([grid[v][valid] for v in feats])
    s = np.full(valid.shape, np.nan)
    s[valid] = model.predict_proba(Xs)[:, 1]
    surfaces[name] = s
    print(f"{name}: p range {np.nanmin(s):.3f}-{np.nanmax(s):.3f}, "
          f"median {np.nanmedian(s):.3f}")

# ---- figure -----------------------------------------------------------------
# The logistic fit is close to separable, so the surface is sharply bounded and most
# cells sit near 0 or 1. The p = 0.5 contour is therefore drawn explicitly - it is the
# readable content of the map, and it keeps the saturated fill from hiding the edge.
qk = df[df.is_qingke_candidate]
fig, axes = plt.subplots(2, 1, figsize=(7.2, 6.9))
ORDER = [("climate_only", "(a) Climate-only model (no elevation)"),
         ("with_elevation", "(b) Model including elevation (near-tautological)")]
NORM = mcolors.Normalize(0, 1)
lon_v, lat_v = fc.lons_lats(BOX, NCOLS, NROWS)

for ax, (name, title) in zip(axes, ORDER):
    im = ax.imshow(surfaces[name], extent=extent, origin="upper", cmap=fc.PROB_CMAP,
                   norm=NORM, interpolation="nearest", aspect="equal")
    im.cmap.set_bad("#f0efec")
    ax.contour(lon_v, lat_v, surfaces[name], levels=[0.5], colors=["white"],
               linewidths=1.3, zorder=3)
    ax.scatter(qk.decimalLongitude, qk.decimalLatitude, s=4.5, c=fc.INK, linewidths=0,
               alpha=0.8, zorder=4,
               label=f"Qingke candidates (n={len(qk)})" if name == "climate_only" else None)
    m = mt[mt.model == name].iloc[0]
    ax.set_title(title, loc="left", pad=6)
    ax.text(0.978, 0.975,
            f"repeated CV AUC = {m.auc_repeated_cv:.3f}\n"
            f"blocked (5-deg lon) CV AUC = {m.auc_blocked_cv:.3f}\n"
            f"white line = p = 0.5 boundary",
            transform=ax.transAxes, fontsize=7.5, ha="right", va="top",
            bbox=dict(boxstyle="round,pad=0.3", fc="white", ec=fc.GRID, lw=0.6),
            zorder=5)
    ax.set_xlabel("Longitude (degrees east)")
    ax.set_ylabel("Latitude (degrees north)")
    ax.set_xticks(np.arange(70, 106, 5))
    ax.set_yticks(np.arange(25, 41, 5))
    ax.set_xlim(BOX[0], BOX[2])
    ax.set_ylim(BOX[1], BOX[3])
    ax.grid(False)
    if name == "climate_only":
        ax.legend(loc="lower left", fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.012)
    cb.set_label("Predicted probability of being a Qingke candidate", fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0.5)

plt.tight_layout()
out = os.path.join(fc.OUT, "fig5_suitability.png")
fig.savefig(out, dpi=300)
plt.close(fig)
print("\nwrote", out)
