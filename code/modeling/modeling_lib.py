"""
Shared library for the highland-barley environmental prediction models (TODO stage 5).

Design summary
--------------
Targets
  T1  qingke_global : y = is_qingke_candidate, all points with known elevation (n=45,701; 558 pos)
  T2  highalt_global: y = elev >= 2000 m,       all points with known elevation (n=45,701; 1,048 pos)
  T3  plateau       : y = is_qingke_candidate, Tibetan-plateau rectangle subset (n=1,002; 558 pos)

Feature sets
  A  : wc_bio1..19 + elev                                  (20 vars, global coverage)
  A' : wc_bio1..19                                         (19 vars, elev removed -> ablation)
  B  : elev + slope + aspect(sin/cos) + wc_bio1..19 + sg_* (16 soil)   (regional only)
  B' : slope + aspect(sin/cos) + wc_bio1..19 + sg_*        (elev removed -> ablation)

Blocking
  Spatial blocks are k-means clusters on standardised (lon, lat) of each target's analysis set.
  Two scales are used:
    G (global)  : k=10 over all analysis points -> used for T1/T2.  Because all 558 qingke
                  candidates sit in a single contiguous region, the global blocking leaves them
                  in 1-2 blocks; folds without positives are marked invalid.
    P (plateau) : k=10 over the 1,002 plateau points -> used for T3, and for the "within-plateau"
                  version of T1 that resolves the positive region into 9-10 test blocks.
  Cross-validation is LeaveOneBlockOut.  Imputation / standardisation / collinearity pruning are
  re-fitted inside every training fold, and hyper-parameters are tuned with an inner GroupKFold
  over the training blocks only (no spatial leakage).

The F1-optimal decision threshold is taken from the inner-CV out-of-fold probabilities of the
selected configuration (never from the held-out block).
"""

import os
import numpy as np
import pandas as pd
from sklearn.base import clone
from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import (average_precision_score, balanced_accuracy_score,
                             f1_score, roc_auc_score)
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

MATRIX = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed/barley_env_matrix.tsv.gz')
OUTDIR = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed/modeling')

WC = [f"wc_bio{i}" for i in range(1, 20)]
SOIL = [f"sg_{v}_{d}" for v in ("clay", "sand", "silt", "soc", "phh2o", "nitrogen", "bdod", "cec")
        for d in ("0_5", "5_15")]
TERRAIN = ["elev", "slope", "aspect_sin", "aspect_cos"]

FEATURE_SETS = {
    "A": WC + ["elev"],
    "A_noelev": WC,
    "B": TERRAIN + WC + SOIL,
    "B_noelev": ["slope", "aspect_sin", "aspect_cos"] + WC + SOIL,
}

MODEL_NAMES = ["LR", "RF", "XGB"]

# Threading: XGBoost's OpenMP backend thrashes badly with n_jobs=-1 on this 16-core box
# (measured 115 s vs 0.3 s per fit); 4 threads is the sweet spot for both tree libraries.
N_JOBS = 4

# ---------------------------------------------------------------- data


def load_points(verbose=True):
    """Load the environmental matrix and build the derived modelling columns."""
    df = pd.read_csv(
        MATRIX, sep="\t", low_memory=False,
        usecols=["decimalLongitude", "decimalLatitude", "elev", "slope", "aspect",
                 "on_plateau", "is_qingke_candidate", "countryCode", "year", "occurrenceStatus",
                 "scientificName"] + WC + SOIL,
    )
    n_all = len(df)
    df = df[df["occurrenceStatus"].fillna("PRESENT").eq("PRESENT")].copy()
    df.drop(columns=["occurrenceStatus"], inplace=True)
    df["lon"] = df["decimalLongitude"].astype(float)
    df["lat"] = df["decimalLatitude"].astype(float)
    asp = np.deg2rad(df["aspect"].astype(float))
    df["aspect_sin"] = np.sin(asp)
    df["aspect_cos"] = np.cos(asp)
    # aspect is undefined on flat ground -> sin/cos follow the NaN pattern
    df.loc[df["aspect"].isna(), ["aspect_sin", "aspect_cos"]] = np.nan
    if verbose:
        print(f"[data] {len(df)} PRESENT points (excluded {n_all-len(df)} ABSENT), "
              f"proxy={int(df['is_qingke_candidate'].sum())}, plateau={int(df['on_plateau'].sum())}")
    return df


def make_targets(df):
    """Return {target_name: (row_mask, y_values, analysis_frame, feature_sets, block_scale)}."""
    known_elev = df["elev"].notna()
    plateau = df["on_plateau"].astype(bool)
    targets = {}

    m = known_elev.values
    targets["T1_qingke_global"] = dict(
        mask=m, y=df.loc[m, "is_qingke_candidate"].astype(int).values,
        feature_sets=["A", "A_noelev"], scale="G",
        label="T1 qingke candidate (global)")
    targets["T2_highalt_global"] = dict(
        mask=m, y=(df.loc[m, "elev"] >= 2000).astype(int).values,
        feature_sets=["A", "A_noelev"], scale="G",
        label="T2 high altitude >=2000 m (global)")
    targets["T3_qingke_plateau"] = dict(
        mask=plateau.values, y=df.loc[plateau, "is_qingke_candidate"].astype(int).values,
        feature_sets=["B", "B_noelev", "A", "A_noelev"], scale="P",
        label="T3 qingke candidate (plateau rectangle)")
    for t in targets.values():
        t["n"] = int(t["mask"].sum())
        t["n_pos"] = int(t["y"].sum())
    return targets


# ---------------------------------------------------------------- blocks


def spatial_blocks(coords, k=10, seed=42):
    """k-means over standardised (lon, lat)."""
    X = StandardScaler().fit_transform(coords)
    km = KMeans(n_clusters=k, random_state=seed, n_init=10).fit(X)
    return km.labels_


def block_distribution(y, blocks):
    """Positives per block: how many blocks carry positives, max share, top block."""
    s = pd.Series(np.asarray(y), index=np.arange(len(y))).groupby(blocks).sum()
    n_pos = int(np.asarray(y).sum())
    sizes = pd.Series(np.ones(len(y))).groupby(blocks).size()
    out = pd.DataFrame({"block": s.index, "n_points": sizes.values, "n_pos": s.values.astype(int)})
    out["pos_share_pct"] = (out["n_pos"] / max(n_pos, 1) * 100).round(2)
    out = out.sort_values("n_pos", ascending=False).reset_index(drop=True)
    summary = dict(n_blocks=len(out), blocks_with_pos=int((out["n_pos"] > 0).sum()),
                   blocks_with_ge5_pos=int((out["n_pos"] >= 5).sum()),
                   n_pos=n_pos,
                   max_block_share_pct=round(float(out["pos_share_pct"].max()), 2),
                   top_block=int(out.loc[0, "block"]) if len(out) else -1,
                   top_block_n_pos=int(out.loc[0, "n_pos"]) if len(out) else 0)
    return out, summary


# ---------------------------------------------------------------- collinearity pruning


def prune_collinear(X: pd.DataFrame, threshold=0.9):
    """Greedy removal of one variable from every |r| >= threshold pair.

    When a pair exceeds the threshold the variable with the higher mean |r| against all the
    remaining variables is dropped (ties broken towards keeping the earlier column).
    Returns (kept_columns, dropped_columns).
    """
    cols = list(X.columns)
    corr = X.corr(numeric_only=True).abs()
    dropped = []
    while True:
        c = corr.loc[cols, cols].to_numpy(copy=True)
        np.fill_diagonal(c, 0.0)
        if c.size == 0 or not np.isfinite(c).any() or np.nanmax(c) < threshold:
            break
        i, j = np.unravel_index(np.nanargmax(c), c.shape)
        a, b = cols[i], cols[j]
        mean_a = float(np.nanmean(corr.loc[a, cols].values))
        mean_b = float(np.nanmean(corr.loc[b, cols].values))
        drop = b if mean_a >= mean_b else a
        cols.remove(drop)
        dropped.append(drop)
    return cols, dropped


# ---------------------------------------------------------------- models


def make_estimator(name, params, seed=42, n_jobs=N_JOBS):
    """Preprocessing + classifier pipeline (imputation is fitted inside the fold)."""
    if name == "LR":
        clf = LogisticRegression(max_iter=3000, class_weight="balanced",
                                 random_state=seed, **params)
        return Pipeline([("imp", SimpleImputer(strategy="median")),
                         ("sc", StandardScaler()), ("clf", clf)])
    if name == "RF":
        clf = RandomForestClassifier(n_estimators=300, n_jobs=n_jobs, random_state=seed,
                                     class_weight="balanced_subsample", **params)
        return Pipeline([("imp", SimpleImputer(strategy="median")), ("clf", clf)])
    if name == "XGB":
        from xgboost import XGBClassifier
        clf = XGBClassifier(n_estimators=300, learning_rate=0.1, subsample=0.8,
                            colsample_bytree=0.8, tree_method="hist", n_jobs=n_jobs,
                            random_state=seed, eval_metric="logloss", **params)
        return Pipeline([("imp", SimpleImputer(strategy="median")), ("clf", clf)])
    raise ValueError(name)


def param_grid(name, pos_weight):
    """Deliberately small grid: it is searched inside every training fold."""
    if name == "LR":
        return [{"C": c} for c in (0.03, 0.3, 3.0)]
    if name == "RF":
        return [{"max_depth": None, "min_samples_leaf": 1},
                {"max_depth": 12, "min_samples_leaf": 5}]
    if name == "XGB":
        return [{"max_depth": 4, "scale_pos_weight": pos_weight},
                {"max_depth": 6, "scale_pos_weight": pos_weight}]
    raise ValueError(name)


def metrics_at(y_true, proba, threshold):
    pred = (proba >= threshold).astype(int)
    out = {"f1": f1_score(y_true, pred, zero_division=0),
           "bacc": balanced_accuracy_score(y_true, pred)}
    if len(np.unique(y_true)) > 1:
        out["auc"] = roc_auc_score(y_true, proba)
        out["pr_auc"] = average_precision_score(y_true, proba)
    else:
        out["auc"] = np.nan
        out["pr_auc"] = np.nan
    return out


def best_f1_threshold(y_true, proba):
    """Threshold maximising F1 on the supplied (in-training) probabilities.

    Vectorised: sorts once, then sweeps candidate cut-offs with cumulative sums, so the cost is
    O(n log n) instead of O(n * #unique).
    """
    y = np.asarray(y_true).astype(int)
    p = np.asarray(proba, dtype=float)
    ok = np.isfinite(p)
    y, p = y[ok], p[ok]
    if len(y) == 0 or y.sum() == 0:
        return 0.5, -1.0
    order = np.argsort(p, kind="mergesort")
    p, y = p[order], y[order]
    n, total_pos = len(y), int(y.sum())
    csum = np.concatenate([[0], np.cumsum(y)])   # csum[i] = positives with p < p[i]
    # a cut-off only changes between distinct scores -> evaluate the first index of each tie group
    idx = np.flatnonzero(np.r_[True, p[1:] != p[:-1]])
    below = csum[idx]                            # positives predicted negative at t = p[idx]
    tp = total_pos - below                       # rule: predicted positive iff proba >= t
    pred_pos = n - idx
    fp = pred_pos - tp
    fn = below
    with np.errstate(invalid="ignore", divide="ignore"):
        f1 = 2 * tp / (2 * tp + fp + fn)
    f1 = np.nan_to_num(f1)
    k = idx[int(np.argmax(f1))]
    return float(p[k]), float(f1[int(np.argmax(f1))])
