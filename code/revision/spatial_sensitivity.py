"""Matched climate-product and distance-thinning sensitivity analyses.

Uses the original plateau labels and outer blocks. Preprocessing, class weights
and AP-based hyperparameter selection are fitted inside each inner training fold.
This is a new sensitivity run, not a replacement for archived ROC-AUC tuning.
"""
import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.metrics.pairwise import haversine_distances
from sklearn.model_selection import GroupKFold
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler


class CorrelationFilter(BaseEstimator, TransformerMixin):
    """Choose columns on each training partition, after median imputation."""
    def fit(self, X, y=None):
        corr = pd.DataFrame(X).corr().abs().to_numpy()
        keep = list(range(X.shape[1]))
        while len(keep) > 1:
            current = corr[np.ix_(keep, keep)].copy()
            np.fill_diagonal(current, 0)
            if not np.isfinite(current).any() or np.nanmax(current) < .9:
                break
            i, j = np.unravel_index(np.nanargmax(current), current.shape)
            a, b = keep[i], keep[j]
            drop = b if np.nanmean(corr[a, keep]) >= np.nanmean(corr[b, keep]) else a
            keep.remove(drop)
        self.keep_ = keep
        return self

    def transform(self, X):
        return np.asarray(X)[:, self.keep_]


def estimator(name, params, y):
    prep = [("impute", SimpleImputer(strategy="median", keep_empty_features=True))]
    if name == "LR":
        prep += [("correlation", CorrelationFilter()), ("scale", StandardScaler())]
        classifier = LogisticRegression(max_iter=3000, class_weight="balanced", random_state=42, **params)
    elif name == "RF":
        classifier = RandomForestClassifier(n_estimators=300, n_jobs=4,
            class_weight="balanced_subsample", random_state=42, **params)
    else:
        from xgboost import XGBClassifier
        weight = float((y == 0).sum() / max((y == 1).sum(), 1))
        classifier = XGBClassifier(n_estimators=300, learning_rate=.1, subsample=.8,
            colsample_bytree=.8, tree_method="hist", n_jobs=4, random_state=42,
            eval_metric="logloss", scale_pos_weight=weight, **params)
    return Pipeline(prep + [("classifier", classifier)])


def grid(name):
    return ({"LR": [{"C": c} for c in [.03, .3, 3.]],
             "RF": [{"max_depth": None, "min_samples_leaf": 1},
                    {"max_depth": 12, "min_samples_leaf": 5}],
             "XGB": [{"max_depth": 4}, {"max_depth": 6}]}[name])


def tune(name, X, y, blocks):
    groups = np.unique(blocks)
    if len(groups) < 2:
        return grid(name)[0], "first_configuration_no_inner_folds"
    folds = list(GroupKFold(n_splits=min(3, len(groups))).split(X, y, blocks))
    scored = []
    for params in grid(name):
        scores = []
        for train, test in folds:
            if len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2:
                continue
            fit = estimator(name, params, y[train])
            fit.fit(X[train], y[train])
            scores.append(average_precision_score(y[test], fit.predict_proba(X[test])[:, 1]))
        scored.append((np.mean(scores) if scores else -np.inf, params))
    best = max(scored, key=lambda pair: pair[0])
    return best[1], "inner_grouped_AP" if np.isfinite(best[0]) else "first_configuration_no_valid_inner_folds"


def thin(frame, distance_km):
    if distance_km == 0:
        return np.arange(len(frame))
    priorities = [hashlib.sha256(f"42:{x}".encode()).hexdigest() for x in frame.gbif_key]
    order = np.argsort(priorities)
    radians = np.deg2rad(frame[["decimalLatitude", "decimalLongitude"]].to_numpy())
    distances = haversine_distances(radians) * 6371.0088
    kept = []
    for index in order:
        if not kept or np.all(distances[index, kept] >= distance_km):
            kept.append(index)
    return np.sort(kept)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--evidence", type=Path, required=True)
    ap.add_argument("--output", type=Path, required=True)
    args = ap.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.evidence / "gis/processed/barley_env_matrix.tsv.gz", sep="\t", low_memory=False)
    frame = frame[frame.occurrenceStatus.fillna("PRESENT").eq("PRESENT") & frame.on_plateau.astype(bool)].copy().reset_index(drop=True)
    original_blocks = np.load(args.evidence / "gis/processed/modeling/blocks_T3_qingke_plateau_present.npy")
    wc, ch = [f"wc_bio{i}" for i in range(1, 20)], [f"ch_bio{i}" for i in range(1, 20)]
    common = frame[wc + ch].notna().all(axis=1).to_numpy()
    all_rows, predictions, retained_rows = [], [], []
    scenarios = [("climate_product", product, 0., np.flatnonzero(common)) for product in ["WorldClim", "CHELSA"]]
    scenarios += [("spatial_thinning", "WorldClim", float(d), thin(frame, d)) for d in [0, 5, 10, 25]]
    started = time.monotonic()
    for analysis, product, distance, retained in scenarios:
        sub = frame.iloc[retained].copy()
        X = sub[wc if product == "WorldClim" else ch].to_numpy()
        y = sub.is_qingke_candidate.astype(int).to_numpy()
        blocks = original_blocks[retained]
        for local, index in enumerate(retained):
            retained_rows.append(dict(analysis=analysis, climate_product=product, thinning_km=distance,
                                      original_row=int(index), gbif_key=str(frame.iloc[index].gbif_key),
                                      block=int(blocks[local]), label=int(y[local])))
        for name in ["LR", "RF", "XGB"]:
            for block in np.unique(original_blocks):
                train, test = np.flatnonzero(blocks != block), np.flatnonzero(blocks == block)
                row = dict(analysis=analysis, climate_product=product, thinning_km=distance,
                           model=name, fold=int(block), cohort_n=len(sub), positives=int(y.sum()),
                           n_train=len(train), n_test=len(test), n_pos_train=int(y[train].sum()),
                           n_pos_test=int(y[test].sum()), valid=False, seed=42)
                if not len(test) or len(np.unique(y[train])) < 2 or len(np.unique(y[test])) < 2 or y[train].sum() < 5 or y[test].sum() < 5:
                    row.update(skip_reason="empty_or_single_class_or_fewer_than_5_positives")
                else:
                    params, tuning = tune(name, X[train], y[train], blocks[train])
                    fit = estimator(name, params, y[train])
                    fit.fit(X[train], y[train])
                    prob = fit.predict_proba(X[test])[:, 1]
                    row.update(valid=True, skip_reason="", ap=average_precision_score(y[test], prob),
                               roc_auc=roc_auc_score(y[test], prob), test_prevalence=float(y[test].mean()),
                               parameters=json.dumps(params), tuning=tuning)
                    for local, proba in zip(test, prob):
                        predictions.append(dict(analysis=analysis, climate_product=product, thinning_km=distance,
                                                model=name, fold=int(block), gbif_key=str(sub.iloc[local].gbif_key),
                                                label=int(y[local]), probability=float(proba)))
                all_rows.append(row)
            print(analysis, product, distance, name, f"{time.monotonic()-started:.0f}s", flush=True)
            pd.DataFrame(all_rows).to_csv(args.output / "spatial_sensitivity_folds.tsv", sep="\t", index=False)
    result = pd.DataFrame(all_rows)
    scored = result[result.valid]
    summary = scored.groupby(["analysis", "climate_product", "thinning_km", "model"]).agg(
        cohort_n=("cohort_n", "first"), positives=("positives", "first"),
        scorable_folds=("fold", "count"), ap_mean=("ap", "mean"), ap_sd=("ap", "std"),
        roc_auc_mean=("roc_auc", "mean"), roc_auc_sd=("roc_auc", "std"),
        test_prevalence_min=("test_prevalence", "min"), test_prevalence_max=("test_prevalence", "max")
    ).reset_index()
    summary.to_csv(args.output / "spatial_sensitivity_summary.tsv", sep="\t", index=False)
    pd.DataFrame(predictions).to_csv(args.output / "spatial_sensitivity_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(retained_rows).to_csv(args.output / "spatial_sensitivity_record_manifest.tsv", sep="\t", index=False)
    import sklearn, scipy, xgboost
    metadata = dict(seed=42, outer_blocks="archived plateau blocks; never redefined after thinning",
        common_climate_cohort_n=int(common.sum()), plateau_n=len(frame),
        spatial_thinning="greedy haversine distance; SHA256(42:gbif_key) ordering, independent of label",
        tuning="three grouped inner folds, AP, training-only preprocessing and class weights",
        minimum_positives_train_test=5, python=platform.python_version(),
        numpy=np.__version__, pandas=pd.__version__, scipy=scipy.__version__,
        sklearn=sklearn.__version__, xgboost=xgboost.__version__, elapsed_seconds=time.monotonic()-started)
    (args.output / "spatial_sensitivity_configuration.json").write_text(json.dumps(metadata, indent=2))
    print(summary.to_string(index=False), flush=True)


if __name__ == "__main__":
    main()
