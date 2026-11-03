#!/usr/bin/env python3
"""Record which variables the |r| >= 0.9 pruning removes from each feature set.

Written to the output dir as collinearity_pruned_full.tsv; the CV itself re-runs the same
procedure inside every training fold (see collinearity_dropped.tsv).
"""
import os
import sys

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from modeling_lib import FEATURE_SETS, OUTDIR, load_points, make_targets, prune_collinear


def main():
    df = load_points()
    tg = make_targets(df)
    rows = []
    for tname in ("T1_qingke_global", "T2_highalt_global", "T3_qingke_plateau"):
        sub = df[tg[tname]["mask"]]
        for fs in tg[tname]["feature_sets"]:
            X = sub[FEATURE_SETS[fs]]
            kept, dropped = prune_collinear(X)
            print(f"{tname} | {fs}: {len(X.columns)} -> {len(kept)} vars; dropped: {dropped}")
            rows.append(dict(target=tname, feature_set=fs, n_input=len(X.columns),
                             n_kept=len(kept), kept=";".join(kept), dropped=";".join(dropped)))
    pd.DataFrame(rows).to_csv(os.path.join(OUTDIR, "collinearity_pruned_full.tsv"),
                              sep="\t", index=False)
    print("[write] collinearity_pruned_full.tsv")


if __name__ == "__main__":
    main()
