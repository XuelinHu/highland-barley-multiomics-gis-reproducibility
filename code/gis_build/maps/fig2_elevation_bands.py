"""Figure 2 - 500 m elevation bands: records per band, Qingke-candidate share per
band, and the elevation distribution of the 558 Qingke candidates.

Panel (a) carries a second (right) y-axis, as specified. To keep the two scales
readable rather than misleading, each axis is colour-keyed to the single mark class
it describes (left = grey record bars, right = red Qingke share), the twin axis draws
no gridlines, and the three mark classes carry an explicit legend.

Output: gis/processed/maps/fig2_elevation_bands.png  (300 dpi)
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

import fig_common as fc

fc.use_style()

df = fc.load_matrix(usecols=["elev", "occurrenceStatus", "on_plateau", "is_qingke_candidate"])
q = df[df.is_qingke_candidate]

# 500 m bands: <1000, 1000-1500, ..., 4000-4500, >4500
EDGES = [-np.inf, 1000, 1500, 2000, 2500, 3000, 3500, 4000, 4500, np.inf]
LABELS = ["<1000", "1000-1500", "1500-2000", "2000-2500", "2500-3000",
          "3000-3500", "3500-4000", "4000-4500", ">4500"]
CENTERS = np.array([750, 1250, 1750, 2250, 2750, 3250, 3750, 4250, 4750], dtype=float)

band = pd.cut(df.elev, bins=EDGES, labels=LABELS, right=False)
qband = pd.cut(q.elev, bins=EDGES, labels=LABELS, right=False)
n_all = band.value_counts().reindex(LABELS).fillna(0).to_numpy()
n_q = qband.value_counts().reindex(LABELS).fillna(0).to_numpy()
share = np.divide(n_q, n_all, out=np.zeros_like(n_q, dtype=float),
                  where=n_all > 0) * 100.0

tbl = pd.DataFrame({"band": LABELS, "n_records": n_all.astype(int),
                    "n_qingke": n_q.astype(int),
                    "qingke_share_pct": np.round(share, 2)})
print(tbl.to_string(index=False))
print(f"total records in bands: {int(n_all.sum())}  (of {len(df)}; "
      f"{int(df.elev.isna().sum())} without elevation)")
print(f"Qingke median elev: {q.elev.median():.0f} m, mean {q.elev.mean():.0f} m, "
      f"range {q.elev.min():.0f}-{q.elev.max():.0f} m")

fig, (axa, axb) = plt.subplots(2, 1, figsize=(7.0, 7.4), sharex=True,
                               gridspec_kw=dict(height_ratios=[1.35, 1.0], hspace=0.13))

# ---------------------------------------------------------------- (a) counts + share
BW = 440.0   # bar width in metres (< 500 -> a 2 px surface gap between neighbours)
axa.bar(CENTERS, n_all, width=BW, color=fc.GREY, edgecolor=fc.SURFACE,
        linewidth=0.8, label=f"All records (n={int(n_all.sum()):,})", zorder=2)
axa.bar(CENTERS, n_q, width=BW * 0.52, color=fc.RED, edgecolor=fc.SURFACE,
        linewidth=0.6, label=f"Qingke candidates (n={int(n_q.sum())})", zorder=3)

axa.set_yscale("log")
axa.set_ylim(1, 2.2e5)
axa.set_ylabel("Number of records (log scale)", color=fc.INK2)
axa.tick_params(axis="y", colors=fc.INK2)
axa.grid(True, axis="y")
axa.set_title("(a) Records per 500 m elevation band and Qingke-candidate share",
              loc="left", pad=8)

axr = axa.twinx()
axr.grid(False)
axr.plot(CENTERS, share, color=fc.RED, linewidth=2.0, marker="o", markersize=5.5,
         markerfacecolor=fc.RED, markeredgecolor=fc.SURFACE, markeredgewidth=1.0,
         zorder=4, label="Qingke share (right axis)")
axr.set_ylim(-4, 122)
axr.set_ylabel("Share of records that are Qingke candidates (%)", color=fc.RED)
axr.tick_params(axis="y", colors=fc.RED)
# direct labels on the share curve instead of a value on every point
for x, s, nn in zip(CENTERS, share, n_q):
    if s > 0:
        axr.annotate(f"{s:.0f}%", (x, s), textcoords="offset points", xytext=(0, 8),
                     ha="center", fontsize=7.5, color=fc.RED)

h1, l1 = axa.get_legend_handles_labels()
h2, l2 = axr.get_legend_handles_labels()
axa.legend(h1 + h2, l1 + l2, loc="upper left", fontsize=8, labelspacing=0.4)
axa.axvline(fc.ELEV_MIN_QINGKE, color=fc.INK2, linewidth=0.9, linestyle=(0, (4, 3)),
            zorder=1)
axa.annotate("2000 m: candidate-definition threshold", xy=(2000, 3.0e4),
             xytext=(2030, 1.1e5), fontsize=7, color=fc.INK2, ha="left", va="top")

# ---------------------------------------------------------------- (b) Qingke histogram
bins = np.arange(2000, 6001, 250)
axb.hist(q.elev, bins=bins, color=fc.RED, edgecolor=fc.SURFACE, linewidth=0.7, zorder=2)
axb.axvline(q.elev.median(), color=fc.INK, linewidth=1.6, zorder=3)
axb.annotate(f"median {q.elev.median():.0f} m", xy=(q.elev.median(), 0.94),
             xycoords=("data", "axes fraction"),
             xytext=(6, 0), textcoords="offset points", fontsize=8, color=fc.INK,
             ha="left", va="top")
axb.set_xlabel("Elevation (m a.s.l.)")
axb.set_ylabel("Qingke candidates (count)")
axb.set_title(f"(b) Elevation distribution of the {len(q)} Qingke candidates",
              loc="left", pad=8)
axb.set_xlim(0, 6100)
axb.set_xticks(np.arange(0, 6001, 1000))
axb.grid(True, axis="y")

plt.tight_layout()
out = os.path.join(fc.OUT, "fig2_elevation_bands.png")
fig.savefig(out, dpi=300)
plt.close(fig)
print("wrote", out)
