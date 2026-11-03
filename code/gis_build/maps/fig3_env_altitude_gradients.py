"""Figure 3 - environmental gradients along elevation.

Three panels: mean annual temperature (wc_bio1), annual precipitation (wc_bio12) and
maximum temperature of the warmest month (wc_bio5) against elevation.

Marks: global records (recessive grey), plateau-box records (blue), Qingke candidates
(red). The solid line is the binned median of the plateau-box records, the dashed
lines the 10th/90th percentiles, the wash the 10-90% envelope.

Per-1000 m slopes are ordinary least squares fits over 2000-5000 m; precipitation is
fitted on log10(precipitation) so the slope converts to %/km. Slopes are reported both
for the plateau-box subset (annotated on the panels - the highland gradient of
interest) and for all global records (printed for the report).

Output: gis/processed/maps/fig3_env_altitude_gradients.png  (300 dpi)
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import stats

import fig_common as fc

fc.use_style()

PANELS = [("wc_bio1", "Mean annual temperature (degC)", "temperature", (0.03, 0.03)),
          ("wc_bio12", "Annual precipitation (mm)", "precipitation", (0.97, 0.97)),
          ("wc_bio5", "Max. temperature of warmest month (degC)", "temperature", (0.03, 0.03))]

df = fc.load_matrix(usecols=["decimalLongitude", "decimalLatitude", "occurrenceStatus", "elev",
                             "on_plateau", "is_qingke_candidate",
                             "wc_bio1", "wc_bio5", "wc_bio12"])
glob = df[df.elev.notna() & df.wc_bio1.notna()]
plat = glob[glob.on_plateau]
qk = glob[glob.is_qingke_candidate]
print(f"global with elev+climate n={len(glob)}  plateau box n={len(plat)}  qingke n={len(qk)}")

# ---- 2000-5000 m windows used for the slope fits ----------------------------
FIT_LO, FIT_HI = 2000.0, 5000.0


def fit_slope(sub, var, transform="linear"):
    """OLS of `var` (optionally log10-transformed) on elevation, 2000-5000 m."""
    s = sub[(sub.elev >= FIT_LO) & (sub.elev <= FIT_HI)][["elev", var]].dropna()
    s = s[s[var] > 0] if transform == "log" else s
    if len(s) < 10:
        return None
    y = np.log10(s[var].to_numpy()) if transform == "log" else s[var].to_numpy()
    x = s.elev.to_numpy()
    r = stats.linregress(x, y)
    if transform == "log":
        per_k = (10 ** (r.slope * 1000) - 1) * 100.0      # %/km (exact, not log-approx)
    else:
        per_k = r.slope * 1000.0                          # units/km
    return dict(n=len(s), per_km=per_k, se_per_km=(r.stderr * 1000.0 * 100.0
                                                  if transform == "log"
                                                  else r.stderr * 1000.0),
                r2=r.rvalue ** 2, p=r.pvalue, slope=r.slope)


rows = []
print("\n=== 2000-5000 m OLS slopes ===")
for var, label, kind, _ in PANELS:
    for name, sub, trans in (("plateau_box", plat, "log" if kind == "precipitation" else "linear"),
                             ("global", glob, "log" if kind == "precipitation" else "linear")):
        f = fit_slope(sub, var, trans)
        if f is None:
            continue
        unit = "%/km" if kind == "precipitation" else "degC/km"
        rows.append(dict(variable=var, subset=name, n=f["n"],
                         slope_per_km=round(f["per_km"], 2),
                         se=round(f["se_per_km"], 3),
                         ci_lo=round(f["per_km"] - 1.96 * f["se_per_km"], 2),
                         ci_hi=round(f["per_km"] + 1.96 * f["se_per_km"], 2),
                         r2=round(f["r2"], 3), unit=unit))
        print(f"{var:10s} {name:12s} n={f['n']:6d}  {f['per_km']:8.2f} {unit} "
              f"(se {f['se_per_km']:.3f})  R2={f['r2']:.3f}  p={f['p']:.2e}")
stats_tbl = pd.DataFrame(rows)

# ---- figure -----------------------------------------------------------------
fig, axes = plt.subplots(1, 3, figsize=(7.5, 3.5))
BIN = 250.0
edges = np.arange(0, 6000 + BIN, BIN)
centers = 0.5 * (edges[:-1] + edges[1:])

for ax, (var, ylabel, kind, anchor) in zip(axes, PANELS):
    # global records - recessive context
    ax.scatter(glob.elev, glob[var], s=1.8, c=fc.GREY, alpha=0.5, linewidths=0,
               rasterized=True, zorder=1)
    # plateau-box records
    ax.scatter(plat.elev, plat[var], s=5, c=fc.BLUE, alpha=0.55, linewidths=0,
               rasterized=True, zorder=2)
    # Qingke candidates
    ax.scatter(qk.elev, qk[var], s=8, c=fc.RED, alpha=0.9, edgecolors="white",
               linewidths=0.25, rasterized=True, zorder=3)

    # binned quantiles of the plateau-box records (pooled summary, ink not series hue)
    idx = np.digitize(plat.elev, edges) - 1
    keep = (idx >= 0) & (idx < len(centers))
    tmp = pd.DataFrame({"b": idx[keep], "v": plat[var].to_numpy()[keep]})
    gb = tmp.groupby("b")["v"]
    med, lo, hi, n = gb.median(), gb.quantile(0.10), gb.quantile(0.90), gb.size()
    ok = n >= 15
    if ok.any():
        x = centers[med.index[ok]]
        ax.fill_between(x, lo[ok], hi[ok], color=fc.INK, alpha=0.07, zorder=4,
                        linewidth=0)
        ax.plot(x, hi[ok], color=fc.INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=5)
        ax.plot(x, lo[ok], color=fc.INK2, linewidth=1.0, linestyle=(0, (4, 3)), zorder=5)
        ax.plot(x, med[ok], color=fc.INK, linewidth=1.8, zorder=6)

    # slope annotation (plateau-box subset = the gradient the panel is about)
    f = fit_slope(plat, var, "log" if kind == "precipitation" else "linear")
    unit = "%/km" if kind == "precipitation" else "degC/km"
    txt = (f"{f['per_km']:+.2f} {unit}\n"
           f"(2000-5000 m, R2={f['r2']:.2f}, n={f['n']})")
    ax.annotate(txt, xy=anchor, xycoords="axes fraction", fontsize=7.5,
                color=fc.INK, ha="left" if anchor[0] < 0.5 else "right",
                va="bottom" if anchor[1] < 0.5 else "top",
                bbox=dict(boxstyle="round,pad=0.32", fc="white", ec=fc.GRID, lw=0.6))

    ax.set_xlabel("Elevation (m a.s.l.)")
    ax.set_ylabel(ylabel)
    ax.set_xlim(0, 6000)
    ax.set_xticks(np.arange(0, 6001, 1500))
    ax.grid(True)

axes[0].set_title("(a) Mean annual temperature", loc="left", pad=6)
axes[1].set_title("(b) Annual precipitation", loc="left", pad=6)
axes[2].set_title("(c) Warmest-month maximum", loc="left", pad=6)
# Panel B focuses the precipitation gradient on the high-elevation range;
# panel A (and panel C) retain the full low-to-high elevation context.
axes[1].set_xlim(1800, 6000)
axes[1].set_xticks(np.arange(1800, 6001, 1000))
axes[1].set_ylim(bottom=0)

h = [plt.Line2D([], [], marker="o", ls="none", ms=3.5, mfc=fc.GREY, mec="none",
                label=f"All records (n={len(glob):,})"),
     plt.Line2D([], [], marker="o", ls="none", ms=4, mfc=fc.BLUE, mec="none",
                label=f"Plateau-box records (n={len(plat):,})"),
     plt.Line2D([], [], marker="o", ls="none", ms=4.5, mfc=fc.RED, mec="white", mew=0.7,
                label=f"Qingke candidates (n={len(qk)})"),
     plt.Line2D([], [], color=fc.INK, lw=1.8, label="Median (plateau box)"),
     plt.Line2D([], [], color=fc.INK2, lw=1.0, ls=(0, (4, 3)), label="10th / 90th pct")]
fig.legend(handles=h, loc="lower center", ncol=5, fontsize=7.5,
           bbox_to_anchor=(0.5, -0.055), columnspacing=1.2, handletextpad=0.4)

plt.tight_layout()
out = os.path.join(fc.OUT, "fig3_env_altitude_gradients.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
fig.savefig(os.path.join(fc.OUT, "fig3_env_altitude_gradients.pdf"), bbox_inches="tight")
plt.close(fig)
print("\nwrote", out)
stats_tbl.to_csv(os.path.join(fc.OUT, "fig3_altitude_slopes.tsv"), sep="\t", index=False)
print("wrote", os.path.join(fc.OUT, "fig3_altitude_slopes.tsv"))
