"""Figure 1 - global barley occurrence records, coloured by elevation, with the
Qinghai-Tibet Plateau box and a bottom-right inset zoom showing the Qingke candidates.

Layout note: the inset is placed in a dedicated band *below* the map rather than
overlaying the map corner. A bottom-right overlay would cover Australia, which holds
1,743 of the 46,102 records; the only data-free bottom-right window in the map itself
(lon 60-115E, lat 60-25S, 1 record) is too small to read. `tight_layout` is not used
here because it is incompatible with manually positioned axes + inset_axes; a fixed
rect + `bbox_inches="tight"` gives a deterministic, unclipped layout instead.

Output: gis/processed/maps/fig1_sample_distribution.png  (300 dpi)
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Rectangle

import fig_common as fc

fc.use_style()

df = fc.load_matrix()
g, plateau, qingke = fc.subsets(df)

has_elev = g.elev.notna()
n_elev, n_noelev = int(has_elev.sum()), int((~has_elev).sum())
print(f"global n={len(g)}  with elev={n_elev}  no elev={n_noelev}")
print(f"plateau box n={len(plateau)}  qingke n={len(qingke)}")

# 360 deg lon x 145 deg lat at equal aspect -> map rect must be 2.483:1 or it letterboxes
FIG_W, FIG_H = 8.0, 5.6
MAP_H = 0.80 * FIG_W / ((180 - -180) / (85 - -60)) / FIG_H
INSET_W = 0.545
INSET_H = INSET_W * FIG_W / ((105 - 70) / (40 - 25)) / FIG_H
fig = plt.figure(figsize=(FIG_W, FIG_H))

# =============================== main global map ==============================
ax = fig.add_axes([0.065, 0.970 - MAP_H, 0.80, MAP_H])

sc = ax.scatter(g.loc[has_elev, "decimalLongitude"], g.loc[has_elev, "decimalLatitude"],
                c=g.loc[has_elev, "elev"], cmap=fc.ELEV_CMAP,
                vmin=0, vmax=5000, s=2.2, alpha=0.55, linewidths=0,
                rasterized=True, zorder=2)
# records with no elevation (WorldClim land-mask / ocean cells) - shown, not hidden
ax.scatter(g.loc[~has_elev, "decimalLongitude"], g.loc[~has_elev, "decimalLatitude"],
           s=2.2, c=fc.GREY, alpha=0.9, linewidths=0, rasterized=True, zorder=2)

lon0, lat0, lon1, lat1 = fc.PLATEAU
ax.add_patch(Rectangle((lon0, lat0), lon1 - lon0, lat1 - lat0, fill=False,
                       edgecolor=fc.INK, linewidth=1.4, zorder=5))
ax.annotate("Qinghai-Tibet Plateau box\n70-105E, 25-40N",
            xy=(lon1, lat1), xytext=(lon1 + 2, lat1 + 7),
            fontsize=8, color=fc.INK, ha="left", va="bottom", zorder=6,
            arrowprops=dict(arrowstyle="-", color=fc.INK, linewidth=0.9))

ax.set_xlim(-180, 180)
ax.set_ylim(-60, 85)
ax.set_xticks(np.arange(-180, 181, 60))
ax.set_yticks(np.arange(-60, 86, 30))
ax.set_xlabel("Longitude (degrees east)")
ax.set_ylabel("Latitude (degrees north)")
ax.set_aspect("equal", adjustable="box")
ax.grid(True)
ax.set_title("(a) Global barley (Hordeum vulgare) occurrence records by elevation",
             loc="left", pad=8)

cb = fig.colorbar(sc, ax=ax, fraction=0.028, pad=0.012, extend="max")
cb.set_label("Elevation (m a.s.l.)", fontsize=9)
cb.ax.tick_params(labelsize=8)
cb.outline.set_linewidth(0.6)

# identity is never colour-alone: one legend covering all mark classes
h = [plt.Line2D([], [], marker="o", ls="none", ms=4.5, mfc=plt.cm.viridis(0.72),
                mec="none", label=f"Records with elevation (n={n_elev:,})"),
     plt.Line2D([], [], marker="o", ls="none", ms=4.5, mfc=fc.GREY, mec="none",
                label=f"No elevation data (n={n_noelev:,})"),
     plt.Line2D([], [], marker="o", ls="none", ms=5.5, mfc=fc.RED, mec="white", mew=0.8,
                label=f"Qingke candidates (n={len(qingke)})")]
# legend lives in the free band under the map (left of the inset), not over the data
fig.legend(handles=h, loc="lower left", bbox_to_anchor=(0.065, 0.115),
           fontsize=8, labelspacing=0.5, title="Mark classes", title_fontsize=8)


# =============================== bottom-right inset ===========================
axi = fig.add_axes([0.400, 0.075, INSET_W, INSET_H])
axi.set_facecolor("#f7f7f4")
axi.plot([], [])
axi.scatter(plateau.decimalLongitude, plateau.decimalLatitude, s=7,
            c=fc.BLUE, alpha=0.8, linewidths=0, zorder=2)
axi.scatter(qingke.decimalLongitude, qingke.decimalLatitude, s=11, c=fc.RED,
            edgecolors="white", linewidths=0.35, zorder=3)
axi.set_xlim(lon0, lon1)
axi.set_ylim(lat0, lat1)
axi.set_xticks([70, 75, 80, 85, 90, 95, 100, 105])
axi.set_yticks([25, 30, 35, 40])
axi.tick_params(labelsize=7.5)
axi.grid(True)
axi.set_aspect("equal", adjustable="box")
for s in axi.spines.values():
    s.set_edgecolor(fc.INK)
    s.set_linewidth(0.9)
axi.set_xlabel("Longitude (degrees east)", fontsize=8)
axi.set_ylabel("Latitude (degrees north)", fontsize=8)
axi.set_title("(b) Inset: records inside the plateau box", fontsize=9.5,
              fontweight="bold", loc="left", pad=5)
hi = [plt.Line2D([], [], marker="o", ls="none", ms=5, mfc=fc.BLUE, mec="none",
                 label=f"Plateau-box records (n={len(plateau)})"),
      plt.Line2D([], [], marker="o", ls="none", ms=5, mfc=fc.RED, mec="white", mew=0.8,
                 label=f"Qingke candidates (n={len(qingke)})")]
axi.legend(handles=hi, loc="upper left", fontsize=7.5, labelspacing=0.3)

out = os.path.join(fc.OUT, "fig1_sample_distribution.png")
fig.savefig(out, dpi=300, bbox_inches="tight")
plt.close(fig)
print("wrote", out)
