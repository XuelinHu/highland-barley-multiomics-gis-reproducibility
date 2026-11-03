"""Figure 4 - plateau-region environmental surfaces (mean annual temperature, annual
precipitation, elevation) with the 558 Qingke candidates overlaid.

Rasters are read windowed and decimated straight from the global WorldClim 2.1 30
arc-second GeoTIFFs (43200 x 21600 px) with rasterio: only the 70-105E / 25-40N
window is touched, area-averaged onto a 0.05 deg grid (700 x 300). Peak memory stays
well under 100 MB of float32 instead of the ~3.7 GB a full-array read would need.

WorldClim nodata is -3.4e38 (float) / -32768 (int16); the land mask is taken from
wc_bio1 and applied to every layer so coastal averages do not bleed ocean into land.

Output: gis/processed/maps/fig4_env_maps.png  (300 dpi)
"""

from __future__ import annotations

import os

import matplotlib.pyplot as plt
import numpy as np
from matplotlib import colors as mcolors

import fig_common as fc

fc.use_style()

NROWS, NCOLS = 300, 700            # 0.05 deg grid over 25 rows of latitude / 70 of longitude
SHAPE = (NROWS, NCOLS)
P = fc.PLATEAU

df = fc.load_matrix(usecols=["decimalLongitude", "decimalLatitude", "occurrenceStatus", "is_qingke_candidate"])
qk = df[df.is_qingke_candidate]
print(f"Qingke candidates overlaid: {len(qk)}")

# ---- read the three layers, decimated to the plotting grid -------------------
land, extent = fc.read_window(os.path.join(fc.WC, "wc2.1_30s_bio_1.tif"), P, SHAPE)
land_ok = np.isfinite(land) & (np.abs(land) < 1e30)
print("land cells in window:", int(land_ok.sum()), "/", land_ok.size)

# Temperature is diverging around the physically meaningful 0 degC midpoint; the
# precipitation field is strongly right-skewed (12-8,485 mm) so it is shown on a
# power-law stretch; elevation uses the perceptual-uniform viridis ramp.
TEMP_NORM = mcolors.TwoSlopeNorm(vmin=-20, vcenter=0, vmax=30)
# precipitation spans 12-8,485 mm (right-skewed by ~2 orders of magnitude), so a log
# scale is the scale on which the plateau's 200-800 mm range is actually legible
PRECIP_NORM = mcolors.LogNorm(vmin=10, vmax=8500)

layers = []
for fname, cmap, label, unit, norm in [
        ("wc2.1_30s_bio_1.tif", "RdBu_r", "Mean annual temperature", "degC", TEMP_NORM),
        ("wc2.1_30s_bio_12.tif", "YlGnBu", "Annual precipitation", "mm", PRECIP_NORM),
        ("wc2.1_30s_elev.tif", fc.ELEV_CMAP, "Elevation", "m a.s.l.", None)]:
    arr, _ = fc.read_window(os.path.join(fc.WC, fname), P, SHAPE)
    arr = np.where(np.abs(arr) > 1e30, np.nan, arr)
    if fname.endswith("elev.tif"):
        arr = np.where(arr < -100, np.nan, arr)     # no sub-sea-level land in window
    layers.append((fc.mask_nodata(arr, land_mask=land_ok), cmap, label, unit, norm))
    print(f"{fname}: min={np.nanmin(arr):.1f} max={np.nanmax(arr):.1f}")

# ---- figure -----------------------------------------------------------------
fig, axes = plt.subplots(3, 1, figsize=(7.2, 10.0))
for i, (ax, (arr, cmap, label, unit, norm), tag) in enumerate(zip(axes, layers, "abc")):
    im = ax.imshow(arr, extent=extent, origin="upper", cmap=cmap, norm=norm,
                   interpolation="nearest", aspect="equal")
    im.cmap.set_bad("#f0efec")
    ax.scatter(qk.decimalLongitude, qk.decimalLatitude, s=4.5, c=fc.INK,
               linewidths=0, alpha=0.85, zorder=3,
               label=f"Qingke candidates (n={len(qk)})" if i == 0 else None)
    ax.set_title(f"({tag}) {label} ({unit})", loc="left", pad=6)
    ax.set_xlabel("Longitude (degrees east)")
    ax.set_ylabel("Latitude (degrees north)")
    ax.set_xticks(np.arange(70, 106, 5))
    ax.set_yticks(np.arange(25, 41, 5))
    ax.set_xlim(P[0], P[2])
    ax.set_ylim(P[1], P[3])
    ax.grid(False)
    if i == 0:
        ax.legend(loc="lower left", fontsize=7.5)
    cb = fig.colorbar(im, ax=ax, fraction=0.030, pad=0.012,
                      extend="both" if norm is PRECIP_NORM else "neither")
    cb.set_label(f"{label} ({unit})" + ("  [log scale]" if norm is PRECIP_NORM
                                        else ""), fontsize=8)
    cb.ax.tick_params(labelsize=7.5)
    cb.outline.set_linewidth(0.5)

plt.tight_layout()
out = os.path.join(fc.OUT, "fig4_env_maps.png")
fig.savefig(out, dpi=300)
plt.close(fig)
print("wrote", out)
