"""Shared loading / styling helpers for the TODO stage-6 paper figures.

Data (read-only)  : /ds2/workspace/ai/highland-barley-omics-gis/gis/processed/
Raster inputs     : /ds2/workspace/ai/highland-barley-omics-gis/gis/raw/{worldclim,chelsa}/
Figure output     : /ds2/workspace/ai/highland-barley-omics-gis/gis/processed/maps/

Palette: categorical slots 1 (blue #2a78d6) and 8 (red #e34948) from the validated
reference palette - the blue/red pair passes the CVD (dE 21.6 protan) and
normal-vision (dE 32.3) all-pairs gates with >=3:1 contrast on the light surface.
Elevation uses the perceptual-uniform `viridis` ramp (sequential magnitude).
Grey marks are deliberately recessive *context*, not an identity series.
"""

from __future__ import annotations

import os

import matplotlib as mpl
import numpy as np
import pandas as pd

# ---------------------------------------------------------------- paths
PROC = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed')
RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
OUT = os.path.join(PROC, "maps")
MATRIX = os.path.join(PROC, "barley_env_matrix.tsv.gz")
QINGKE = os.path.join(PROC, "qingke_plateau_subset.tsv.gz")

WC = os.path.join(RAW, "worldclim")          # wc2.1_30s_*.tif  (global, 30 arc-sec)
CHELSA = os.path.join(RAW, "chelsa")         # CHELSA_bio*_1981-2010_V.2.1.tif

# ---------------------------------------------------------------- study region
PLATEAU = (70.0, 25.0, 105.0, 40.0)          # lon_min, lat_min, lon_max, lat_max
PLATEAU_LABEL = "Qinghai-Tibet Plateau box (70-105E, 25-40N)"
ELEV_MIN_QINGKE = 2000                        # candidate definition threshold (m)

# ---------------------------------------------------------------- palette
BLUE = "#2a78d6"        # categorical slot 1 - plateau-box records
RED = "#e34948"         # categorical slot 8 - Qingke candidates
GREY = "#c9c8c2"        # recessive context - global records
INK = "#0b0b0b"         # text / primary
INK2 = "#52514e"        # secondary text
GRID = "#d8d7d2"        # hairline grid
SURFACE = "#ffffff"

ELEV_CMAP = "viridis"
PROB_CMAP = "YlOrRd"    # warm = high suitability

CAPTION_STYLE = dict(fontsize=9, color=INK2)

# ---------------------------------------------------------------- style


def use_style() -> None:
    """Publication rcParams: 300 dpi, sans, recessive grid, no chartjunk."""
    mpl.rcParams.update({
        "figure.dpi": 300,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
        "savefig.facecolor": SURFACE,
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "font.family": "sans-serif",
        "font.sans-serif": ["Liberation Sans", "Nimbus Sans", "DejaVu Sans"],
        "font.size": 9,
        "axes.titlesize": 10,
        "axes.labelsize": 9,
        "axes.titleweight": "bold",
        "axes.edgecolor": INK2,
        "axes.linewidth": 0.7,
        "axes.labelcolor": INK,
        "axes.grid": True,
        "axes.axisbelow": True,
        "grid.color": GRID,
        "grid.linewidth": 0.6,
        "grid.linestyle": "-",
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.labelsize": 8,
        "ytick.labelsize": 8,
        "xtick.direction": "out",
        "ytick.direction": "out",
        "legend.frameon": False,
        "legend.fontsize": 8,
        "legend.handletextpad": 0.5,
        "legend.borderaxespad": 0.3,
        "lines.linewidth": 2.0,
        "lines.solid_capstyle": "round",
    })


def panel_tag(ax, letter: str) -> None:
    """Bold (a)/(b)/(c) tag in the top-left corner, outside the data."""
    ax.text(0.0, 1.045, f"({letter})", transform=ax.transAxes,
            fontsize=10, fontweight="bold", color=INK, va="bottom", ha="left")


# ---------------------------------------------------------------- data


MATRIX_COLS = [
    "decimalLongitude", "decimalLatitude", "year", "basisOfRecord", "countryCode", "occurrenceStatus",
    "on_plateau", "is_qingke_candidate", "elev", "elev_source", "slope", "aspect",
    "wc_bio1", "wc_bio5", "wc_bio6", "wc_bio12", "wc_bio15",
]


def load_matrix(cols=None, usecols=None) -> pd.DataFrame:
    """Load the cleaned 46,102-point environment matrix (only the needed columns)."""
    df = pd.read_csv(MATRIX, sep="\t", usecols=usecols or cols or MATRIX_COLS,
                     low_memory=False)
    if "occurrenceStatus" in df.columns:
        df = df[df["occurrenceStatus"].fillna("PRESENT").eq("PRESENT")].copy()
    for flag in ("is_qingke_candidate", "on_plateau"):
        if flag in df.columns:
            df[flag] = df[flag].fillna(False).astype(bool)
    return df


def subsets(df: pd.DataFrame):
    """(global, plateau_box, qingke) row masks used by every figure."""
    qingke = df["is_qingke_candidate"].to_numpy()
    plateau = df["on_plateau"].to_numpy()
    return df, df[plateau], df[qingke]


def load_qingke() -> pd.DataFrame:
    return pd.read_csv(QINGKE, sep="\t", low_memory=False)


# ---------------------------------------------------------------- rasters


def read_window(path, bounds=PLATEAU, out_shape=None, resampling=None):
    """Memory-safe windowed read, optionally decimated to `out_shape` (rows, cols).

    Returns (array, extent) where extent = (lon_min, lon_max, lat_min, lat_max)
    exactly matching the window, so `imshow(..., extent=extent, origin="upper")`
    georeferences correctly. Decimating with an area-average resampling keeps the
    output grid aligned to the window bounds (cell i spans window_min + i*step).
    """
    import rasterio
    from rasterio.enums import Resampling
    from rasterio.windows import bounds as win_bounds
    from rasterio.windows import from_bounds

    if resampling is None:
        resampling = Resampling.average
    with rasterio.open(path) as src:
        win = from_bounds(*bounds, transform=src.transform)
        arr = src.read(1, window=win, out_shape=out_shape, resampling=resampling)
        b = win_bounds(win, transform=src.transform)
    extent = (b[0], b[2], b[1], b[3])       # lon_min, lon_max, lat_min, lat_max
    return arr.astype("float32"), extent


def mask_nodata(arr: np.ndarray, land_mask: np.ndarray | None = None) -> np.ma.MaskedArray:
    """Mask WorldClim (|v|>1e30) / SRTM (-32768) nodata, plus optional land mask."""
    bad = ~np.isfinite(arr) | (np.abs(arr) > 1e30) | (arr <= -30000.0)
    if land_mask is not None:
        bad |= ~land_mask
    return np.ma.masked_array(arr, mask=bad)


def lons_lats(bounds=PLATEAU, ncols=700, nrows=300):
    """Cell-centre coordinate vectors for a bounds-aligned regular grid."""
    lon_min, lat_min, lon_max, lat_max = bounds
    lon = lon_min + (np.arange(ncols) + 0.5) * (lon_max - lon_min) / ncols
    # raster row 0 is the northernmost band -> descending latitudes
    lat = lat_max - (np.arange(nrows) + 0.5) * (lat_max - lat_min) / nrows
    return lon, lat
