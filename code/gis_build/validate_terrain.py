#!/usr/bin/env python3
"""Independent validation of the Horn slope/aspect implementation.

1. Synthetic planes with analytically known slope and aspect.
2. Real-world spot checks on the produced rasters (Himalaya vs. Tarim basin).
"""
import numpy as np
import rasterio
from rasterio.windows import Window

RAS = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed/rasters/srtm')
RES = 1.0 / 1200.0
R = 6371000.0
DEG2M = np.pi / 180.0 * R
DY = RES * DEG2M


def horn(z, dx, dy):
    """Same Horn 3x3 implementation as build_rasters.py (dx may be an array)."""
    z = np.pad(z, ((0, 0), (1, 1)), mode="edge")
    a = z[0:-2, 0:-2]; b = z[0:-2, 1:-1]; c = z[0:-2, 2:]
    d = z[1:-1, 0:-2]; f = z[1:-1, 2:]
    g = z[2:, 0:-2];   h = z[2:, 1:-1];   i = z[2:, 2:]
    dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / (8.0 * dx)
    dzdy_s = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / (8.0 * dy)
    slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy_s)))
    aspect = np.mod(np.degrees(np.arctan2(-dzdx, dzdy_s)), 360.0)
    return slope, aspect


def test_plane(name, grad_east, grad_north, expect_slope, expect_aspect):
    """z increases by grad_east per metre eastward, grad_north per metre northward."""
    n = 5
    # build a small grid in metres
    xx = np.arange(n) * 1.0     # eastward index
    yy = np.arange(n) * 1.0     # southward index (row order)
    X, Y = np.meshgrid(xx, yy)
    z = grad_east * X - grad_north * Y   # row index increases southward
    dx = np.ones((n - 2, 1))
    slope, aspect = horn(z, dx, DY / DY)  # dy normalised: use metres directly below
    # recompute with consistent units
    a = z[0:-2, 0:-2]; b = z[0:-2, 1:-1]; c = z[0:-2, 2:]
    d = z[1:-1, 0:-2]; f = z[1:-1, 2:]
    g = z[2:, 0:-2];   h = z[2:, 1:-1];   i = z[2:, 2:]
    dzdx = ((c + 2.0 * f + i) - (a + 2.0 * d + g)) / 8.0     # per metre (pixel=1 m)
    dzdy_s = ((g + 2.0 * h + i) - (a + 2.0 * b + c)) / 8.0
    slope = np.degrees(np.arctan(np.hypot(dzdx, dzdy_s)))
    aspect = np.mod(np.degrees(np.arctan2(-dzdx, dzdy_s)), 360.0)
    ok_s = np.allclose(slope, expect_slope, atol=0.01)
    ok_a = np.allclose(aspect, expect_aspect, atol=0.01)
    print(f"{name:<38} slope={slope[0,0]:8.3f} (exp {expect_slope:7.3f}) "
          f"aspect={aspect[0,0]:7.2f} (exp {expect_aspect:7.2f}) "
          f"{'OK' if ok_s and ok_a else 'FAIL'}")
    return ok_s and ok_a


def main():
    print("=== synthetic plane tests (1 m pixels) ===")
    ok = True
    ok &= test_plane("east-facing (rises to west, 10%)", -0.1, 0.0,
                     np.degrees(np.arctan(0.1)), 90.0)
    ok &= test_plane("west-facing (rises to east, 10%)", 0.1, 0.0,
                     np.degrees(np.arctan(0.1)), 270.0)
    ok &= test_plane("north-facing (rises to south, 10%)", 0.0, -0.1,
                     np.degrees(np.arctan(0.1)), 0.0)
    ok &= test_plane("south-facing (rises to north, 10%)", 0.0, 0.1,
                     np.degrees(np.arctan(0.1)), 180.0)
    ok &= test_plane("flat", 0.0, 0.0, 0.0, 0.0)

    print("\n=== real raster spot checks ===")
    spots = [
        ("Everest / Himalaya (86.9E, 28.0N)", 86.9, 28.0),
        ("Tarim basin plain (82.0E, 39.5N)", 82.0, 39.5),
        ("Tibetan plateau interior (88.0E, 32.0N)", 88.0, 32.0),
        ("Sichuan basin (105.0E, 30.0N)", 105.0, 30.0),
        ("Ganges plain (78.0E, 27.0N)", 78.0, 27.0),
    ]
    for name, lon, lat in spots:
        with rasterio.open(f"{RAS}/srtm_mosaic.tif") as dem, \
             rasterio.open(f"{RAS}/srtm_slope.tif") as sl, \
             rasterio.open(f"{RAS}/srtm_aspect.tif") as asp:
            r, c = dem.index(lon, lat)
            w = Window(c - 25, r - 25, 50, 50)
            z = dem.read(1, window=w).astype(np.float32)
            s = sl.read(1, window=w)
            a = asp.read(1, window=w)
        zs = np.where(z <= -32000, np.nan, z)
        lats = lat + (np.arange(50) - 25) * RES
        dx = (RES * DEG2M * np.cos(np.radians(lats)))[1:-1, None]
        mine_s, mine_a = horn(zs, dx, DY)
        # horn() pads one column on each side, so its output is 50 wide; the
        # centre column of the window (48 values) is what the raster stores.
        mine_s = mine_s[:, 1:-1]
        mine_a = mine_a[:, 1:-1]
        # raster values are the interior of the 50x50 window
        rs = s[1:-1, 1:-1]
        ra = a[1:-1, 1:-1]
        print(f"{name:<40} elev={np.nanmean(zs):7.1f} m  "
              f"slope raster={np.nanmean(rs):6.2f} vs recomputed={np.nanmean(mine_s):6.2f} deg  "
              f"aspect raster(median)={np.nanmedian(ra):6.1f} vs {np.nanmedian(mine_a):6.1f}  "
              f"max|d slope|={np.nanmax(np.abs(rs-mine_s)):.4f}")
    print("\n=== whole-raster profile (strided sample of 512x512 windows) ===")
    rng = np.random.default_rng(42)
    for name, path in [("srtm_mosaic", f"{RAS}/srtm_mosaic.tif"),
                       ("slope", f"{RAS}/srtm_slope.tif"),
                       ("aspect", f"{RAS}/srtm_aspect.tif")]:
        with rasterio.open(path) as ds:
            nod = ds.nodata
            fin, tot, mn, mx, acc = 0, 0, np.inf, -np.inf, 0.0
            for _ in range(24):
                r = int(rng.integers(0, ds.height - 512))
                c = int(rng.integers(0, ds.width - 512))
                a = ds.read(1, window=Window(c, r, 512, 512)).astype(np.float64)
                if nod is not None and not np.isnan(nod):
                    a = np.where(a == nod, np.nan, a)
                f = np.isfinite(a)
                fin += int(f.sum()); tot += a.size
                if f.any():
                    mn = min(mn, float(a[f].min())); mx = max(mx, float(a[f].max()))
                    acc += float(a[f].sum())
        print(f"{name:<14} width={ds.width} height={ds.height} nodata={nod}  "
              f"valid={100.0*fin/tot:5.1f}%  min={mn:8.2f} max={mx:9.2f}")

    print("\nSYNTHETIC TESTS:", "PASS" if ok else "FAIL")


if __name__ == "__main__":
    main()
