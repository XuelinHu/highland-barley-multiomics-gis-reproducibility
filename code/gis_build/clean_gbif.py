#!/usr/bin/env python3
"""Clean the GBIF Hordeum vulgare occurrence download.

Every step is counted and written to processed/cleaning_counts.json.
Output: processed/points_clean.tsv.gz  (cleaned points, no environment yet)
"""
import json
import os
import time

import numpy as np
import pandas as pd

RAW = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/raw')
PROC = __import__('os').path.join(__import__('os').environ.get('BARLEY_DATA_ROOT', 'evidence'), 'gis/processed')
IN = os.path.join(RAW, "gbif", "gbif_hordeum_vulgare.tsv")
OUT = os.path.join(PROC, "points_clean.tsv.gz")

KEEP_BASIS = [
    "PRESERVED_SPECIMEN",
    "HUMAN_OBSERVATION",
    "OBSERVATION",
    "MATERIAL_SAMPLE",
    "LIVING_SPECIMEN",
    "MACHINE_OBSERVATION",
]

META_COLS = [
    "key", "gbifID", "scientificName", "species", "taxonKey",
    "decimalLatitude", "decimalLongitude", "coordinateUncertaintyInMeters",
    "year", "month", "day", "basisOfRecord", "occurrenceStatus",
    "countryCode", "country", "stateProvince", "locality",
    "elevation", "elevationAccuracy", "recordedBy",
    "institutionCode", "collectionCode", "catalogNumber", "datasetKey",
    "publishingCountry", "establishmentMeans", "occurrenceID", "issues", "license",
]

counts = []


def step(name, n_before, n_after, extra=""):
    rec = {"step": name, "before": int(n_before), "after": int(n_after),
           "removed": int(n_before - n_after), "note": extra}
    counts.append(rec)
    print(f"{name:<42} {n_before:>7} -> {n_after:>7}  (-{n_before-n_after})  {extra}", flush=True)
    return n_after


def main():
    t0 = time.time()
    df = pd.read_csv(IN, sep="\t", dtype=str, keep_default_na=False, na_values=[""],
                     low_memory=False)
    n0 = len(df)
    print(f"raw records: {n0}", flush=True)
    df["decimalLatitude"] = pd.to_numeric(df["decimalLatitude"], errors="coerce")
    df["decimalLongitude"] = pd.to_numeric(df["decimalLongitude"], errors="coerce")
    df["year_num"] = pd.to_numeric(df["year"], errors="coerce")
    df["elevation_recorded"] = pd.to_numeric(df["elevation"], errors="coerce")
    df["coord_uncertainty_m"] = pd.to_numeric(
        df["coordinateUncertaintyInMeters"], errors="coerce")

    # --- step 1: coordinates ------------------------------------------------
    n = len(df)
    bad_missing = df.decimalLatitude.isna() | df.decimalLongitude.isna()
    bad_zero = (df.decimalLatitude == 0) | (df.decimalLongitude == 0)
    bad_range = ~(df.decimalLatitude.between(-90, 90)
                  & df.decimalLongitude.between(-180, 180))
    keep = ~(bad_missing | bad_zero | bad_range)
    df = df[keep].copy()
    n = step("1. coords valid (non-missing, non-zero, in range)", n, len(df),
             f"missing={int(bad_missing.sum())} zero={int(bad_zero.sum())} "
             f"out_of_range={int(bad_range.sum())}")

    # --- step 2: deduplicate on 4-decimal rounded coordinates ---------------
    # keep the most precise record (lowest coordinateUncertainty) at each
    # rounded coordinate, tie-broken by the most recent year
    df["lon_r"] = df.decimalLongitude.round(4)
    df["lat_r"] = df.decimalLatitude.round(4)
    df = df.sort_values(
        by=["coord_uncertainty_m", "year_num"],
        ascending=[True, False], na_position="last", kind="stable")
    df = df.drop_duplicates(subset=["lon_r", "lat_r"], keep="first").copy()
    n = step("2. unique coordinates (rounded to 4 decimals)", n, len(df),
             "kept most precise record per coordinate")

    # --- step 3: year -------------------------------------------------------
    n_missing_year = int(df.year_num.isna().sum())
    n_pre1950 = int((df.year_num < 1950).sum())
    bad_year = df.year_num.isna() | (df.year_num < 1950)
    df = df[~bad_year].copy()
    n = step("3. year present and >= 1950", n, len(df),
             f"missing_year={n_missing_year} before_1950={n_pre1950}")

    # --- step 4: basisOfRecord ---------------------------------------------
    basis_all = df.basisOfRecord.value_counts(dropna=False).to_dict()
    dropped_basis = {k: int(v) for k, v in basis_all.items() if k not in KEEP_BASIS}
    keep = df.basisOfRecord.isin(KEEP_BASIS)
    df = df[keep].copy()
    n = step("4. basisOfRecord in allowed set", n, len(df),
             f"allowed={KEEP_BASIS} dropped={dropped_basis}")

    # --- assemble output ----------------------------------------------------
    df["year"] = df.year_num.astype("Int64")
    out = pd.DataFrame({
        "gbif_key": df["key"].astype(str),
        "gbifID": df["gbifID"].astype(str),
        "scientificName": df["scientificName"],
        "lon": df["lon_r"].astype(float),
        "lat": df["lat_r"].astype(float),
        "year": df["year"],
        "month": pd.to_numeric(df["month"], errors="coerce").astype("Int64"),
        "day": pd.to_numeric(df["day"], errors="coerce").astype("Int64"),
        "basisOfRecord": df["basisOfRecord"],
        "occurrenceStatus": df["occurrenceStatus"],
        "countryCode": df["countryCode"],
        "country": df["country"],
        "stateProvince": df["stateProvince"],
        "locality": df["locality"],
        "elevation_recorded_m": df["elevation_recorded"],
        "elevation_accuracy_m": pd.to_numeric(df["elevationAccuracy"], errors="coerce"),
        "coord_uncertainty_m": df["coord_uncertainty_m"],
        "recordedBy": df["recordedBy"],
        "institutionCode": df["institutionCode"],
        "collectionCode": df["collectionCode"],
        "catalogNumber": df["catalogNumber"],
        "datasetKey": df["datasetKey"],
        "publishingCountry": df["publishingCountry"],
        "establishmentMeans": df["establishmentMeans"],
        "occurrenceID": df["occurrenceID"],
        "issues": df["issues"],
        "license": df["license"],
    })
    os.makedirs(PROC, exist_ok=True)
    out.to_csv(OUT, sep="\t", index=False, compression="gzip")
    print(f"wrote {OUT}  rows={len(out)}", flush=True)

    # --- diagnostics --------------------------------------------------------
    stats = {
        "raw_records": int(n0),
        "clean_records": int(len(out)),
        "steps": counts,
        "basis_of_record_raw": {str(k): int(v) for k, v in basis_all.items()},
        "basis_of_record_dropped": dropped_basis,
        "allowed_basis": KEEP_BASIS,
        "occurrence_status": {str(k): int(v)
                              for k, v in out.occurrenceStatus.value_counts().items()},
        "establishment_means": {str(k): int(v) for k, v in
                                out.establishmentMeans.fillna("NA").value_counts().items()},
        "year_range": [int(out.year.min()), int(out.year.max())],
        "year_median": int(out.year.median()),
        "year_hist": {str(k): int(v) for k, v in
                      (out.year // 10 * 10).value_counts().sort_index().items()},
        "top_countries": {str(k): int(v) for k, v in
                          out.country.value_counts().head(25).items()},
        "scientificName_top": {str(k): int(v) for k, v in
                               out.scientificName.value_counts().head(15).items()},
        "lon_min": float(out.lon.min()), "lon_max": float(out.lon.max()),
        "lat_min": float(out.lat.min()), "lat_max": float(out.lat.max()),
        "n_on_plateau_box": int(((out.lon.between(70, 105)) & (out.lat.between(25, 40))).sum()),
    }
    with open(os.path.join(PROC, "cleaning_counts.json"), "w") as fh:
        json.dump(stats, fh, indent=1, ensure_ascii=False)
    print(json.dumps({k: v for k, v in stats.items()
                      if k not in ("steps", "top_countries", "scientificName_top",
                                   "basis_of_record_raw", "year_hist")},
                     indent=1, ensure_ascii=False), flush=True)
    print(f"done in {time.time()-t0:.1f}s", flush=True)


if __name__ == "__main__":
    main()
