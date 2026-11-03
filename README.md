---
license: other
license_name: mixed-source-terms
license_link: https://github.com/XuelinHu/highland-barley-multiomics-gis-reproducibility/blob/main/SOURCES.md
language:
- en
- zh
tags:
- biology
- barley
- reproducibility
pretty_name: Barley drought candidates and high-elevation occurrence environments
---

# Highland barley reproducibility archive

Public companion to **Drought-responsive barley protein candidates and high-elevation occurrence environments: parallel analyses of public data**. Molecular samples and GIS occurrences are unpaired. The archive does not establish gene–environment associations, adaptive function or habitat suitability.

GitHub: https://github.com/XuelinHu/highland-barley-multiomics-gis-reproducibility

Hugging Face: https://huggingface.co/datasets/xuelinhu/highland-barley-multiomics-gis-reproducibility

Version: v1.0.1. Exact files: SHA256SUMS.json. Original source base commit and subsequent revision scope: release_metadata.json. No archive DOI is claimed unless one is registered.

## Contents

- `evidence/processed/reanalysis/`: drought normalized expression, original fitted R objects, contrast outputs, design tables, probe/gene consensus, candidate mapping, historical fragment annotations and selected structures/network inputs.
- `evidence/processed/mapping/`: platform-to-MorexV3 mappings.
- `evidence/gis/processed/`: retained GBIF points (no climate values), saved spatial blocks and historical fold/summary metrics.
- `results/`: 240 audited fold rows, assignments, panel manifest and full-length FASTA, corrected cold counts/permutations, WGCNA memberships, slope intervals, spatial/climate sensitivities, per-record predictions and vector figures.
- `code/`: portable historical workflows plus reviewer-requested analyses. Python/R/shell paths use `BARLEY_DATA_ROOT` (default `evidence`). Historical baseline tuning differs from the corrected sensitivity workflow; see the configuration JSON.
- `sources/` and `SOURCES.md`: precise original access/query information, versions, licences and download scripts.

## Verify the release without external data

Use an existing environment containing numpy, pandas, scipy and scikit-learn. The authors used Conda `pyg` locally and the recorded R/Linux environments for the source runs; version records are in results/ and evidence/. No installation is needed if these packages already exist.

```sh
conda run -n pyg python code/revision/verify_public_release.py --root .
```

This verifies every archive checksum, recomputes sensitivity AP/ROC-AUC from saved held-out predictions, checks fold class counts and disjoint partitions, and checks the corrected cold-count and sequence/model invariants. It does not masquerade as rerunning unavailable external services.

## Reproduce the molecular and sensitivity analyses

R requires limma and WGCNA; historical versions are retained in sessionInfo records. Use an existing appropriate R installation.

```sh
Rscript code/revision/export_molecular_inputs.R evidence rerun/molecular
Rscript code/revision/module_membership_sensitivity.R evidence rerun/wgcna
```

For source-data-dependent checks, first retrieve the licensed climate/terrain/soil inputs using sources/download_scripts and the original report. Set `BARLEY_DATA_ROOT` to the absolute evidence directory. The point snapshot is fixed; a new GBIF query is unnecessary. `code/gis_build/build_rasters.py` assembles terrain/soil rasters (requires GDAL), and `code/gis_build/extract_env.py` samples them at those points. Raster acquisition is approximately 33 GB. Do not publish the rebuilt WorldClim point matrix without permission. Paths in the distributed Python/R/shell source have been adapted from the original machine root; those edits do not alter analytical rules.

```sh
python code/gis_build/build_rasters.py
python code/gis_build/extract_env.py
python code/revision/complete_revision_analysis.py --evidence evidence --output rerun/completed
python code/revision/spatial_sensitivity.py --evidence evidence --output rerun/spatial
python code/revision/redraw_revision_figures.py --evidence evidence --results results --output rerun/figures
```

The full audit requires the newly reconstructed climate matrix. Raw-array preprocessing uses `code/molecular/reanalyse_drought.R` and additionally requires the GEO raw arrays/phenotype metadata in its documented paths. The supplied fitted R objects and normalized matrices support downstream reproduction without raw-array downloads. The 107 historical protein-annotation inputs were only 60-residue prefixes: use `results/top100_proteins_full_length.fasta` for a future full-length run, not `legacy_60aa_annotation_input.fasta`. An EBI full-length submission was rejected under a temporary fair-use restriction, so full-length annotation/disorder is an explicitly uncompleted external step.

## Corrections and interpretation

3,076 direction-consistent probes include 1,568 uniquely mapped probes representing 1,197 genes, 54 ambiguous probes and 1,454 without a qualified hit. GSE41513 is 28/33 gene–contrast pairs involving 24 genes, not 49/61 unique genes; its four fixed-margin p-values equal 1. GSE47882 retains 71.7–84.1% concordance, subject to exchangeability/co-regulation limitations. DHN1 is corrected to DHN2. The ten-protein panel is purposive with retrospective inclusion rationales, not an unbiased or prespecified top-ten selection. No wet-lab validation or geographically matched molecular samples are added. High spatial AP classifies an elevation-defined geographic proxy among sampled occurrences.

## Licensing and citation

Code is MIT. Project-derived tables/maps are CC BY 4.0 subject to source terms. Third-party observations, references and models retain their original licences; mixed GBIF records include noncommercial restrictions. See SOURCES.md, record-level licence columns and CITATION.cff. Private manuscript/reviewer correspondence and private Git history are excluded.
