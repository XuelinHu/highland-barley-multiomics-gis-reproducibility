"""Recover and verify reviewer-requested evidence without changing historical outputs.

conda run -n pyg python paper/revision2/scripts/complete_revision_analysis.py \
  --evidence .tmp/revision3/evidence --output paper/revision2/analysis/completed
"""
import argparse
import gzip
import hashlib
import json
import platform
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fasta(path):
    opener = gzip.open if path.suffix == ".gz" else open
    records, name, seq = {}, None, []
    with opener(path, "rt") as stream:
        for line in stream:
            if line.startswith(">"):
                if name is not None:
                    records[name] = "".join(seq).rstrip("*")
                name, seq = line[1:].split()[0], []
            elif line.strip():
                seq.append(line.strip())
    if name is not None:
        records[name] = "".join(seq).rstrip("*")
    return records


def write_table(frame, path):
    frame.to_csv(path, sep="\t", index=False, float_format="%.12g")


def pdb_sequence(path):
    aa = dict(zip(
        "ALA ARG ASN ASP CYS GLN GLU GLY HIS ILE LEU LYS MET PHE PRO SER THR TRP TYR VAL".split(),
        "ARNDCQEGHILKMFPSTWYV"))
    result, scores, seen = [], [], set()
    for line in path.read_text().splitlines():
        key = line[21:27]
        if line.startswith("ATOM") and line[12:16].strip() == "CA" and key not in seen:
            seen.add(key)
            result.append(aa.get(line[17:20].strip(), "X"))
            scores.append(float(line[60:66]))
    return "".join(result), scores


def model_cohorts(root, out, frame):
    model = root / "gis/processed/modeling"
    points = frame[frame.occurrenceStatus.fillna("PRESENT").eq("PRESENT")].copy()
    soil = [c for c in points if c.startswith("sg_")]
    wc = [f"wc_bio{i}" for i in range(1, 20)]
    subsets = {
        "T1_qingke_global": points[points.elev.notna()],
        "T2_highalt_global": points[points.elev.notna()],
        "T3_qingke_plateau": points[points.on_plateau.astype(bool)],
    }
    cv = pd.read_csv(model / "cv_metrics_present.tsv", sep="\t")
    rows, assignments = [], []
    for target, sub in subsets.items():
        sub = sub.reset_index(drop=True)
        y = (sub.elev.ge(2000) if target.startswith("T2") else sub.is_qingke_candidate).astype(int).to_numpy()
        blocks = np.load(model / f"blocks_{target}_present.npy")
        assert len(sub) == len(blocks), target
        for i in range(len(sub)):
            assignments.append(dict(target=target, row_in_target=i, gbif_key=str(sub.iloc[i].gbif_key),
                                    block=int(blocks[i]), label=int(y[i])))
        for _, original in cv[cv.target.eq(target)].iterrows():
            test = blocks == original.fold
            train = ~test
            expected = [int(train.sum()), int(test.sum()), int(y[train].sum()), int(y[test].sum())]
            actual = [int(original[c]) for c in ["n_train", "n_test", "n_pos_train", "n_pos_test"]]
            assert expected == actual, (target, original.fold, actual, expected)
            record = original.to_dict()
            record.update(cohort_n=len(sub), cohort_positive_n=int(y.sum()),
                          n_neg_train=int(train.sum() - y[train].sum()),
                          n_neg_test=int(test.sum() - y[test].sum()),
                          test_prevalence=float(y[test].mean()),
                          train_soil16_observed=int(sub.loc[train, soil].notna().all(axis=1).sum()),
                          test_soil16_observed=int(sub.loc[test, soil].notna().all(axis=1).sum()),
                          train_wc19_observed=int(sub.loc[train, wc].notna().all(axis=1).sum()),
                          test_wc19_observed=int(sub.loc[test, wc].notna().all(axis=1).sum()))
            rows.append(record)
    result = pd.DataFrame(rows)
    write_table(result, out / "model_fold_effective_n.tsv")
    write_table(pd.DataFrame(assignments), out / "model_fold_assignments.tsv")
    valid = result[result.valid.astype(bool)]
    summary = valid.groupby(["target", "feature_set", "model"], sort=False).agg(
        cohort_n=("cohort_n", "first"), positives=("cohort_positive_n", "first"),
        scorable_folds=("fold", "count"), n_train_min=("n_train", "min"),
        n_train_max=("n_train", "max"), n_test_min=("n_test", "min"),
        n_test_max=("n_test", "max"), test_prevalence_min=("test_prevalence", "min"),
        test_prevalence_max=("test_prevalence", "max"), ap_mean=("pr_auc", "mean"),
        ap_sd=("pr_auc", "std"), roc_auc_mean=("auc", "mean"), roc_auc_sd=("auc", "std")
    ).reset_index()
    old = pd.read_csv(model / "metrics_summary_present.tsv", sep="\t")
    for _, row in summary.iterrows():
        matched = old[old.target.eq(row.target) & old.feature_set.eq(row.feature_set) & old.model.eq(row.model)].iloc[0]
        assert np.isclose(row.ap_mean, matched.pr_auc_mean, atol=1e-10)
        assert np.isclose(row.roc_auc_mean, matched.auc_mean, atol=1e-10)
    write_table(summary, out / "model_summary_verified.tsv")
    return {"fold_rows": len(result), "all_fold_counts_verified": True,
            "all_aggregate_metrics_reproduced": True}


def structural_panel(root, out):
    base = root / "processed/reanalysis"
    struct = base / "structures"
    targets = json.loads((struct / "targets.json").read_text())
    full = fasta(root / "reference_pep.fa.gz")
    old = fasta(base / "top100_proteins.fasta")
    ranked = pd.read_csv(base / "drought_consensus_genes.tsv", sep="\t")
    ranked["consensus_rank"] = np.arange(1, len(ranked) + 1)
    ranked = ranked.set_index("gene_id")
    summary = pd.read_csv(struct / "structure_plddt_summary.tsv", sep="\t").set_index("gene_id")
    reasons = {
        "HORVU.MOREX.r3.4HG0386900": "Annotation-discrepancy case: CHP01589-domain evidence versus a historical NAC label",
        "HORVU.MOREX.r3.5HG0524540": "NAC/NAM transcription-factor comparator",
        "HORVU.MOREX.r3.3HG0318380": "LRR-related receptor-like protein and potentially incomplete gene-model case",
        "HORVU.MOREX.r3.1HG0046190": "Upregulated aquaporin TIP3-1 family example",
        "HORVU.MOREX.r3.6HG0596770": "Downregulated aquaporin PIP2-2 family example",
        "HORVU.MOREX.r3.2HG0099870": "Heat-shock HSP90 family example",
        "HORVU.MOREX.r3.5HG0507410": "Dehydrin DHN2 family example with an identical-sequence AFDB entry",
        "HORVU.MOREX.r3.1HG0080020": "LEA/HVA1 family example with an identical-sequence AFDB entry",
        "HORVU.MOREX.r3.7HG0637810": "Highest-ranked LEA_1 candidate; comparison of low-confidence protein models",
        "HORVU.MOREX.r3.2HG0197000": "Second LEA_1 candidate for comparison within the same family",
    }
    rows, published_fasta, hashes = [], [], []
    for gene, target in targets.items():
        seq, protein = target["seq"], target["protein_id"]
        assert full[protein] == seq, protein
        rank = int(ranked.loc[gene, "consensus_rank"])
        assert rank <= 100, (gene, rank)
        model = summary.loc[gene]
        if str(model.source).startswith("AFDB:"):
            ac = model.source.split(":")[1].split()[0]
            pdb = struct / "afdb" / f"AF-{ac}-F1.pdb"
        else:
            pdb = struct / "pdb" / f"{gene}_colabfold_rank001.pdb"
        modeled, scores = pdb_sequence(pdb)
        assert modeled == seq, (gene, "PDB/input sequence mismatch")
        assert abs(np.mean(scores) - model["mean"]) < 0.02
        rows.append(dict(gene_id=gene, protein_id=protein, consensus_rank=rank,
                         direction=ranked.loc[gene, "direction"], mean_abs_logFC=ranked.loc[gene, "mean_abs_logFC"],
                         label=model.label.replace("Dehydrin DHN1", "Dehydrin DHN2"),
                         length=len(seq), sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                         source=model.source, model_file=str(pdb.relative_to(struct)),
                         model_sha256=sha256(pdb), mean_plddt=float(np.mean(scores)),
                         selection_type="purposive panel from top100; not top10 by rank",
                         inclusion_rationale=reasons[gene]))
    write_table(pd.DataFrame(rows), out / "protein_panel_manifest.tsv")
    for protein, previous_seq in old.items():
        assert protein in full, protein
        seq = full[protein]
        assert seq.startswith(previous_seq), protein
        published_fasta.append(f">{protein}\n" + "\n".join(seq[i:i+60] for i in range(0, len(seq), 60)))
        hashes.append(dict(protein_id=protein, historical_annotation_length=len(previous_seq),
                           full_length=len(seq), sequence_sha256=hashlib.sha256(seq.encode()).hexdigest(),
                           historical_input_is_prefix=True))
    (out / "top100_proteins_full_length.fasta").write_text("\n".join(published_fasta) + "\n")
    write_table(pd.DataFrame(hashes), out / "protein_sequence_audit.tsv")
    return {"panel_size": len(rows), "reference_and_pdb_sequences_verified": True,
            "panel_ranks": [r["consensus_rank"] for r in rows],
            "historical_annotation_truncation": sorted(set(len(s) for s in old.values())),
            "full_sequences_recovered": len(hashes)}


def gradient_intervals(out, frame):
    present = frame.occurrenceStatus.fillna("PRESENT").eq("PRESENT")
    rows = []
    for name, mask in [("plateau_2000_5000m", frame.on_plateau.astype(bool)),
                       ("global_2000_5000m", pd.Series(True, index=frame.index))]:
        sub = frame[present & mask & frame.elev.between(2000, 5000)]
        for variable in ["wc_bio1", "wc_bio5", "wc_bio12"]:
            x, y = sub.elev.to_numpy() / 1000, sub[variable].to_numpy()
            good = np.isfinite(x) & np.isfinite(y)
            if variable == "wc_bio12":
                good &= y > 0
            x, y = x[good], y[good]
            log = variable == "wc_bio12"
            if log:
                y = np.log10(y)
            fit = stats.linregress(x, y)
            delta = stats.t.ppf(.975, len(x) - 2) * fit.stderr
            native = np.array([fit.slope, fit.slope - delta, fit.slope + delta])
            values = 100 * np.expm1(np.log(10) * native) if log else native
            rows.append(dict(variable=variable, subset=name, n=len(x), slope_per_km=values[0],
                             ci95_lower=values[1], ci95_upper=values[2],
                             native_slope=fit.slope, native_se=fit.stderr,
                             native_ci95_lower=native[1], native_ci95_upper=native[2],
                             r_squared=fit.rvalue**2, reported_unit="percent/km" if log else "degC/km",
                             uncertainty="nominal OLS; spatial dependence and raster error not included"))
    write_table(pd.DataFrame(rows), out / "gradient_confidence_intervals.tsv")
    return [r for r in rows if r["subset"].startswith("plateau")]


def exploratory_rnaseq(root, out):
    reference_ids, legacy_aliases = set(), {}
    with gzip.open(root / "reference.gff3.gz", "rt") as stream:
        for line in stream:
            if line.startswith("#"):
                continue
            fields = line.rstrip().split("\t")
            if len(fields) != 9 or fields[2] != "gene":
                continue
            attrs = dict(x.split("=", 1) for x in fields[8].split(";") if "=" in x)
            gene = attrs.get("ID", "").removeprefix("gene:")
            reference_ids.add(gene)
            for key in ["Alias", "Name"]:
                for alias in attrs.get(key, "").split(","):
                    if alias:
                        legacy_aliases.setdefault(alias, set()).add(gene)
    rows, examples = [], []
    for contrast in ["W8_drought_vs_ck", "W24_drought_vs_ck", "mean_drought_vs_ck"]:
        table = pd.read_csv(root / f"processed/reanalysis/GSE225255_{contrast}_limma.tsv", sep="\t")
        significant = table[table["adj.P.Val"].lt(.05) & table.logFC.abs().ge(1)]
        ids = table.gene_id.astype(str)
        direct = ids.isin(reference_ids)
        aliases = ids.map(lambda x: len(legacy_aliases.get(x, set())))
        rows.append(dict(dataset="GSE225255", contrast=contrast, tested_features=len(table),
                         fdr05_abslogfc1=len(significant), up=int(significant.logFC.gt(0).sum()),
                         down=int(significant.logFC.lt(0).sum()), exact_v3_gene_ids=int(direct.sum()),
                         unique_v3_name_or_alias_matches=int(aliases.eq(1).sum()),
                         ambiguous_v3_alias_matches=int(aliases.gt(1).sum()),
                         identifier_audit_date="2026-11-03", model="limma on deposited log2(FPKM+1); exploratory"))
        for identifier in ids.head(8):
            examples.append(dict(contrast=contrast, native_id=identifier,
                                 direct_morexv3_match=identifier in reference_ids,
                                 gff_alias_matches=";".join(sorted(legacy_aliases.get(identifier, set()))),
                                 decision="excluded without a validated cross-version bridge or legacy sequence"))
    write_table(pd.DataFrame(rows), out / "GSE225255_exploratory_summary.tsv")
    write_table(pd.DataFrame(examples), out / "GSE225255_identifier_examples.tsv")
    return rows


def sources(root, out, frame):
    ids = frame[["gbif_key", "gbifID", "datasetKey", "license", "decimalLongitude", "decimalLatitude",
                 "year", "basisOfRecord", "occurrenceStatus", "on_plateau", "is_qingke_candidate"]].copy()
    write_table(ids, out / "occurrence_record_manifest.tsv")
    counts = frame.groupby(["datasetKey", "license"], dropna=False).size().rename("records").reset_index()
    counts["source_url"] = "https://www.gbif.org/dataset/" + counts.datasetKey
    write_table(counts, out / "GBIF_contributing_datasets.tsv")
    meta = dict(source="GBIF occurrence search API; not a DOI-assigned GBIF download",
                query=dict(taxonKey=2706056, hasCoordinate=True, hasGeospatialIssue=False),
                retrieval_date="2026-09-11", original_global_records=93536,
                original_regional_records=9302, retained_matrix_records=len(frame),
                snapshot_sha256=sha256(root / "gis/processed/barley_env_matrix.tsv.gz"),
                doi_status="No historical download DOI exists; a GBIF derived-dataset DOI requires registration",
                derived_dataset_guidance="https://www.gbif.org/derived-dataset/about")
    (out / "source_provenance.json").write_text(json.dumps(meta, indent=2))


def cold_direction_tests(root, out):
    """Preserve both direction marginals within each shared gene set."""
    base = root / "processed/reanalysis"
    drought = pd.read_csv(base / "drought_consensus_genes.tsv", sep="\t")[["gene_id", "direction"]]
    rows, pairs = [], []
    rng = np.random.default_rng(42)
    for dataset, platform, contrasts in [
        ("GSE41513", "GPL16130", ["shoot_3H", "shoot_6H", "root_3H", "root_6H"]),
        ("GSE47882", "GPL17286", ["d1_vs_d0", "d7_vs_d0", "d21_vs_d0"]),
    ]:
        mapping = pd.read_csv(root / f"processed/mapping/{platform}_morex_mapping.tsv", sep="\t", dtype=str)
        mapping = mapping[mapping.status.eq("unique_gene")].set_index("probe_id").gene_ids
        for contrast in contrasts:
            table = pd.read_csv(base / f"cold/{dataset}_{contrast}_limma.tsv", sep="\t", index_col=0)
            table["gene_id"] = table.index.astype(str).map(mapping)
            sig = table[table["adj.P.Val"].lt(.05) & table.logFC.abs().ge(1)].dropna(subset=["gene_id"])
            sig = sig.iloc[np.argsort(-sig.logFC.abs().to_numpy(), kind="stable")].drop_duplicates("gene_id")
            shared = sig.merge(drought, on="gene_id", how="inner")
            n = len(shared)
            a = int(shared.direction.eq("up").sum())
            b = int(shared.logFC.gt(0).sum())
            up_up = int((shared.direction.eq("up") & shared.logFC.gt(0)).sum())
            same = n - a - b + 2 * up_up
            permuted = n - a - b + 2 * rng.hypergeometric(a, n - a, b, size=10000)
            p = (1 + int((permuted >= same).sum())) / 10001
            rows.append(dict(dataset=dataset, contrast=contrast, compared_unique_genes=n,
                             concordant=same, drought_up=a, cold_up=b, up_up=up_up,
                             concordance_fraction=same/n, permutation_expected_concordance=float(permuted.mean()),
                             permutation_p_one_sided=p, exact_fixed_margin_p=stats.hypergeom.sf(up_up-1, n, a, b),
                             permutations=10000, seed=42))
            for _, item in shared.iterrows():
                pairs.append(dict(dataset=dataset, contrast=contrast, gene_id=item.gene_id,
                                  drought_direction=item.direction, cold_logFC=item.logFC,
                                  cold_FDR=item["adj.P.Val"], concordant=(item.direction == "up") == (item.logFC > 0)))
    result = pd.DataFrame(rows)
    order = np.argsort(result.permutation_p_one_sided.to_numpy())
    ranked_p = result.permutation_p_one_sided.to_numpy()[order]
    q = np.minimum.accumulate((ranked_p * len(ranked_p) / np.arange(1, len(ranked_p)+1))[::-1])[::-1]
    result["permutation_BH_q"] = np.nan
    result.loc[order, "permutation_BH_q"] = np.minimum(q, 1)
    write_table(result, out / "cold_direction_fixed_margin_tests.tsv")
    pair_frame = pd.DataFrame(pairs)
    write_table(pair_frame, out / "cold_direction_gene_contrast_pairs.tsv")
    return dict(contrasts=len(rows), gse41513_gene_contrast_pairs=int(pair_frame.dataset.eq("GSE41513").sum()),
                gse41513_distinct_genes=int(pair_frame[pair_frame.dataset.eq("GSE41513")].gene_id.nunique()),
                note="Gene-label exchangeability remains an approximation; co-regulation is not removed")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    frame = pd.read_csv(args.evidence / "gis/processed/barley_env_matrix.tsv.gz", sep="\t", low_memory=False)
    report = dict(
        model_cohorts=model_cohorts(args.evidence, args.output, frame),
        proteins=structural_panel(args.evidence, args.output),
        gradients=gradient_intervals(args.output, frame),
        gse225255=exploratory_rnaseq(args.evidence, args.output),
        cold=cold_direction_tests(args.evidence, args.output),
        runtime=dict(python=platform.python_version(), numpy=np.__version__, pandas=pd.__version__),
    )
    sources(args.evidence, args.output, frame)
    (args.output / "evidence_verification.json").write_text(json.dumps(report, indent=2))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
