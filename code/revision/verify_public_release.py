"""Independent consistency checks for the public numerical archive."""
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd
from sklearn.metrics import average_precision_score, roc_auc_score
from sklearn.metrics.pairwise import haversine_distances


def main():
    p=argparse.ArgumentParser();p.add_argument('--root',type=Path,default=Path('.'));a=p.parse_args()
    checks=json.loads((a.root/'SHA256SUMS.json').read_text())
    for name,digest in checks.items():
        assert hashlib.sha256((a.root/name).read_bytes()).hexdigest()==digest,name
    result=a.root/'results'
    predictions=pd.read_csv(result/'spatial_sensitivity_predictions.tsv',sep='\t')
    folds=pd.read_csv(result/'spatial_sensitivity_folds.tsv',sep='\t')
    keys=['analysis','climate_product','thinning_km','model','fold']
    verified=0
    for key,group in predictions.groupby(keys):
        mask=np.ones(len(folds),dtype=bool)
        for k,v in zip(keys,key):mask &= folds[k].eq(v).to_numpy()
        row=folds[mask].iloc[0]
        assert len(group)==row.n_test and group.label.sum()==row.n_pos_test
        assert group.gbif_key.nunique()==len(group)
        assert np.isclose(average_precision_score(group.label,group.probability),row.ap,atol=1e-10)
        assert np.isclose(roc_auc_score(group.label,group.probability),row.roc_auc,atol=1e-10)
        verified+=1
    assignments=pd.read_csv(result/'model_fold_assignments.tsv',sep='\t')
    historic=pd.read_csv(result/'model_fold_effective_n.tsv',sep='\t')
    for row in historic.itertuples():
        sub=assignments[assignments.target.eq(row.target)]
        train=sub[sub.block.ne(row.fold)];test=sub[sub.block.eq(row.fold)]
        assert set(train.gbif_key).isdisjoint(set(test.gbif_key))
        assert (len(train),len(test),train.label.sum(),test.label.sum())==(row.n_train,row.n_test,row.n_pos_train,row.n_pos_test)
    cold=pd.read_csv(result/'cold_direction_fixed_margin_tests.tsv',sep='\t')
    first=cold[cold.dataset.eq('GSE41513')]
    assert first.compared_unique_genes.sum()==33 and first.concordant.sum()==28
    assert first.permutation_p_one_sided.eq(1).all()
    panel=pd.read_csv(result/'protein_panel_manifest.tsv',sep='\t')
    assert len(panel)==10 and panel.gene_id.nunique()==10
    for row in panel.itertuples():
        model=a.root/'evidence/processed/reanalysis/structures'/row.model_file.replace('\\','/')
        assert hashlib.sha256(model.read_bytes()).hexdigest()==row.model_sha256
    audit=pd.read_csv(result/'protein_sequence_audit.tsv',sep='\t')
    assert len(audit)==107
    retained=pd.read_csv(result/'spatial_sensitivity_record_manifest.tsv',sep='\t')
    points=pd.read_csv(result/'occurrence_record_manifest.tsv',sep='\t').set_index('gbif_key')
    for distance,group in retained[retained.analysis.eq('spatial_thinning')].groupby('thinning_km'):
        if distance<=0:continue
        ids=group.gbif_key
        coordinates=np.deg2rad(points.loc[ids,['decimalLatitude','decimalLongitude']].to_numpy())
        distances=haversine_distances(coordinates)*6371.0088
        np.fill_diagonal(distances,np.inf)
        assert distances.min()>=distance-1e-6,(distance,distances.min())
    print(json.dumps({'checksums_verified':len(checks),'sensitivity_metric_folds_recomputed':verified,
        'historical_fold_counts_verified':len(historic),'cold_count_correction_verified':True,'model_hashes_verified':10}))


if __name__=='__main__':main()
