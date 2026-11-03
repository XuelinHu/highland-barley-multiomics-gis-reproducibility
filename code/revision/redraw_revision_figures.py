"""Regenerate manuscript figures from recovered numerical evidence (PDF + 600-dpi PNG)."""
import argparse
from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm, TwoSlopeNorm
from matplotlib.patches import Rectangle
from scipy import stats


def save(fig, output, name):
    fig.savefig(output / (name + '.pdf'), bbox_inches='tight')
    fig.savefig(output / (name + '.png'), dpi=600, bbox_inches='tight')
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument('--evidence', type=Path, required=True)
    p.add_argument('--results', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    p.add_argument('--maps-only', action='store_true')
    a = p.parse_args()
    a.output.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update({'font.family':'serif', 'font.size':11, 'axes.labelsize':11,
        'axes.titlesize':12, 'xtick.labelsize':9, 'ytick.labelsize':9,
        'legend.fontsize':9, 'figure.facecolor':'white', 'axes.facecolor':'white',
        'pdf.fonttype':42, 'axes.spines.top':False, 'axes.spines.right':False})
    frame = pd.read_csv(a.evidence / 'gis/processed/barley_env_matrix.tsv.gz', sep='\t', low_memory=False)
    frame = frame[frame.occurrenceStatus.fillna('PRESENT').eq('PRESENT')]
    proxy = frame[frame.is_qingke_candidate.astype(bool)]
    plateau = frame[frame.on_plateau.astype(bool)]
    blue, orange, gray = '#0072B2', '#D55E00', '#B0B0B0'
    if a.maps_only:
        import rasterio
        from rasterio.windows import from_bounds
        from rasterio.enums import Resampling
        fig, axes = plt.subplots(3, 1, figsize=(8, 11), constrained_layout=True)
        for i, (ax, fn, title, cmap, norm) in enumerate(zip(axes,
            ['wc2.1_30s_bio_1.tif','wc2.1_30s_bio_12.tif','wc2.1_30s_elev.tif'],
            ['Mean annual temperature (°C)','Annual precipitation (mm)','Elevation (m)'],
            ['RdBu_r','YlGnBu','viridis'], [TwoSlopeNorm(vmin=-20,vcenter=0,vmax=30),LogNorm(10,8500),None])):
            with rasterio.open(a.evidence / 'gis/raw/worldclim' / fn) as ds:
                arr = ds.read(1, window=from_bounds(70,25,105,40,ds.transform),
                              out_shape=(300,700), resampling=Resampling.average, masked=True)
            im = ax.imshow(arr, extent=(70,105,25,40), origin='upper', cmap=cmap, norm=norm)
            ax.scatter(proxy.decimalLongitude, proxy.decimalLatitude, s=5, c='black', alpha=.75,
                       label='High-elevation proxy (n=558)')
            ax.set(xlim=(70,105),ylim=(25,40),xlabel='Longitude (°E)',ylabel='Latitude (°N)',title=f'({chr(65+i)}) {title}')
            fig.colorbar(im,ax=ax,shrink=.85,pad=.02)
            if i == 0: ax.legend(loc='lower left')
        save(fig,a.output,'fig4_env_maps')
        return
    base = a.evidence / 'processed/reanalysis'
    for dataset in ['GSE103278','GSE128048']:
        v = pd.read_csv(a.results / (dataset+'_volcano_values.tsv'),sep='\t')
        significant = v.omnibus_FDR.lt(.05)
        fig,ax=plt.subplots(figsize=(6,4.8),constrained_layout=True)
        ax.scatter(v.mean_logFC, -np.log10(v.omnibus_FDR.clip(lower=1e-300)),s=3,
                   c=np.where(significant,np.where(v.mean_logFC.ge(0),orange,blue),gray),alpha=.45,rasterized=True)
        ax.axhline(-np.log10(.05),color='black',lw=.7,ls='--')
        ax.set(xlabel='Mean log₂ fold change',ylabel='−log₁₀ omnibus FDR',title=dataset)
        save(fig,a.output,dataset+'_volcano')
    probes = pd.read_csv(base/'drought_consensus_probes.tsv',sep='\t')
    assert len(probes)==3076
    top=pd.read_csv(base/'top100_drought_candidates.tsv',sep='\t')
    genes=set(top.gene_id)
    selected=probes.gene_ids.fillna('').apply(lambda s:bool(set(str(s).split(';')) & genes))
    fig,ax=plt.subplots(figsize=(6.4,5),constrained_layout=True)
    ax.scatter(probes.gse103_logFC,probes.gse128_logFC,s=9,c=gray,alpha=.5,rasterized=True)
    ax.scatter(probes.loc[selected,'gse103_logFC'],probes.loc[selected,'gse128_logFC'],s=16,c=orange,label='Probes mapping to top-100 genes')
    ax.axhline(0,lw=.6,c='black');ax.axvline(0,lw=.6,c='black')
    ax.set(xlabel='GSE103278 mean log₂ fold change',ylabel='GSE128048 mean log₂ fold change')
    ax.legend(loc='upper left');save(fig,a.output,'consensus_effects')
    fig,axes=plt.subplots(1,2,figsize=(10,4),constrained_layout=True)
    for ax,df,title in [(axes[0],frame,'(A) Global occurrences'),(axes[1],plateau,'(B) Plateau window')]:
        sc=ax.scatter(df.decimalLongitude,df.decimalLatitude,c=df.elev,s=3 if ax is axes[0] else 10,cmap='viridis',vmin=0,vmax=6000,rasterized=True)
        ax.set(xlabel='Longitude (°E)',ylabel='Latitude (°N)',title=title)
    axes[0].add_patch(Rectangle((70,25),35,15,fill=False,ec=orange,lw=1.2))
    axes[1].scatter(proxy.decimalLongitude,proxy.decimalLatitude,s=14,facecolor='none',edgecolor=orange,lw=.6,label='High-elevation proxy')
    axes[1].set(xlim=(70,105),ylim=(25,40));axes[1].legend(loc='lower left')
    fig.colorbar(sc,ax=axes,label='Elevation (m)',shrink=.8);save(fig,a.output,'fig1_sample_distribution')
    edges=np.arange(-500,6501,500);centers=(edges[:-1]+edges[1:])/2
    total=np.histogram(frame.elev.dropna(),edges)[0];nproxy=np.histogram(proxy.elev,edges)[0]
    fig,axes=plt.subplots(3,1,figsize=(7,8),constrained_layout=True)
    axes[0].bar(centers,total,width=440,color=blue);axes[0].set(ylabel='Occurrence count',yscale='log',title='(A) All records with elevation')
    axes[1].bar(centers,np.divide(nproxy,total,out=np.zeros_like(total,dtype=float),where=total>0)*100,width=440,color=orange)
    axes[1].set(ylabel='Proxy share (%)',title='(B) Proxy fraction in each elevation band')
    axes[2].hist(proxy.elev,bins=edges,color=orange,edgecolor='white');axes[2].set(ylabel='Proxy count',xlabel='Elevation (m)',title='(C) High-elevation proxy')
    save(fig,a.output,'fig2_elevation_bands')
    fig,axes=plt.subplots(1,3,figsize=(12,4),constrained_layout=True)
    ci=pd.read_csv(a.results/'gradient_confidence_intervals.tsv',sep='\t')
    ci=ci[ci.subset.eq('plateau_2000_5000m')].set_index('variable')
    for i,(ax,col,label) in enumerate(zip(axes,['wc_bio1','wc_bio12','wc_bio5'],['Mean annual temperature (°C)','Annual precipitation (mm)','Warmest-month maximum (°C)'])):
        for df,color,name in [(frame,gray,'Global'),(plateau,blue,'Plateau'),(proxy,orange,'Proxy')]:
            ax.scatter(df.elev,df[col],s=3,alpha=.25,c=color,label=name,rasterized=True)
        fit=plateau[plateau.elev.between(2000,5000)].dropna(subset=[col])
        x=fit.elev.to_numpy()/1000;y=fit[col].to_numpy()
        if col=='wc_bio12':y=np.log10(y);ax.set_yscale('log');ax.set_xlim(1800,6100)
        else:ax.set_xlim(-500,6100)
        res=stats.linregress(x,y);xx=np.linspace(2,5,100)
        yy=res.intercept+res.slope*xx
        sy=np.sqrt(np.sum((y-(res.intercept+res.slope*x))**2)/(len(x)-2))
        delta=stats.t.ppf(.975,len(x)-2)*sy*np.sqrt(1/len(x)+(xx-x.mean())**2/np.sum((x-x.mean())**2))
        ax.plot(xx*1000,10**yy if col=='wc_bio12' else yy,c='black',lw=1.5)
        ax.fill_between(xx*1000,10**(yy-delta) if col=='wc_bio12' else yy-delta,10**(yy+delta) if col=='wc_bio12' else yy+delta,color='black',alpha=.15)
        r=ci.loc[col];unit='%/km' if col=='wc_bio12' else '°C/km'
        ax.set(xlabel='Elevation (m)',ylabel=label,title=f'({chr(65+i)}) n={len(x)}; R²={res.rvalue**2:.3f}')
        ax.text(.03,.03,f'{r.slope_per_km:.2f} {unit}\n95% CI [{r.ci95_lower:.2f}, {r.ci95_upper:.2f}]',transform=ax.transAxes,fontsize=9,bbox=dict(fc='white',ec='none',alpha=.85))
    axes[0].legend(loc='upper right',markerscale=2);save(fig,a.output,'fig3_env_altitude_gradients')
    for dataset in ['GSE103278','GSE128048']:
        tab=pd.read_csv(base/'wgcna'/('module_trait_'+dataset+'.tsv'),sep='\t')
        corr=tab.pivot(index='module',columns='trait',values='cor');q=tab.pivot(index='module',columns='trait',values='fdr')
        fig,ax=plt.subplots(figsize=(5.4,7),constrained_layout=True)
        im=ax.imshow(corr,vmin=-1,vmax=1,cmap='RdBu_r',aspect='auto')
        ax.set_xticks(range(len(corr.columns)),corr.columns);ax.set_yticks(range(len(corr)),corr.index)
        for iy in range(len(corr)):
            for ix in range(len(corr.columns)):
                value=corr.iloc[iy,ix];ax.text(ix,iy,f'{value:.2f}'+('*' if q.iloc[iy,ix]<.05 else ''),ha='center',va='center',fontsize=9,color='white' if abs(value)>.65 else 'black')
        ax.set(title=dataset+'; * FDR < 0.05',ylabel='Module');fig.colorbar(im,ax=ax,label='Pearson correlation',shrink=.6)
        save(fig,a.output,'module_trait_heatmap_'+dataset)
    panel=pd.read_csv(a.results/'protein_panel_manifest.tsv',sep='\t')
    fig,axes=plt.subplots(2,5,figsize=(13,5.8),constrained_layout=True)
    for ax,(_,row) in zip(axes.flat,panel.iterrows()):
        t=pd.read_csv(base/'structures'/('plddt_'+row.gene_id+'.tsv'),sep='\t')
        ax.plot(t.residue_index,t.plddt,c=blue,lw=.9);ax.axhline(70,color=gray,lw=.8,ls='--')
        ax.set(ylim=(0,100),xlabel='Residue',ylabel='pLDDT',title=row.label.split(' (')[0]+'\nrank '+str(row.consensus_rank))
        ax.tick_params(labelsize=8)
    save(fig,a.output,'structure_plddt_summary')


if __name__=='__main__':main()
