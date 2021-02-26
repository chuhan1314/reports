import stlearn as st
import csv
import numpy as np
import pandas as pd
import scanpy as sc

st.settings.set_figure_params(dpi=600)
# Reading data
data = st.Read10X(path="E:/data/work/kongjianshujv/e13_beirui/out")
#st.add.image(adata=data, imgpath="E:/data/work/kongjianshujv/e13_beirui/out/spatial/tissue_hires_image.png", library_id="out",visium=True)
# Save raw_count
data.layers["raw_count"] = data.X
# Preprocessing
st.pp.filter_genes(data,min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.layers["normal_count"] = data.X
# Keep raw data
data.raw = data
st.pp.scale(data)

st.add.labels(data, 'label_transfer_bc.csv')
st.pl.cluster_plot(data,use_label="predictions", name='label_transfer', output='.',spot_size=20)

#区域
st.tl.cci.het.count(data, use_clustering='label_transfer')
st.pl.het_plot(data, use_het='cci_het', name='het_louvain', output='.',spot_size=20)
data.uns["lr"] = ['Tgfb1_Tgfbr1_Tgfbr2']
st.tl.cci.base.lr(adata=data)
st.pl.het_plot(data, use_het='cci_lr', data_alpha=0.7, name='cci_lr', output='.',spot_size=20)
st.tl.cci.merge(data, use_lr='cci_lr', use_het='cci_het')
st.pl.het_plot(data, use_het='merged', data_alpha=0.7, name='merged', output='.',spot_size=20)
st.tl.cci.permutation(data, n_pairs=200)
st.pl.het_plot(data, use_het='merged_pvalues', data_alpha=0.7, name='permutation', output='.',spot_size=20)
st.pl.het_plot(data, use_het='merged_sign', data_alpha=0.7, name='final', output='.',spot_size=20)

#spot内
st.tl.cci.het.count(data, use_clustering='label_transfer', distance=0)# distance=0表示点内
st.pl.het_plot(data, use_het='cci_het', name='het_louvain', output='.',spot_size=20)
data.uns["lr"] = ['Tgfb1_Tgfbr1_Tgfbr2']
st.tl.cci.base.lr(adata=data, distance=0)
st.pl.het_plot(data, use_het='cci_lr', data_alpha=0.7, name='cci_lr', output='.',spot_size=20)
st.tl.cci.merge(data, use_lr='cci_lr', use_het='cci_het')
st.pl.het_plot(data, use_het='merged', data_alpha=0.7, name='merged', output='.',spot_size=20)
st.tl.cci.permutation(data, n_pairs=200, distance=0)
st.pl.het_plot(data, use_het='merged_pvalues', data_alpha=0.7, name='permutation', output='.',spot_size=20)
st.pl.het_plot(data, use_het='merged_sign', data_alpha=0.7, name='final', output='.',spot_size=20)


