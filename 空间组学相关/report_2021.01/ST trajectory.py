import stlearn as st
import csv
import numpy as np
import pandas as pd
import scanpy as sc

st.settings.set_figure_params(dpi=600)
# Reading data
data = st.Read10X(path="E:/data/work/kongjianshujv/e13_beirui/out")
# Save raw_count
sc.pl.highest_expr_genes(data, n_top=20)
data.layers["raw_count"] = data.X
# Preprocessing
st.pp.filter_genes(data,min_cells=3)
st.pp.normalize_total(data)
st.pp.log1p(data)
data.layers["normal_count"] = data.X
# Keep raw data
data.raw = data
st.pp.scale(data)
data.layers["scale_count"] = data.X

# Run PCA
st.em.run_pca(data,n_comps=50,random_state=0)
# Tiling image
st.pp.tiling(data,out_path="tiling",crop_size = 40)
# Using Deep Learning to extract feature
st.pp.extract_feature(data)
# Apply stSME spatial-PCA option
st.spatial.morphology.adjust(data,use_data="X_pca",radius=50,method="mean")
st.pp.neighbors(data,n_neighbors=50,use_rep='X_pca_morphology',random_state=0)
st.tl.clustering.louvain(data,random_state=0)
#sc.tl.umap(data)
#sc.pl.umap(data, color='louvain')
st.pl.cluster_plot(data,use_label="louvain",tissue_alpha=1,spot_size=5,show_legend=True)
st.pl.cluster_plot(data,use_label="louvain",tissue_alpha=1,spot_size=5,list_cluster=[2],show_legend=True)

data.uns["iroot"] = 1427
st.spatial.trajectory.pseudotime(data,eps=50,use_rep="X_pca",use_sme=True)
if 1!=1:
    x=pd.read_csv("./meta.csv")
    dct_data = pd.DataFrame(x.loc[:, :],)
    z=pd.DataFrame(dct_data.loc[5])
    z=z.T
    z=z.drop(['Unnamed: 0'],axis=1)
    z.index=['seurat']
    y=data.obs
    y=y.T
    y = y.append(z)
    y=y.T
    data.obs=y
st.spatial.trajectory.pseudotimespace_global(data,use_label="louvain",list_cluster=[0,2])
st.pl.cluster_plot(data,use_label="louvain",show_trajectory=True,list_cluster=[0,2],show_subcluster=False,spot_size=15)
st.pl.trajectory.tree_plot(data)
st.spatial.trajectory.detect_transition_markers_clades(data,clade=2,use_raw_count=False,cutoff_spearman=0.3)
st.pl.trajectory.transition_markers_plot(data,top_genes=30,trajectory="clade_2")
st.pl.gene_plot(data,genes="Cxcl12",list_clusters=[0,2])

data.uns["iroot"] = np.flatnonzero(data.obs["louvain"]  == str(3))[50]
st.spatial.trajectory.pseudotime(data,eps=50,use_rep="X_pca")
st.pl.non_spatial_plot(data,use_label="louvain")
st.pl.trajectory.pseudotime_plot(data,list_cluster="all",show_graph=True,node_alpha=1,tissue_alpha=1,edge_alpha=0.1,node_size=5)
st.spatial.trajectory.pseudotimespace_local(data,use_label="louvain",cluster=2)
st.pl.subcluster_plot(data,use_label="louvain",cluster=2,tissue_alpha=1)
st.pl.trajectory.local_plot(data,use_cluster=2,branch_alpha=0.2,reverse=True)

st.spatial.trajectory.pseudotimespace_global(data,use_label="louvain",list_cluster=[0,2])
st.pl.cluster_plot(data,use_label="louvain",show_trajectory=True,list_cluster=[0,2],show_subcluster=False)
