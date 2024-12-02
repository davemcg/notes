# scanpy 1.6 notes

# set parallelization globally
`sc.settings.n_jobs = 8`

# import h5ad anndata
```
import scanpy as sc
adata = sc.read_h5ad(
    'your_anndata.h5ad')
``` 
# make bigger plot
```
sc.set_figure_params(figsize=[20,20], dpi = 200)
```
# scatter plot with custom coords
```
sc.pl.scatter(adata, size = 5, basis = 'scviumap', color = ['batch'])
```

# scatter plot (gene) to file
```
# sc.pp.log1p(adata) 
sc.pl.scatter(adata, size = 5, basis = 'umap', color = ['ENSG00000163914'], save = 'rho.png',use_raw=False, show = False)
```

# extract meta to pandas DF and output to csv
```
import pandas as pd
obsm_data=pd.DataFrame(adata.obsm['X_umap'])
obsm_data.to_csv("umap.csv", sep=",")
```

# Remove cells with certain critera in `obs`
```
adata_new = adata[~adata.obs['CellType'].isin(['NA', 'Doublet','Doublets']),:]
```

# In R, build a anndata object from a seurat obj
```
# you should set up a conda env with scanpy 
# see scanpy guides for help
Sys.setenv(RETICULATE_PYTHON = "/data/mcgaugheyd/conda/envs/scanpy/bin/python")
library(reticulate)
sc <- import("scanpy")
library(Seurat)
library(Matrix)
exprs <- GetAssayData(seurat)
meta <- seurat[[]]
feature_meta <- GetAssay(seurat)[[]]
reduction <- 'PCA'
embedding <- Embeddings(seurat, reduction)
# if you are clustering then you can
# build anndata obj with a subset of the genes counts
# as they aren't used in clustering
# and it makes the obj creation far faster with less mem usage
# adata_seurat = sc$AnnData(X = t(exprs[1:10,]), obs = meta, var = feature_meta[1:10,])
adata_seurat = sc$AnnData(X = t(exprs), obs = meta, var = feature_meta)
adata_seurat$obsm$update(X_pca = embedding)
```

# Build Count Matrix
```
# first subset to certain cells (rods)
rod_adata = adata[adata.obs['MajorCellType'] == 'rod']
rod_df = pd.DataFrame(data=rod_adata.X.toarray(), index=rod_adata.obs_names, columns=rod_adata.var_names)
rod_df.to_csv('~/Desktop/rod_df.csv.gz')
```

# Counts across multiple columns (pandas)
```
adata.obs[['Manual_CellType', 'Machine_CellType']].apply(pd.Series.value_counts)
```
