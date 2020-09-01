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
# extract meta to pandas DF and output to csv
```
import pandas as pd
obsm_data=pd.DataFrame(adata.obsm['X_umap'])
obsm_data.to_csv("umap.csv", sep=",")
```

# In R, build a anndata object from a seurat obj
```
# you need to set up your own conda env with scanpy installed within it
# see scanpy guides for help
Sys.setenv(RETICULATE_PYTHON = "/data/mcgaugheyd/conda/envs/scanpy/bin/python")
library(reticulate)
sc <- import("scanpy")
exprs <- GetAssayData(seurat)
meta <- seurat[[]]
feature_meta <- GetAssay(seurat)[[]]
reduction <- 'PCA'
embedding <- Embeddings(seurat, reduction)
# if you are clustering then you can
# build anndata obj with a subset of the genes counts
# as they aren't used in clustering
# and it makes the obj creation far faster with less mem usage
# adata_seurat = sc$AnnData(X = t(as.matrix(exprs[1:10,])), obs = meta, var = feature_meta[1:10,])
adata_seurat = sc$AnnData(X = t(as.matrix(exprs)), obs = meta, var = feature_meta)
```
