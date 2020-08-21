# scanpy 1.6 notes

# import h5ad anndata
import scanpy as sc
adata = sc.read_h5ad(
    'your_anndata.h5ad')
    
# make bigger plot
sc.set_figure_params(figsize=[20,20], dpi = 200)

# scatter plot with custom coords
sc.pl.scatter(adata, basis = 'scviumap', color = ['batch'])
