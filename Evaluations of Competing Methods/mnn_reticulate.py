from anndata import AnnData
import os
import scanpy as sc
import numpy as np
import pandas as pd
from scipy.io import mmread
from CarDEC.CarDEC_utils import normalize_scanpy

def get_hvgs(array, batch_vec, gene_vec):
    adata = AnnData(array)
    
    n, p = adata.shape
    
    if n != len(batch_vec):
        adata = adata.T
        
    adata.obs['batch'] = batch_vec
    adata.var.index = gene_vec
    
    adata = normalize_scanpy(adata, batch_key = "batch", n_high_var = 2000, LVG = True)
    
    return adata.var.index[adata.var["Variance Type"] == "HVG"].values

def read_pancreas(path):
    pathlist = os.listdir(path)
    adata = sc.read(os.path.join(path, pathlist[0]))
          
    for i in range(1,len(pathlist)):
        adata = adata.concatenate(sc.read(os.path.join(path, pathlist[i])))

    sc.pp.filter_cells(adata, min_genes = 200)
    sc.pp.filter_genes(adata, min_cells = 30)
    mito_genes = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    
    notmito_genes = [not x for x in mito_genes]
    adata = adata[:,notmito_genes]
    del adata.obs['batch']
    
    return adata.X.toarray(), adata.obs, list(adata.var.index)

def read_retina(path):
    adata = sc.read(os.path.join(path, 'matrix.mtx')).T
    genes_file = pd.read_csv(os.path.join(path, 'genes.tsv'), sep='\t')
    barcodes_file = pd.read_csv(os.path.join(path, 'barcodes.tsv'), sep='\t')
    
    adata.var_names = genes_file["genename"]
    adata.obs_names = barcodes_file["cellname"]
    adata.var['genenames'] = genes_file["genename"].values
    adata.obs['celltype'] = barcodes_file["celltype"].values
    adata.obs['celltype_com'] = barcodes_file["celltype_com"].values
    adata.obs['celltypeID'] = barcodes_file["celltypeID"].values
    adata.obs['BatchID'] = barcodes_file["BatchID"].values
    
    sc.pp.filter_cells(adata, min_genes=200)
    sc.pp.filter_genes(adata, min_cells=30)
    
    mito_genes = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata = adata[adata.obs['n_genes'] < 2500, :]
    adata = adata[adata.obs['percent_mito'] < 0.05, :]
    
    return adata.X.toarray(), adata.obs, list(adata.var.index)

def read_pbmc(path):
    adata = AnnData(mmread(os.path.join(path, "counts.umi.txt"))).T
    adata.X = adata.X.toarray()
    adata.obs.index = pd.read_csv(os.path.join(path, 'cells.umi.new.txt'), sep = '\t')['cellname']
    adata.var.index = pd.read_csv(os.path.join(path, 'genes.txt'), sep = '\t')['scientific']
    
    meta_data = pd.read_csv(os.path.join(path, 'meta.txt'), sep = '\t')
    meta_data = meta_data.iloc[1:,]
    meta_data.index = meta_data['NAME'].values
    meta_data = meta_data.iloc[[x in adata.obs.index for x in meta_data.index],]
    adata = adata[[x in meta_data.index for x in adata.obs.index],]
    adata.obs = meta_data
    
    sc.pp.filter_cells(adata, min_genes = 200)
    sc.pp.filter_genes(adata, min_cells = 30)
    
    adata = adata[adata.obs["CellType"].values != "Unassigned"]
    
    return adata.X.toarray(), adata.obs, list(adata.var.index)

def read_macaque(path):
    adata = sc.read(path)
    sc.pp.filter_cells(adata, min_genes = 200)
    sc.pp.filter_genes(adata, min_cells = 30)

    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata = adata[adata.obs['n_genes'] < 2500, :]
    
    return adata.X.toarray(), adata.obs, list(adata.var.index)

def read_cortex(path):
    adata = AnnData(mmread(os.path.join(path, "count.umis.txt"))).T
    urls_ = ["meta_combined.txt", "genes.txt", "cell.names.new.txt"]
    urls_ = [os.path.join(path, x) for x in urls_]
    metadata, genes, cells = [pd.read_csv(url_, sep = '\t') for url_ in urls_]
    genes, cells = list(genes['scientific']), list(cells['cellname'])
    indices = [x in list(metadata["NAME"]) for x in cells]
    adata.obs.index, adata.var.index = cells, genes
    adata.X = adata.X.toarray()
    adata = adata[indices]
    adata.obs = metadata
    
    sc.pp.filter_cells(adata, min_genes = 200)
    sc.pp.filter_genes(adata, min_cells = 30)
    
    adata = adata[adata.obs["CellType"].values != "Unassigned"]
    
    return adata.X, adata.obs, list(adata.var.index)

def read_liver(path):
    adata = sc.read_mtx(os.path.join(path, 'matrix.mtx')).T
    genes_file = pd.read_csv(os.path.join(path, 'genes.tsv'), sep='\t')
    barcodes_file = pd.read_csv(os.path.join(path, 'barcodes.tsv'), sep='\t')

    adata.var.index = genes_file["genename"]
    adata.obs.index = barcodes_file["cellname"]
    adata.obs = barcodes_file
        
    sc.pp.filter_cells(adata, min_genes = 200)
    mito_genes = adata.var_names.str.startswith('mt-')
    adata.obs['percent_mito'] = np.sum(
        adata[:, mito_genes].X, axis=1).A1 / np.sum(adata.X, axis=1).A1
    adata.obs['n_counts'] = adata.X.sum(axis=1).A1
    adata = adata[adata.obs['percent_mito'] < 0.2, :]
    sc.pp.filter_genes(adata, min_cells = 30)

    return adata.X, adata.obs, list(adata.var.index)