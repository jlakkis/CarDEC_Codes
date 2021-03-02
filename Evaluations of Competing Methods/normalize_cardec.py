import numpy as np
from scipy.sparse import issparse

import scanpy as sc
from anndata import AnnData


def normalize_scanpy(adata, batch_key = None, n_high_var = 1000, LVG = True, 
                     normalize_samples = True, log_normalize = True, 
                     normalize_features = True):
    """ This function preprocesses the raw count data.
    
    
    Arguments:
    ------------------------------------------------------------------
    - adata: `anndata.AnnData`, the annotated data matrix of shape (n_obs, n_vars). Rows correspond to cells and columns to genes.
    - batch_key: `str`, string specifying the name of the column in the observation dataframe which identifies the batch of each cell. If this is left as None, then all cells are assumed to be from one batch.
    - n_high_var: `int`, integer specifying the number of genes to be idntified as highly variable. E.g. if n_high_var = 2000, then the 2000 genes with the highest variance are designated as highly variable.
    - LVG: `bool`, Whether to retain and preprocess LVGs.
    - normalize_samples: `bool`, If True, normalize expression of each gene in each cell by the sum of expression counts in that cell.
    - log_normalize: `bool`, If True, log transform expression. I.e., compute log(expression + 1) for each gene, cell expression count.
    - normalize_features: `bool`, If True, z-score normalize each gene's expression.
    
    Returns:
    ------------------------------------------------------------------
    - adata: `anndata.AnnData`, the annotated data matrix of shape (n_obs, n_vars). Contains preprocessed data.
    """
    
    n, p = adata.shape
    sparsemode = issparse(adata.X)
    
    if batch_key is not None:
        batch = list(adata.obs[batch_key])
        batch = convert_vector_to_encoding(batch)
        batch = np.asarray(batch)
        batch = batch.astype('float32')
    else:
        batch = np.ones((n,), dtype = 'float32')
        norm_by_batch = False
        
    sc.pp.filter_genes(adata, min_counts=1)
    sc.pp.filter_cells(adata, min_counts=1)
        
    count = adata.X.copy()
        
    if normalize_samples:
        out = sc.pp.normalize_total(adata, inplace = False)
        obs_ = adata.obs
        var_ = adata.var
        adata = None
        adata = AnnData(out['X'])
        adata.obs = obs_
        adata.var = var_
        
        size_factors = out['norm_factor'] / np.median(out['norm_factor'])
        out = None
    else:
        size_factors = np.ones((adata.shape[0], ))
        
    if not log_normalize:
        adata_ = adata.copy()
    
    sc.pp.log1p(adata)
    
    if n_high_var is not None:
        sc.pp.highly_variable_genes(adata, inplace = True, min_mean = 0.0125, max_mean = 3, min_disp = 0.5, 
                                          n_bins = 20, n_top_genes = n_high_var, batch_key = batch_key)
        
        hvg = adata.var['highly_variable'].values
        
        if not log_normalize:
            adata = adata_.copy()

    else:
        hvg = [True] * adata.shape[1]
        
    if normalize_features:
        batch_list = np.unique(batch)

        if sparsemode:
            adata.X = adata.X.toarray()

        for batch_ in batch_list:
            indices = [x == batch_ for x in batch]
            sub_adata = adata[indices]
            
            sc.pp.scale(sub_adata)
            adata[indices] = sub_adata.X
        
        adata.layers["normalized input"] = adata.X
        adata.X = count
        adata.var['Variance Type'] = [['LVG', 'HVG'][int(x)] for x in hvg]
            
    else:
        if sparsemode:   
            adata.layers["normalized input"] = adata.X.toarray()
        else:
            adata.layers["normalized input"] = adata.X
            
        adata.var['Variance Type'] = [['LVG', 'HVG'][int(x)] for x in hvg]
        
    if n_high_var is not None:
        del_keys = ['dispersions', 'dispersions_norm', 'highly_variable', 'highly_variable_intersection', 'highly_variable_nbatches', 'means']
        del_keys = [x for x in del_keys if x in adata.var.keys()]
        adata.var = adata.var.drop(del_keys, axis = 1)
            
    y = np.unique(batch)
    num_batch = len(y)
    
    adata.obs['size factors'] = size_factors.astype('float32')
    adata.obs['batch'] = batch
    adata.uns['num_batch'] = num_batch
    
    if sparsemode:
        adata.X = adata.X.toarray()
        
    if not LVG:
        adata = adata[:, adata.var['Variance Type'] == 'HVG']
        
    return adata

def convert_string_to_encoding(string, vector_key):
    """A function to convert a string to a numeric encoding"""
    return np.argwhere(vector_key == string)[0][0]

def convert_vector_to_encoding(vector, printkey = False):
    """A function to convert a vector of strings to a dense numeric encoding"""
    vector_key = np.unique(vector)
    if printkey:
        print(pd.Series(vector_key))
    
    vector_strings = list(vector)
    vector_num = [convert_string_to_encoding(string, vector_key) for string in vector_strings]
    
    return vector_num