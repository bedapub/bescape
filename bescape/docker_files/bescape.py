#!/usr/bin/env python
# coding: utf-8

# import os
# import numpy as np
# import pandas as pd
# from sklearn.linear_model import LinearRegression
# from sklearn.model_selection import KFold
# from sklearn.linear_model import Ridge
# from sklearn.linear_model import Lasso
# from sklearn.model_selection import GridSearchCV
# from sklearn.model_selection import ShuffleSplit
# from sklearn.pipeline import Pipeline
# from sklearn.preprocessing import StandardScaler
# from sklearn.svm import NuSVR, LinearSVR
# from sklearn.ensemble import RandomForestRegressor
# 

# In[1]:


import os
import errno
import numpy as np
import pandas as pd
from sklearn.linear_model import LinearRegression, ElasticNet, Lasso,  BayesianRidge, LassoLarsIC
from sklearn.kernel_ridge import KernelRidge
from sklearn.model_selection import KFold, cross_val_score, train_test_split, ShuffleSplit,GridSearchCV
from sklearn.base import BaseEstimator, TransformerMixin, RegressorMixin, clone
from sklearn.pipeline import Pipeline, make_pipeline
from sklearn.preprocessing import StandardScaler, RobustScaler
from sklearn.metrics import mean_squared_error
from sklearn.svm import NuSVR, LinearSVR
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor


# In[2]:


from _helper_deconv import load_gep, load_input_batch, export_predictions, get_dir_input


# ### Load DREAM input and our GEPs
# Need to specify:
# * in_dir
# * gold_standards_file
# * gold_standards_dir
# * signature_path
# * drop['others']?

# In[8]:


input_df = load_input_batch()
gep = load_gep()


# In[9]:


dataset_names = input_df['dataset.name']
scales = input_df['scale']
c_types = input_df['cancer.type']
#native_probe = input_df['native.probe.type']
#expression_files = input_df['hugo.expr.file']
expression_files = input_df['ensg.expr.file']
expression_paths = []

for file in expression_files:
    dir_input = get_dir_input()
    expression_paths.append(os.path.join(dir_input, file))
    

def load_expression_file(expression_path):
    expression_df = pd.read_csv(expression_path, index_col = 'Gene')
    return expression_df


# In[10]:


# TODO
# deal with duplicate genes in sample bulk and in signature genes
# Note: do not merge samples - not all genes are intersecting across samples


# ### nu-SVR without CV:

# In[54]:



def intersect_genes(bulk_sample, signature_matrix):
    """Find intersecting subset of genes between sample and signature matrix"""
    #TODO add warning if only few genes are common
    
    # drop rows with nan
    bulk_sample_notna = bulk_sample[pd.notna(bulk_sample.iloc[:,0:1]).any(axis=1)]
    
    idx = bulk_sample_notna.index & signature_matrix.index
    return bulk_sample_notna.loc[idx], signature_matrix.loc[idx]


def build_model(bulk_sample_subset, signature_matrix_subset):
    """Create and fit a regression model to bulk RNA (y) and a signature matrix (X)
    
    y = w * X | w = 'predicted cell type fractions'
    
    
    Args:
        bulk_sample_subset (pandas.DataFrame): bulk RNA
        signatre_matrix_subset (pandas.DataFrame): GEP
    
    Returns:
        Instance of the regression model
    
    """
    
    pipe = Pipeline(steps=[
        ('scale', StandardScaler()),
        ('svr', NuSVR())
    ])
    
    cv = ShuffleSplit(test_size=0.01, n_splits=1) # no CV - we want to minimize training error and not generalization error
    parameters = {'svr__nu' : [0.25, 0.5, 0.75], 'svr__C' : [1e-2,1e-1,1], 'svr__kernel' : ['linear'], 'svr__verbose' : [True]}
    grid = GridSearchCV(pipe, param_grid = parameters, cv = cv, scoring = 'neg_mean_squared_error', verbose=10, n_jobs=-1)
    grid.fit(signature_matrix_subset, bulk_sample_subset.values.ravel())

    return grid

# following methods are temporary
def rm_neg(coefs):
    """Removes negative coefficients from the weight vector
    
    Args:
        coefs (numpy.array): coefficients from the regression model, representing cell fractions
    
    Returns:
        numpy.array with non-negative values
    
    """
    c = coefs.copy()
    c[c < 0] = 0 
    return c

def rm_small(coefs):
    """Remove small weight vector coefficents based on an arbitrary cutoff"""
    # improvement: use standard deviation
    c = coefs.copy()
    c[c < 0.00001] = 0
    return c

def rescale (coefs):
    """Rescale the weight vector so that the total sums up to 1"""
    c = coefs.copy()
    scale_factor = 1.0 / (np.sum(c))
    c *= scale_factor
    return c

def append_missing_celltype(df, dataset_name, sample, cell_type, prediction=0):
    """Appends missing cell type to the output
    
    Useful until we are able to integrate signatures with a full set of cell types
    
    Args:
        df (pandas.DataFrame): result df to be appended to
        dataset_name (string)
        sample (string)
        cell_type (string)
    
    Returns:
        pandas.DataFrame with appended celltype and its predicted proportion set to 0
    
    """
    
    res_df = df.copy()
    append_df = { "dataset.name" : dataset_name, 
                 "sample.id" : sample, 
                 "cell.type" : cell_type, 
                 "prediction" : prediction}
    res_df = res_df.append(append_df, ignore_index=True)
    return res_df

def linearize(bulk, scale):
    """Takes anti-log transformation on log scaled data
    
    Some data inputs from the DREAM challenge are log2 or log10 scaled and need to be linearized.
    
    Args:
        bulk (pandas.DataFrame): Bulk expression matrix
        scale (String): The scale of the expression data (i.e., Log2, Log10, Linear)
    Returns:
        Linearized bulk expression matrix
    
    """
    
    if scale == 'Log2':
        return 2 ** bulk
    elif scale == 'Log10':
        return 10 ** bulk
    else:
        return bulk
    
def to_coarse(result_df):
    """Converts fine grained cell type predictions to coarse grained. 
    
    This is done by simply summing up the fractions.
    
    Args:
        result_df (pandas.DataFrame): Result table with predicted fine grained cell type fractions as defined in the DREAM challenge
    Returns:
        pandas.DataFrame result table with coarse-grained cell type predictions as
        
    """
    
    dict_coarse = {'memory.B.cells' : 'B.cells', 
               'naive.B.cells' : 'B.cells',
               'memory.CD4.T.cells' : 'CD4.T.cells',
               'naive.CD4.T.cells' : 'CD4.T.cells',
               'regulatory.T.cells' : 'CD4.T.cells',
               'memory.CD8.T.cells' : 'CD8.T.cells',
               'naive.CD8.T.cells' : 'CD8.T.cells',
               'monocytes' : 'monolytic.lineage',
               'myeloid.dendritic.cells' : 'monolytic.lineage',
               'macrophages' : 'monolytic.lineage',
              }

    coarse_df = result_df.replace({'cell.type' : dict_coarse})
    coarse_df = coarse_df.groupby(['dataset.name', 'sample.id', 'cell.type']).sum().reset_index()
    return coarse_df
    


# In[55]:


# for nu-SVR
result_df = pd.DataFrame( columns = ['dataset.name', 'sample.id', 'cell.type', 'prediction'])
param_df = pd.DataFrame(columns =['dataset.name', 'sample.id'])
for dataset_name, c_type, scale, expression_path in list(zip(dataset_names, c_types, scales, expression_paths)):
    samples_df = load_expression_file(expression_path)
    samples_df = linearize(samples_df, scale)
    
   
    for sample in samples_df:
        print('Deconvoluting dataset {:6} and sample [name: {:3}] [{} out of {}]'.format(dataset_name, sample, samples_df.columns.get_loc(sample) + 1, samples_df.shape[1]))
        
        bulk_sample_subset, signature_matrix_subset = intersect_genes(samples_df[[sample]], gep)
        grid = build_model(bulk_sample_subset, signature_matrix_subset)
        estimator = grid.best_estimator_.named_steps['svr']
        fractions = estimator.coef_[0]
        
        fractions = rm_neg(fractions)
        fractions = rescale(fractions)
        
        out_df = pd.DataFrame( columns = ['dataset.name', 'sample.id', 'cell.type', 'prediction'])
        out_df['cell.type'] = signature_matrix_subset.columns
        out_df['prediction'] = fractions
        out_df['dataset.name'] = dataset_name
        out_df['sample.id'] = sample
        
        # TODO remove once integrated GEPs
        if c_type != 'BRCA':
            #out_df = append_missing_celltype(out_df, dataset_name, sample, cell_type="neutrophils")
            #mono = out_df.loc[out_df['cell.type']=='macrophages', 'prediction'].values + out_df.loc[out_df['cell.type']=='myeloid.dendritic.cells', 'prediction'].values
            #out_df = append_missing_celltype(out_df, dataset_name, sample, cell_type="monocytes", prediction = mono[0])
            pass
        
        
        
        result_df = result_df.append(out_df, ignore_index=True)
        export_predictions(result_df)

