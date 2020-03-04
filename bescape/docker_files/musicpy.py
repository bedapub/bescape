import os
import pandas as pd
import subprocess
from _helper_deconv import export_predictions, get_eset_gep_file, get_eset_gep_file

def music(gep, bulk_rna, sep='\t'):
    """R wrapper function. Uses subprocess to call shell commands form python. Calls script_music.R file
        to deconvolute. 

        The input files are R classes. The user will either have directly provided these, or they will have
        been generated from AnnData files by bescape conversion functions.
    Args:
        gep (eset): eset expression file (R class) containing single cell annotations that must come 
        from at least two different subjects for music to build its own GEP
        bulk_rna (eset): eset expression file (R class) with bulk RNA expressions. Can be one or multiple samples
        to be deconvoluted.
    """
    if not os.path.isfile(gep):
        raise FileNotFoundError('Gene expression profile '+gep+' cannot be found.')
    if not os.path.isfile(bulk_rna):
        raise FileNotFoundError('Bulk RNA file '+ bulk_rna +' cannot be found.')
    print('Starting deconvolution using MuSiC. This process may take a few moments')
    print('-----------------------------------------------------------------------')
    cmd = 'Rscript /app/script_music.R ' + gep + ' ' + bulk_rna + ' temp.txt'
    p = subprocess.run([cmd], shell=True, stdout=None, stderr=None)
    
    out = pd.read_csv('temp.txt',sep='\t',header=0,index_col=0)
    os.remove('temp.txt')

    print('Deconvolution using MuSiC finished successfully')

    return out


file_annot = get_gep_eset_file()
dir_input = get_input_eset_file()

pred = music(gep=file_annot, bulk_rna=dir_input)
export_predictions(pred_df=pred, header=True, index=True)
