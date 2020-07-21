import os
import pandas as pd
import subprocess
from _helper_deconv import  get_gep_eset_file, get_dir_input, get_dir_output


def music(gep, bulk_rna, out_dir):
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
        raise FileNotFoundError(
            'Gene expression profile '+gep+' cannot be found.')
    if not os.path.isfile(os.path.join(bulk_rna, 'input.csv')):
        raise FileNotFoundError(
            'Input file: input.csv cannot be found in the input folder')
    cmd = 'Rscript /app/script_music.R ' + gep + ' ' + bulk_rna + ' ' + out_dir
    p = subprocess.run([cmd], shell=True, stdout=None, stderr=None)

    print('Deconvolution using MuSiC finished')
    return (None)


file_annot = get_gep_eset_file()
dir_input = get_dir_input()
dir_output = get_dir_output()

pred = music(gep=file_annot, bulk_rna=dir_input, out_dir=dir_output)
#export_predictions(pred_df=pred, header=True, index=True)
