import os
import argparse
import pandas as pd
import subprocess
from _helper_deconv import export_predictions, get_gep_eset_multiple, get_input_eset_file, get_dir_input, get_dir_output


def scdc(gep, bulk_rna, out_dir, sep='\t'):
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

    if not os.path.isfile(os.path.join(bulk_rna, 'input.csv')):
        raise FileNotFoundError(
            'Input file: input.csv cannot be found in the input folder')
    parser = argparse.ArgumentParser(description='internal parser for docker scdc module')
    parser.add_argument("--celltypevar", default='cluster')
    parser.add_argument("--celltypesel", nargs="+", default=["alpha","beta","delta","gamma","acinar","ductal"])
    parser.add_argument("--samplevar", default='sample')

    args = parser.parse_args()
    celltypevar = args.celltypevar
    samplevar = args.samplevar
    celltypesel_list = args.celltypesel
    celltypesel = ", ".join(":{0}:".format(i) for i in celltypesel_list) # result: c(:celltype_a:, :celltype_b:, ...) the colons will get replaced by single quotes (') in the R script. This silly workaround is because the python subprocess /bin/bash removes the single quotes.
    celltypesel = celltypesel.replace('^', ' ')
    
    cmd = ('Rscript /app/scdc_input.R '
           "'c(" + gep +")' " 
           + bulk_rna + ' '
           + out_dir
           + " 'c(" + celltypesel + ")' "
           + celltypevar + ' '
           + samplevar
    )

    p = subprocess.run([cmd], shell=True, stdout=None, stderr=None, executable='/bin/bash')

#    out = pd.read_csv('temp.txt', sep='\t', header=0, index_col=0)
#    os.remove('temp.txt')

    print('Deconvolution using SCDC finished')
    return (None)


file_annot = get_gep_eset_multiple()
dir_input = get_dir_input()
dir_output = get_dir_output()

pred = scdc(gep=file_annot, bulk_rna=dir_input, out_dir=dir_output)
#export_predictions(pred_df=pred, header=True, index=True)
