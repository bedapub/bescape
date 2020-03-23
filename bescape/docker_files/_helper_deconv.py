import os
import tempfile
import pandas as pd
"""These static variables are set up for singularity/docker - the container has a defined /app folder that is the working directory for the deconvolution module contained in it"""
DIR_INPUT = "/app/input"
DIR_OUTPUT = "/app/output"
DIR_GEP = "/app/gep"

FILE_INPUT = "input.csv"
FILE_OUTPUT = "predictions.csv"


def load_input_batch():
    """Load input data from DIR_INPUT mounted inside the Docker/Singularity container

    """
    in_path = os.path.join(DIR_INPUT, FILE_INPUT)
    input_df = pd.read_csv(in_path)
    return input_df


def load_gep():
    """Load Gene expression profile

    Args:
        gep_file (.csv): path to the gene expression profile .csv file. 
            First column has to be ENSEMBL gene names, first row are cell type labels. 
            Gene expression should be linear, not log-transformed. 
 
    Returns:
        pandas.DataFrame with GEP

    """

    f = os.listdir(DIR_GEP)[0]
    path_gep = os.path.join(DIR_GEP, f)
    if os.path.isfile(path_gep):
        return pd.read_csv(path_gep, index_col=0)


def load_input_sc():
    """Load input for singe-cell based deconvolution methods ['music', 'scdc']"""
    f = os.listdir(DIR_INPUT)[0]
    file_input  = os.path.join(DIR_INPUT, f)
    return None


def get_input_eset_file():
    """Load bulk rna file path to be deconvoluted. This is for R based methods ['music', 'scdc'] - file as eset class"""
    f = os.listdir(DIR_INPUT)[0]
    path_input = os.path.join(DIR_INPUT, f)

    return path_input


def get_gep_eset_file():
    """Load eset single cell expression file path for R based deconvolution methods ['music', 'scdc'].
       It only passes the full path of the gep eset file to the R-script. """
    try:
        f = os.listdir(DIR_GEP)[0]
    except FileNotFoundError:
        print('GEP directory is empty, need to provide an annotation file')
        
    file_gep = os.path.join(DIR_GEP, f)

    return file_gep


def export_predictions(pred_df, filename=FILE_OUTPUT, header=True, index=False):

    out_path = os.path.join(DIR_OUTPUT, FILE_OUTPUT)
    pred_df.to_csv(out_path, header=header, index=index)


def get_dir_input():
    """Returns the static input directory path inside the Docker/Singularity container"""
    return DIR_INPUT
