import os
from tempfile import TemporaryDirectory
from shutil import copy
from spython.main import Client

dirname = os.path.dirname(__file__)

PATH_SINGULARITY = os.path.join(dirname, './bescape_singularity.sif')
PATH_DOCKER = os.path.join(dirname, '')


class Bescape:
    img_dir_input = ':/app/input'
    img_dir_output = ':/app/output'
    img_dir_gep = ':/app/gep'

    def __init__(self):
        Client.load(PATH_SINGULARITY)

    def deconvolute_gep(self,
                        file_gep,
                        dir_input,
                        dir_output,
                        method='bescape'):
        """Perform deconvolution using a Gene expression profile as the basis vector
        Calls a Singularity/Docker container which contains the required Python and R dependencies 
        for all supported deconvolution methods. 

        Args:
            file_gep (.csv): path to the gene expression profile .csv file. First column has to be 
                ENSEMBL gene names, first row are cell type labels. Gene expression should be linear, 
                not log-transformed. 
            dir_input (String): absolute path to the directory  containing the input.csv file. 
                The correct format of input.csv is described in the tutorial.
            dir__output (String): absolute path for desired output directory
            method (String): desired method to be used for deconvolution. 
                Current options: ['bescape', 'epic', 'music', 'scdc']

        Returns:
            A predictions.csv files with cell type proportions
        """

        if method == 'bescape':
            with TemporaryDirectory() as tmp_gep_dir:
                copy(file_gep, tmp_gep_dir)
                dir_gep = os.path.join(tmp_gep_dir, self.img_dir_gep)
                dir_input = os.path.join(dir_input, self.img_dir_input)
                dir_output = os.path.join(dir_output, self.img_dir_output)
                param_bind = [dir_input, dir_output, dir_gep]
                cmd = ['python3', '/app/bescape.py']
                Client.execute(cmd, bind=param_bind)
        else:
            raise Exception('Selected method not supported: ', method)
        return None

    def deconvolute_sc(self,
                       dir_annot,
                       dir_input,
                       dir_output,
                       method='music'):
        """Perform deconvolution using single cell (sc) annotations as the basis vector.
        User provides single cell annotations (either as scanpy AnnData obj or as eset RDS obj). 
        The deconvolution method will generate its own GEP based on the provided annotations.

        Args:
            dir_annot (String) && method=='music': absolute path to the directory  containing only 
                ONE cell annotation file. This should either be an ExpressionSet class saved as .RDS, 
                or an AnnData object saved as .h5ad. If AnnData is provided, it will be converted to 
                an ExpressionSet.
            dir_annot (String) && method=='scdc': absolute path to the directory  containing MULTIPLE
                cell annotation files. These should either be an ExpressionSet class saved as .RDS, 
                or an AnnData object saved as .h5ad. If AnnData is provided, it will be converted to 
                an ExpressionSet.

        """

        if method == 'music':
                dir_gep = os.path.join(dir_annot, self.img_dir_gep)
                dir_input = os.path.join(dir_input, self.img_dir_input)
                dir_output = os.path.join(dir_output, self.img_dir_output)
                param_bind = [dir_input, dir_output, dir_gep]
                cmd = ['python3', '/app/musicpy.py']
                Client.execute(cmd, bind=param_bind)
        elif method == 'scdc':
                dir_gep = os.path.join(dir_annot, self.img_dir_gep)
                dir_input = os.path.join(dir_input, self.img_dir_input)
                dir_output = os.path.join(dir_output, self.img_dir_output)
                param_bind = [dir_input, dir_output, dir_gep]
                cmd = ['python3', '/app/scdc.py']
                Client.execute(cmd, bind=param_bind)
        else:
            raise Exception('Selected method not supported: ', method)
