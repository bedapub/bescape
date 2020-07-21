import os
from tempfile import TemporaryDirectory
from shutil import copy


dirname = os.path.dirname(__file__)

PATH_SINGULARITY = os.path.join(dirname, '/bescape_singularity.sif')
PATH_DOCKER = os.path.join(dirname, '')


class Bescape:
    """BESCApe - a cell deconvolution class

    The deconvolution classifier is packed in a container. This class is a wrapper
    that uses either Docker or Singularity as a service to run containers.

    Args:
        service (str): One of ['docker', 'singularity']. Determine which service to use
        docker_image (str): Only for Docker. Provide the URI to the BESCApe deconvolution image on DockerHub. It will first
            search for a local image of same name. If not present, it will pull the image from Dockerhub. Thus, one can build
            a local image from the Bescape source, modify it, and run it in this Bescape class.
        path_singularity (str): Only for Singularity. Provide the absolute path to the Singularity container file for BESCApe.
            If not specified, Bescape will attempt to pull a docker image from DockerHub and build a singularity container from there.
    """
    img_dir_input = '/app/input'
    img_dir_output = '/app/output'
    img_dir_gep = '/app/gep'

    def __init__(self, service='docker', docker_image='bedapub/bescape:latest', path_singularity=None):
        self.service = service
        self.docker_image = docker_image
        self.path_singularity = path_singularity

        if self.service == 'docker':
            try:
                import docker
                self.client = docker.from_env()
                print('Docker client instantiated')
                
                try:
                    self.client.images.get(docker_image)
                    print('Local Docker image found: ', docker_image)
                except docker.errors.ImageNotFound:
                    print('No local Docker image found. Attempting to pull from Dockerhub.')
                    self.client.images.pull(docker_image)
                    print('Docker image loaded: ', docker_image)

            except ModuleNotFoundError:
                print('Python library for the Docker Engine API not installed')
                    
        elif self.service == 'singularity':
            try:
                from spython.main import Client
                self.sclient = Client
                if path_singularity is not None:
                    self.sclient.load(self.path_singularity)
                    if self.sclient is None:
                        raise FileNotFoundError("Singularity container not found. Modify path or set path_singularity to None and specify docker_image to pull and build a singularity container")
                else:
                    dockerhub_image_uri = os.path.join('docker://', docker_image)
                    self.sclient = Client.pull(image=dockerhub_image_uri, name='bescape.sif')
                    if self.sclient is None:
                        raise ValueError("Could not pull the docker image. Singularity container was not created.")
                
                print('Singularity client loaded')
                print('Singularity container loaded: ', self.path_singularity)

            except ModuleNotFoundError:
                print('Python library for the Singularity API not installed')
        else:
            raise ValueError(
                "Selected platform not supported. Chooes either Docker or Singularity")


    def deconvolute_sc(self,
                       dir_annot,
                       dir_input,
                       dir_output,
                       method='music',
                       **kwargs):
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

        if not os.path.exists(dir_output):
            os.makedirs(dir_output)
        if self.service == 'docker':
            vol_dict = {dir_annot.rstrip('/'): {'bind': self.img_dir_gep, 'mode': 'ro'},
                        dir_input.rstrip('/'): {'bind': self.img_dir_input, 'mode': 'ro'},
                        dir_output.rstrip('/'): {'bind': self.img_dir_output, 'mode': 'rw'}}

            if method == 'music':

                c = self.client.containers.run(
                    self.docker_image, command='python3 musicpy.py', volumes=vol_dict, detach=True)
                for line in c.logs(stream=True):
                    print(line.decode().strip())

            elif method == 'scdc':
                celltypevar = 'cluster' #kwargs.get("celltype_var")
                celltypesel = kwargs.get("celltype_sel")
                celltypesel = ["^".join(i.split()) for i in celltypesel]
                samplevar = 'SubjectName' #kwargs.get("sample_var")
                cmd = 'python3 scdcpy.py --celltypevar ' +  celltypevar +  ' --samplevar ' +  samplevar + ' --celltypesel ' + celltypesel
                c = self.client.containers.run(
                    self.docker_image, command=cmd, volumes=vol_dict, detach=True)
                for line in c.logs(stream=True):
                    print(line.decode().strip())
            else:
                raise ValueError("Selected method not supported: %s" % method)

        elif self.service == 'singularity':
            if method == 'music':
                dir_gep = os.path.join(dir_annot.rstrip(
                    '/'), ':/', self.img_dir_gep.lstrip('/'))
                dir_input = os.path.join(dir_input.rstrip(
                    '/'), ':/', self.img_dir_input.lstrip('/'))
                dir_output = os.path.join(dir_output.rstrip(
                    '/'), ':/', self.img_dir_output.lstrip('/'))
                param_bind = [dir_input, dir_output, dir_gep]
                cmd = ['python3', '/app/musicpy.py']
                executor = self.sclient.execute(cmd, bind=param_bind, stream=True)
                for line in executor:
                    print(line)

                
            elif method == 'scdc':
                dir_gep = os.path.join(dir_annot, self.img_dir_gep)
                dir_input = os.path.join(dir_input, self.img_dir_input)
                dir_output = os.path.join(dir_output, self.img_dir_output)
                param_bind = [dir_input, dir_output, dir_gep]
                cmd = ['python3', '/app/scdc.py']
                self.sclient.execute(cmd, bind=param_bind)
                
        else:
            raise Exception('Selected method not supported: ', method)
