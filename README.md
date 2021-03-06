BESCAPE - BESCA Proportion Estimator
=====================================

BESCAPE is a cell deconvolution package. The user can specify a custom basis vector, as well as the preferred deconvolution method. Thus it allows us to detach the deconvolution algorithm from the underlying basis vector it originally comes packaged with. 

Currently supported deconvolution methods:
* [MuSiC](https://github.com/xuranw/MuSiC)
* [SCDC](https://github.com/meichendong/SCDC/)

![summary fig][bescape summary]

[bescape summary]: https://github.com/bedapub/bescape/blob/master/docs/fig/bescape_summary_hires.png "BESCApe summary figure"

## Example Usage
A more detailed guide on example use and integration with [BESCA](https://github.com/bedapub/besca) is available in a [notebook here](docs/bescape_tutorial.ipynb). A tutorial for converting AnnData to ExpressionSet is available [here](https://github.com/bedapub/bescape/blob/master/docs/adata_to_eset.ipynb).

```python
import os
from bescape import Bescape

# docker
deconv = Bescape(service='docker', docker_image='bedapub/bescape:latest')

# Important to specify ABSOLUTE directory paths
wd = os.getcwd()
annot = wd + '/datasets/music/gep/'
inpt = wd + '/datasets/music/input'
output = wd + '/datasets/music/output'
# deconvolute using MuSiC - sc based basis vector
deconv.deconvolute_sc(dir_annot= annot, 
                      dir_input= inpt,
                      dir_output= output, 
                      method='music')
```

## Installation
Install Bescape using `pip`:

```sh
pip install bescape
```

Bescape requires Docker or Singularity service to run. Links to installation instructions:
* [Docker][docker] - once installed make sure Docker is running before using the Bescape package. Either open the Docker Desktop app or run `sudo dockerd` in terminal
* [Singularity][singularity]

[docker]: https://docs.docker.com/get-docker/
[singularity]: https://sylabs.io/guides/3.0/user-guide/installation.html


###
