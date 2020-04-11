BESCAPE - BESCA Proportion Estimator
=====================================

BESCAPE is a cell deconvolution package. The user can specify a custom basis vector, as well as the preferred deconvolution method. Thus it allows us to detach the deconvolution algorithm from the underlying basis vector it originally comes packaged with. 

The module distinguishes between two types of basis vectors as input:
1. Gene Expression Profile (GEP) - generated from single-cell annotations using BESCA.export functions 
2. Single-cell annotation AnnData object - should contain single-cell annotations of multiple samples from which the deconvolution method generates its own GEP

Currently supported deconvolution methods:
* bescape - in-house method based on nu-SVR (CIBERSORT)
* [EPIC](https://github.com/GfellerLab/EPIC)
* [MuSiC](https://github.com/xuranw/MuSiC)
* [SCDC](https://github.com/meichendong/SCDC/)

![summary fig][bescape summary]

[bescape summary]: https://github.com/bedapub/bescape/blob/master/docs/fig/bescape_summary_hires.png "BESCApe summary figure"

## Example Usage
A more detailed guide on example use and integration with [BESCA](https://github.com/bedapub/besca) are documented in the following tutorial: <https://bedapub.github.io/besca/tutorials/bescape.html>

```python
from bescape import Bescape

# Initiate Bescape object with docker set as service and specifying an image to use
# (can be local or on DockerHub)
deconv = Bescape(service='docker', docker_image='phanmir/bescape:0.4')

# deconvolute using MuSiC - single-cell-based basis vector
deconv.deconvolute_sc(dir_annot='./datasets/music/gep/', 
                      dir_input='./datasets/music/input',
                      dir_output='./datasets/music/output', 
                      method='music')
```

## Installation
Install Bescape using `pip`:

```sh
pip install bescape
```

Bescape requires Docker or Singularity service to run. Links to installation instructions:
* [Docker][docker]
* [Singularity][singularity]

[docker]: https://docs.docker.com/get-docker/
[singularity]: https://sylabs.io/guides/3.0/user-guide/installation.html

###
