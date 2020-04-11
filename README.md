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

The package requires a Singularity or Docker platfrom to run.

Details on installation, example use, and integration with BESCA are documented in the following tutorial: <https://bedapub.github.io/besca/tutorials/bescape.html>

![summary fig][bescape summary]

[bescape summary]: https://github.com/bedapub/bescape/blob/master/docs/fig/bescape_summary_hires.png "BESCApe summary figure"
