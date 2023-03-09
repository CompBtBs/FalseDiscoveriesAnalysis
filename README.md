# Adjusting for false discoveries in constraint- and sampling-based differential metabolic flux analysis
Journal of Biomedical Informatics

## Authors
- Bruno G. Galuzzi <bruno.galuzzi@unimib.it> (a, c)
- Luca Milazzo <l.milazzo1@campus.unimib.it> (b)
- Chiara Damiani <chiara.damiani@unimib.it> (a, c)

(a) Department of Biotechnology and Biosciences, University of Milano-Bicocca, Piazza
dell’Ateneo Nuovo, 1, Milan, 20125, Italy

(b) Department of Informatics, Systems, and Communications, University of
Milano-Bicocca, Piazza dell’Ateneo Nuovo, 1, Milan, 20125, Italy

(c) SYSBIO Centre of Systems Biology/ ISBE.IT, Milan, Milan,

## Table of Contents
* [General Information](#general-information)
* [Technologies Used](#technologies-used)
* [Setup and usage](#setup)

## General Information

![alt text](https://github.com/CompBtBs/FalseDiscoveriesAnalysis/blob/main/image.png)

- Different samples of the very same feasible region of a metabolic network 
can produce different marginal flux distributions, with the risk of false discoveries
- For Hit-and-Run strategies, the thinning value has a higher impact on
false discoveries than the sample size.
- Hypothesis test on KL-divergence fully correct for false discoveries
- Sampling the corners of a feasible region with random functions is less
prone to false discoveries and produces marginal flux distributions dif-
ferent from the ones of Hit-and-Run strategies.


## Technologies Used
- JupyterLab - version 3.2.1
- Python - version 3.9.7
- Matlab - version 9.8.0.1396136
- R - version 4.2.2
- CobraToolBox - https://opencobra.github.io/cobratoolbox/stable/installation.html


## Setup and usage
This repository was built in such a manner that allows users to easily reproduce
the entire analysis pipeline. In particular, the code automatically handles
the folders structures and files by using relative paths. It is mandatory to 
follow the instructions reported at the beginning of each file in the "code"
folder in order to correctly execute all the analysis.

