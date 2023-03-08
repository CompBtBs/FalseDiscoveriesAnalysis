# Adjusting for false discoveries in constraint- and sampling-based differential metabolic flux analysis
## Bruno G. Galuzzi, Luca Milazzo, Chiara Damiani
>


## Table of Contents
* [General Information](#general-information)
* [Technologies Used](#technologies-used)
* [Setup and usage](#setup)
* [Project Status](#project-status)
* [Room for Improvement](#room-for-improvement)
* [Acknowledgements](#acknowledgements)
* [Contact](#contact)



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
- Jupyter - version x.x
- Python - version x.x
- Cobra - version x.x


## Setup and usage
This repository was built in such a manner that allows users to easily reproduce
the entire analysis pipeline. In particular, the code automatically handles
the folders structures and files by using relative paths. It is mandatory to 
follow the instructions reported at the beginning of each file in the "code"
folder in order to correctly execute all the analysis.



## Project Status
Project is: _complete_ 

## Room for Improvement


Room for improvement:


## Acknowledgements


## Contact
Created by 

