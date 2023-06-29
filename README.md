# MORE

MORE (Multi-Omics REgulation) is an R package for the application of Generalized Linear Models (GLM) with Elastic Net 
regularization or Partial Least Squares (PLS) to multi-omics data. The MORE method applies GLMs or PLS to model gene expression 
as a function of experimental variables, such as diseases or treatments, and the potential regulators of a given gene.
The aim is to obtain specific candidate regulators for the biological system under study.


## Installing

Currently, the package can be installed directly from GitHub using the `devtools` R package:

    install.packages("devtools")
    devtools::install_github("ConesaLab/MORE")

Before installation, it might be necessary to install the required dependencies:

* pbapply
* glmnet
* igraph
* MASS
* parallel
* psych
* car
* ltm
* ropls
* fastDummies

## Usage

You can find the User Guide for the package in the vignettes folder or access it directly [here](https://github.com/ConesaLab/MORE/blob/master/vignettes/UserGuide.pdf).

MORE can also be run as a Shiny application. [In this directory](https://bitbucket.org/ConesaLab/more/downloads/) 
you will find .zip file containing all the required inputs to run it.

You can find more information on how to run MORE as a Shiny app in the User Guide,
and also in the following [video](https://youtu.be/SSIaeFRNsXg).

