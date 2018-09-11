\name{TestData}
\alias{TestData}
\docType{data}
\title{
Omic Data
}
\description{
This file is composed of four objects and can be used to test the MORE algorithm: gene expression matrix, association matrix, data omics matrix and the experimental design matrix.
}
\usage{data(TestData)}
\format{
  \code{TestData} contains the following four lists:
  \describe{
    \item{\code{GeneExpressionDE}}{It contains the counts of gene expression of each gene for the omic RNA-seq.}
    \item{\code{data.omics}}{It contains three data matrices with the expression values (count values) for each considered omic (ChIP-seq, miRNA-seq, TF).}
    \item{\code{edesign}}{It is the experimental design matrix and it is composed of 8 time points in two conditions, therefore, we ave 16 experimental samples.}
    \item{\code{associations}}{It contains three association matrices for each omic (potential regulators).}
}
}
\details{
The count values for ChIP-seq are categorize, therefore, there are value 1 or 0 (1 if the TF analysed in the ChIP-seq experiment binds to the gene and 0 if not). However, the count values for miRNA-seq and TF omics are not categorize.
}
\source{
\url{https://bitbucket.org/ConesaLab/mosim/}
}
%\references{
%}
\keyword{datasets}
