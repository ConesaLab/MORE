#########################################################################################
######           Functions to integrate omics data using GLMs                      ######
#########################################################################################


## By Sonia & Monica
## 07-July-2016
## Last modified: August 2022


options(stringsAsFactors = FALSE)

# library(maSigPro)  ### ver que pasa si no cargamos maSigPro
library(igraph)
library(MASS)
library(glmnet)
library(psych)
library(car)
library(ltm)
library(ppcor)
library(Hmisc)
library(ropls)
library(fastDummies)

suppressMessages(source("auxFunctions.R"))
suppressMessages(source("MORE_GLM.R"))
suppressMessages(source("ComputeGLM_function.R"))
# Generalized Lineal Model  -----------------------------------------------
# DETAILS
# In the case that the data matrix contains more variables than samples, a stepwise forward is applied, but the experimental
# design variables are kept.
#
# VALUE
# A list with 3 objects:
# - Arguments
# - SummaryTable: A data.frame with one row for each gene supplied for the user. The data.frame has 4 columns:
#   - gene: gene IDs
#   - RP: A tag with 3 possibles values. If the gene is classified corretly in a regulatory program, the RP name,
#     if the gene is wrong classified, "WC", and if the gene is not classified in any of the regulatory program, "NC".
# - FinalResults: A list with one element for each gene studied. Each element of the list is also a list containing:
#  - GLMfinal: an object of class glm  with the best model obtained
#  - GLMorigen: an object of class glm with the starting model if we have enough df for performing at generalized linear models.
#    In other case, GLMfinal=GLMorigen
#  - SummaryStepwise: A named (by methologies stepwise applied) list. The list names could be one or some of the following
#    "stepfor","two.ways.stepfor" or "two.ways.stepback“. It depends on "stepwise" selection. Each one of these objects contain
#    a list with some summary parameters, "sol","coefficients","t.score","variables" and "edesign"
#    - sol, matrix for summary results of the stepwise regression. The following values are given:
#    p-value of the regression ANOVA, R-squared of the model, AIC of the model and p-value of the regression coefficients
#    of the selected variables.
#    - coefficients, regression coefficients for significant regulators
#    - t.score, value of the t statistics of significant regulators
#    - variables: variables in the final model
#    - edesign: matrix of experimental design. ¿LO QUITAMOS DE AQUI?
#  - RegulOrig: a list with 3 objects, with some summary information about original regulators in the model.
#    - OrigRegulators: The start regulators of the model (previous to make any stepwise). Ponemos los mismos que en RegulatorsValue??
#    - RegulatorsValue: The start regulators value of the model after applying CleanPredictors (previous to make any stepwise)
#    - edesign: matrix of experimental design
#    - cor.max: Correlation value to decide when multicollinearity is present.
#    - action: Action to perform when regulators are correlated. If "mean" (default), the average of the correlated regulators is
#              computed. If "random", one of the correlated regulators is randomly selected and the others are discarded.

#' Genome-wide Generalized Linear Models
#'
#'\code{more} fits a regression model (if the selected method is GLM) or a PLS (if the selected method is PLS) for all the genes in the data set to identify
#' the experimental variables and potential regulators that show a significant effect on
#' the expression of each gene.
#'
#' @param GeneExpression Data frame containing gene expression data with genes in rows and
#' experimental samples in columns. Row names must be the gene IDs.
#' @param associations List where each element correspond to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will be the omics. Each element
#' is a data frame with 2 columns (optionally 3) describing the potential interactions between genes
#' and regulators for that omic. First column must contain the genes (or features in
#' GeneExpression object), second column must contain the regulators, and an additional column can
#' be added to describe the type of interaction (e.g., for methylation, if a CpG site is located in
#' the promoter region of the gene, in the first exon, etc.).
#   (optional).
#' @param data.omics List where each element correspond to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will be the omics. Each element
#' is a data matrix with omic regulators in rows and samples in columns.
#' @param edesign Data frame describing the experimental design. Rows must be the samples (columns
#' in GeneExpression) and columns must be the experimental variables to be included in the model
#' (e.g. time, treatment, etc.).
#' @param cont.var Name of the column in edesign that is to be considered as a continuous
#' explanatory variable in the regression model. NULL if edesign does not contain any continuous
#' variables. By default, "Time".
#' @param Res.df Number of degrees of freedom in the residuals. By default, 5. Increasing
#' Res.df will increase the power of the statistical model.
#' @param alfa Significance level. By default, 0.05.
#' @param MT.adjust Multiple testing correction method to be used within the stepwise variable
#' selection procedure. By default, "none". See the different options in ?\code{p.adjust}.
#' @param interactions.exp If TRUE (default), interactions among the experimental variables
#' are included in the model.
#' @param interactions.reg If TRUE (default), interactions between regulators and experimental
#' variables are included in the model.
#' @param min.variation  For numerical regulators, the minimum change that a regulator must present across conditions
#' to keep it in the regression models. For binary regulators, if the proportion of the most repeated value equals or exceeds this value,
#' the regulator will be considered to have low variation and removed from the regression models.
#' @param correlation
#' @param min.obs
#' @param elasticnet  ElasticNet mixing parameter. NULL = No ElasticNet variable selection. A number between 0 and 1 = ElasticNet is applied
#' with this number being the combination between Ridge and Lasso penalization (elasticnet=0 is the ridge penalty, elasticnet=1 is the lasso
#' penalty). The default is 0.5, which gives equal weight to ridge and lasso.
#' @param omic.type 0 = numerical variable (default), 1 = binary (categorical) variable. It should be a vector with length the number of omic types or a number (0 or 1) if all the omics are of the same type.
#' @param edesign.type 0 = numerical variable (default), 1 = binary (categorical) variable. It should be a vector with length the number of experimental design variables or a number (0 or 1) if all the experimental design variables are of the same type.
#' @param method 'GLM' = Apply a GLM with ElasticNet regularization. 'PLS' = Apply PLS model
#' @param filter '1' = Applies MORE original multicollinearity filter method. '2' = Applies MORE multicollinearity filter method based on correlations between different omics (COR). '3' = Applies multicollinearity filter method based on Partial Correlation (PC).
#' @return
#' @export
#'
#' @examples

more <-function(GeneExpression,
                associations,
                data.omics,
                edesign = NULL,
                center = TRUE, scale = FALSE,
                Res.df = 5, epsilon = 0.00001,
                alfa = 0.05, MT.adjust = "none",
                family = negative.binomial(theta=10),
                elasticnet = 0.5,
                interactions.reg = TRUE,
                min.variation = 0,
                correlation = 0.9,
                min.obs = 10,
                omic.type = 0,
                filtermet = 1,
                scaletype = 1,
                edesign.type= 0,
                p.method = 'jack',
                all_settings,
                method  ='GLM'){
  
  if(method=='GLM'){
    
    return(GetGLM(GeneExpression,
                  associations,
                  data.omics,
                  edesign = edesign,
                  center = center, scale = scale,
                  Res.df = Res.df, epsilon = epsilon,
                  alfa = alfa, MT.adjust = MT.adjust,
                  family = family,
                  elasticnet = elasticnet,
                  interactions.reg = interactions.reg,
                  min.variation = min.variation,
                  correlation = correlation,
                  min.obs = min.obs,
                  omic.type = omic.type,
                  filtermet = filtermet,
                  all_settings))
    
  }
  
  if (method == 'PLS'){
    
    return(GetPLS(GeneExpression,
           associations,
           data.omics,
           edesign = edesign,
           center = center, scale = scale,
           epsilon = epsilon,
           alfa = alfa, MT.adjust = MT.adjust,
           interactions.reg = interactions.reg,
           min.variation = min.variation,
           min.obs = min.obs,
           omic.type = omic.type,
           edesign.type = edesign.type,
           scaletype = scaletype,
           p.method =p.method,
           all_settings))
  }
  
}
