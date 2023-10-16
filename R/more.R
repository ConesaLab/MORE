#########################################################################################
######           MORE:      Multi-Omics REgulation                        ######
#########################################################################################


## By Sonia, Monica & Maider
## 15-October-2023
## Last modified: August 2023


options(stringsAsFactors = FALSE)

library(igraph)
library(MASS)
library(glmnet)
library(psych)
library(car)
library(stringr)
library(fastDummies)

#' more: Multi-Omics Regulation
#'
#' \code{more} fits a GLM regression model (when the selected method is GLM) or a PLS model (when the selected method is PLS) for all genes in the dataset to identify
#' the potential regulators that show a significant impact on gene expression under specific experimental conditions.
#'
#' @param GeneExpression Data frame containing gene expression data with genes in rows and
#' experimental samples in columns. Row names must be the gene IDs.
#' @param data.omics List where each element corresponds to a different omic data type to be considered (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will represent the omics, and each element in 
#' the list should be a data matrix with omic regulators in rows and samples in columns.
#' @param associations List where each element corresponds to a different omic data type (miRNAs,
#' transcription factors, methylation, etc.). The names of the list will represent the omics. Each element in 
#' the list should be a data frame with 2 columns (optionally 3), describing the potential interactions between genes
#' and regulators for that omic. First column must contain the genes (or features in
#' GeneExpression object), second column must contain the regulators, and an optional third column can
#' be added to describe the type of interaction (e.g., for methylation, if a CpG site is located in
#' the promoter region of the gene, in the first exon, etc.). If the user lacks prior knowledge of the potential regulators, they can set the parameter to NULL. 
#' In this case, all regulators in \code{\link{data.omics}} will be treated as potential regulators for all genes. In this case, for computational efficiency, it is recommended to use pls2 \code{\link{method}}.
#' Additionally, if the users have prior knowledge for certain omics and want to set other omics to NULL, they can do so.
#' @param edesign Data frame describing the experimental design. Rows must be the samples (columns
#' in \code{\link{GeneExpression}}) and columns must be the experimental variables to be included in the model (e.g. treatment, etc.).
#' @param clinic Data.frame with all clinical variables to consider,with samples in rows and variables in columns.
#' @param clinic.type Vector which indicates the type of data of variables introduced in \code{\link{clinic}}. The user should code as 0 numeric variables and as 1 categorical or binary variables. 
#' By default is set to NULL. In this case, the data type will be predicted automatically. However, the user must verify the prediction and manually input the vector if incorrect.
#' @param center By default TRUE. It determines whether centering is applied to \code{\link{data.omics}}.
#' @param scale By default TRUE. It determines whether scaling is applied to \code{\link{data.omics}}.
#' @param epsilon Convergence threshold for coordinate descent algorithm in elasticnet. Default value, 1e-5.
#' @param alfa Significance level for variable selection in pls1and pls2 \code{\link{method}}. By default, 0.05.
#' @param family Error distribution and link function to be used in the model when \code{\link{method}} glm. By default, gaussian().
#' @param elasticnet ElasticNet mixing parameter. There are three options:
#' - NULL : The parameter is selected from a grid of values ranging from 0 to 1 with 0.1 increments. The chosen value optimizes the mean cross-validated error when optimizing the lambda values.
#' - A number between 0 and 1 : ElasticNet is applied with this number being the combination between Ridge and Lasso penalization (elasticnet=0 is the ridge penalty, elasticnet=1 is the lasso penalty). 
#' - A vector with the mixing parameters to try. The one that optimizes the mean cross-validated error when optimizing the lambda values will be used.
#' By default, NULL.
#' @param interactions.reg If TRUE, the model includes interactions between regulators and experimental variables. By default, TRUE.
#' @param min.variation  For numerical regulators, it specifies the minimum change required across conditions to retain the regulator in 
#' the regression models. In the case of binary regulators, if the proportion of the most common value is equal to or inferior this value, 
#' the regulator is considered to have low variation and will be excluded from the regression models. The user has the option to set a single 
#' value to apply the same filter to all omics, provide a vector of the same length as omics if they want to specify different levels for each omics, 
#' or use 'NA' when they want to apply a minimum variation filter but are uncertain about the threshold. By default, 0.
#' @param min.obs Minimum number of observations required for a gene to be considered. Default, 10.
#' @param col.filter Type of correlation coefficients to use when applying the multicollinearity filter when glm \code{\link{method}} is used. 
#' - cor: Computes the correlation between omics. Pearson correlation between numeric variables, phi coefficient between numeric and binary and biserial correlation between binary variables. 
#' - pcor : Computes the partial correlation.
#' @param correlation  Value to determine the presence of collinearity between two regulators when using the glm \code{\link{method}}. By default, 0.7.
#' @param scaletype Type of scaling to be applied. Three options:
#' - auto : Applies the autoscaling. 
#' - pareto : Applies the pareto scaling. \[ \frac{X_k}{s_k \sqrt[4]{m_b}} \]
#' - block : Applies the block scaling. \[ \frac{X_k}{s_k \sqrt{m_b}} \]
#' considering \(m_b\) the number of variables of the block. By default, auto.
#' @param p.method Type of resampling method to apply for the p-value calculation when pls1 or pls2 \code{\link{method}}. Two options:
#' - jack : Applies Jack-Knife resampling technique.
#' - perm : Applies a resampling technique in which the response variable is permuted 100 times to obtain the distribution of the coefficients and compute then their associated p-value.
#' By default, jack.
#' @param vip Value of VIP above which a variable can be considered significant in addition to the computed p-value in \code{\link{p.method}}. By default, 0.8.
#' @param method Model to be fitted. Two options:
#' - glm : Applies a Generalized Linear Model (GLM) with ElasticNet regularization.
#' - pls1 : Applies a Partial Least Squares (PLS) model, one for each of the genes at \code{\link{GeneExpression}}.
#' - pls2 : Applies a PLS model to all genes at the same time, only possible when \code{\link{associations}}= NULL.
#' By default, glm.
#' @return List containing the following elements:
#' - ResultsPerGene : List with as many elements as genes in \code{\link{GeneExpression}}. For each gene, it includes information about gene values, considered variables, estimated coefficients,
#'                    detailed information about all regulators, and regulators identified as relevant (in glm scenario) or significant (in pls scenarios).
#' - GlobalSummary : List with information about the fitted models, including model metrics, information about regulators, genes without models, regulators, master regulators and hub genes.
#' - Arguments : List containing all the arguments used to generate the models.
#'
#' @examples
#' 
#' more(GeneExpression, associations, data.omics, center = TRUE, scale = TRUE, epsilon = 0.00001, family = gaussian(), elasticnet = NULL, interactions.reg = TRUE,
#'  min.variation = 0,  min.obs = 10, col.filter = 'cor', correlation = 0.7, method  ='glm')
#' 
#' @export

setClass("MORE")

more <-function(GeneExpression,
                associations,
                data.omics,
                edesign = NULL,
                clinic = NULL,
                center = TRUE, scale = FALSE,
                epsilon = 0.00001,
                alfa = 0.05, 
                family = gaussian(),
                elasticnet = 0.5,
                interactions.reg = TRUE,
                min.variation = 0,
                correlation = 0.9,
                min.obs = 10,
                omic.type = 0,
                col.filter = 'cor',
                scaletype = 'auto',
                clinic.type= 0,
                p.method = 'jack',
                vip = vip,
                method  ='glm'){
  
  if(method=='glm'){
    
    return(GetGLM(GeneExpression,
                  associations,
                  data.omics,
                  edesign = edesign,
                  center = center, scale = scale,
                  epsilon = epsilon,
                  alfa = alfa, 
                  family = family,
                  elasticnet = elasticnet,
                  interactions.reg = interactions.reg,
                  min.variation = min.variation,
                  correlation = correlation,
                  min.obs = min.obs,
                  omic.type = omic.type,
                  col.filter = col.filter))
    
  }
  
  if (method == 'pls'){
    
    return(GetPLS(GeneExpression,
           associations,
           data.omics,
           clinic = clinic,
           edesign = edesign,
           center = center, scale = scale,
           epsilon = epsilon,
           alfa = alfa, 
           interactions.reg = interactions.reg,
           min.variation = min.variation,
           min.obs = min.obs,
           omic.type = omic.type,
           clinic.type = clinic.type,
           scaletype = scaletype,
           p.method =p.method, 
           vip = vip))
  }
  
}
