\name{GetGLM}
\alias{GetGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
GLM model for each gene
}
\description{
The \code{GetGLM} function adjusts a generalized linear model (GLM) for each of the genes, so that it can be analyzed which regulators and experimental variables have a significant effect. To get to the final model, variable selection methods are applied (lasso, ridge regression and ElasticNet), stepwise procedure, filter regulators with low variation, missing values and correlated regulators.
}
\usage{
GetGLM(GeneExpression, associations, data.omics, edesign = NULL,
center = TRUE, scale = FALSE, Res.df = 5, epsilon = 0.00001,
alfa = 0.05, MT.adjust = "none", family = negative.binomial(theta=10),
elasticnet = 0.5, stepwise = "backward", interactions.reg = TRUE,
min.variation = 0, correlation = 0.9, min.obs = 10, omic.type = 0)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{GeneExpression}{Matrix or data frame containing gene expression data with genes in rows and experimental samples in columns. The row names must be the gene IDs.}

  \item{associations}{List where each element corresponds to a different omic data type (miRNAs, transcription factors, methylation, etc.). The names of the list will be the omics. Each element is a data frame with two columns (optionally three) describing the potential interactions between genes and regulators for that omic. First column must contain the regulators, second the gene IDs, and an additional column can be added to describe the type of interaction (for example, in methylation data, if a CpG site is located in the promoter region of the gene, in the first exon, etc.).}

  \item{data.omics}{List where each element corresponds to a different omic data type (miRNAs, transcription factors, methylation, etc.). The names of this list will be the omics and each element of the list is a  matrix or data frame with omic regulators in rows and samples in columns.}

  \item{edesign}{Data frame or matrix describing the experimental design. Rows must be the samples, that is, the columns in the \code{GeneExpression}, and columns must be the experimental variables to be included in the model, such as disease, treatment, etc.}

  \item{center}{If TRUE (default), the omic data are centered.}

  \item{scale}{If TRUE, the omic data are scaled. Default value is set to FALSE.}

  \item{Res.df}{Number of degrees of freedom in the residuals. By default, 5. Increasing \code{Res.df} will increase the power of the statistical model and decrease the number of significant predictors.}

  \item{epsilon}{A threshold for the positive convergence tolerance in the GLM model. By default, 0.00001.}

  \item{alfa}{Significance level. By default, 0.05.}

  \item{MT.adjust}{Multiple testing correction method to be used within the stepwise variable selection procedure. By default, "none". You can see the different options in \code{?p.adjust}.}

  \item{family}{Error distribution and link function to be used in the model (see glm for more information). By default, \code{negative.binomial(theta = 10)}.}

  \item{elasticnet}{ElasticNet mixing parameter. Its values can be the following: NULL (No ElasticNet variable selection), value between 0 and 1 (ElasticNet is applied with this number being the combination between ridge and lasso penalization), where value 0 is the ridge penalty and value 1 is the lasso penalty. By default, 0.5.}

  \item{stepwise}{Stepwise variable selection method to be applied. It can be one of: "none", "backward" (by default), "forward", "two.ways.backward" or "two.ways.forward".}

  \item{interactions.reg}{If TRUE (default), MORE allows for interactions between each regulator and the experimental covariate.}

  \item{min.variation}{For numerical regulators, the minimum change that a regulator must present across conditions to keep it in the regression models. For binary regulators, if the proportion of the most repeated value equals or exceeds this value, the regulator will be considered to have low variation and removed from the regression models. By default, its value is 0.}

  \item{correlation}{Correlation threshold (in absolute value) to decide which regulators are correlated, in which case, a representative of the group of correlated regulators is chosen to enter the model. By default, 0.9.}

  \item{min.obs}{Minimun number of observations a gene must have to compute the GLM model. By default, 10.}

  \item{omic.type}{Vector with as many elements as the number on omics, which indicates if the omic values are numeric (0, default) or binary (1). When a single value is given, the type for all the omics is set to that value. By default, 0.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
\code{ResultsPerGene} It is a list that contains the following objects for each gene:
  \item{Y}{Data frame with the response variable values for that gene (y), the values fitted by the model (fitted.y), and the residuals of the GLM (residuals).}

  \item{X}{Data frame with all the predictors included in the final model.}

  \item{coefficients}{Matrix with the estimated coefficients for significant regulators and their p-values.}

  \item{allRegulators}{Data frame with all the initial potential regulators in rows and the following information in columns: gene, regulator, omic, area (the third optional column in associations), filter (if the regulator has been filtered out of the model, this column indicates the reason), and Sig (1 if the regulator is significant and 0 if not). Regarding the filter column, several values are possible:
  \itemize{
    \item MissingValue. If the regulator has been filtered out of the study because it has missing values.
    \item LowVariation. If the regulator has been filtered out of the study because it has lower variability than the threshold set by the user in min.variation parameter.
    \item ElasticNet. If the regulator has been excluded of the model by the ElasticNet variable selection procedure.
    \item Model. When the regulator is included in the initial equation model.
    \item omic_mcX_X. This notation is related to highly correlated regulators and how they are treated to avoid the multicollinearity problem.
  }}

  \item{significantRegulators}{A character vector containing the significant regulators.}

\code{GlobalSummary} List that contains the following elements:
  \item{GoodnessOfFit}{Matrix that collects the model p-value, the final number of degrees of freedom in the residuals, the R-squared value (which for GLMs is defined as the percentage of deviance explained by the model) and the Akaike information criterion (AIC) value for all the genes studied.
}

  \item{ReguPerGene}{Matrix containing, for each omic and gene, the number of initial regulators, the number of regulators included in the initial model and the number of significant regulators.}

  \item{GenesNOmodel}{List of genes for which the final GLM model could not be obtained. There are three possible reasons for that and they are indicated: "Too many missing values", "No predictors after EN" and "GLM error", where EN refers to the ElasticNet variable selection procedure.}

\code{Arguments} List with the arguments used to generate the model: experimental design matrix, minimum degrees of freedom in the residuals, significance level, family distribution, variable selection, etc.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Sonia Tarazona, Blanca Tomas Riquelme
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
\code{\link{TestData}}
}

\examples{
data(TestData)
require(MASS);
require(igraph);
require(car);
require(glmnet);
require(psych)

## Omic type

OmicType = c(1, 0, 0)
names(OmicType) = names(TestData$data.omics)

## GetGLM function

SimGLM = GetGLM(GeneExpression = TestData$GeneExpressionDE,
                associations = TestData$associations,
                data.omics = TestData$data.omics,
                edesign = TestData$edesign[,-1, drop = FALSE],
                Res.df = 10,
                epsilon = 0.00001,
                alfa = 0.05,
                MT.adjust = "fdr",
                family = negative.binomial(theta = 10),
                elasticnet = 0.7,
                stepwise = "two.ways.backward",
                interactions.reg = TRUE,
                correlation = 0.8,
                min.variation = 0,
                min.obs = 10,
                omic.type = OmicType)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
