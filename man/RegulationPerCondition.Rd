\name{RegulationPerCondition}
\alias{RegulationPerCondition}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summary table
}
\description{
The function \code{RegulationPerCondition} returns a summary table that contains for each gene the regulators that are significant for each GLM model. Moreover, it shows the regression coefficients of the model of each regulator in the experimental covariates and, in the case of having correlated regulators, which regulator is the representative of each group of correlated regulators.
}
\usage{
RegulationPerCondition(getGLMoutput, betaTest = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{getGLMoutput}{It is the name of the object where we have saved the execution of the function \code{getGLM}.}

\item{betaTest}{If TRUE (default), a hypothesis test is perfomed on the coefficients of the final GLM model. If the coefficients are significant, their values are displayed. If not, it is indicated with a 0.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
\value{
The output is a summary table that shows for each gene (gene column) the significant regulators (regulator column) with its omic (omic column), area (area column), which is representative (R) of the correlated group (representative column) and the coefficients of the final model GLM (end of the columns). When no information is provided in the representative column, it means that the regulator was not part of a correlated group of regulators. Regulators correlated positively with the representative will have the same coefficients (same sign) than the representative, while negatively correlated regulators will have the same coefficients of the representative but with opposite sign. The final columns correspond to the regression coefficients of each regulator for each experimental group.
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
\code{\link{GetGLM}}
\code{\link{TestData}}
}
\examples{
data(TestData)
require(MASS); require(igraph); require(car); require(glmnet); require(psych)

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

myresults = RegulationPerCondition(SimGLM)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
