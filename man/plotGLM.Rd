\name{plotGLM}
\alias{plotGLM}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot MORE results
}
\description{
MORE includes the possibility of studying graphically the relationship of genes and regulators: given a specific gene and regulator, study the regulators for a given gene, analyze which genes are regulated by a specific regulator...
}
\usage{
plotGLM(GLMoutput, gene, regulator = NULL, reguValues = NULL,
plotPerOmic = FALSE, gene.col = 1, regu.col = NULL, order = TRUE,
xlab = "", cont.var = NULL, cond2plot = NULL,...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
 \item{GLMoutput}{It is the name of the object where we have saved the execution of the function \code{getGLM}.}

\item{gene}{ID of the gene to be plotted.}

\item{regulator}{ID of the regulator to be plotted. If NULL (default value), all the regulators of the gene are plotted.}

\item{reguValues}{Vector containing the values of a regulator that the user can optionally provide. If NULL (default value), these values are taken from \code{GLMoutput} as long as they are available.}

\item{plotPerOmic}{If TRUE, all the significant regulators of the given gene and the same omic are plotted in the same graph. If FALSE (default value), each regulator is plotted in a separate plot.}

\item{gene.col}{Color to plot the gene. By default, 1 (black).}

\item{regu.col}{Color to plot the regulator. If NULL (default), a color will be assigned by the function, that will be different for each regulatory omic.}

\item{order}{If TRUE (default), the values in X axis are ordered.}

\item{xlab}{Label for X axis.}

\item{cont.var}{Vector with length equal to the number of observations in data, which optionally may contain the values of the numerical variable (e.g. time) to be plotted in X axis. By default, NULL.}

\item{cond2plot}{Vector or factor indicating the experimental group of each value to represent. If NULL (default), the labels are taken from the experimental design matrix.}
}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
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

## Given a gene and regulator
plotGLM(GLMoutput = SimGLM,
          gene = "ENSMUSG00000012535",
          regulator = "6_29609160_29609425",
          plotPerOmic = FALSE,
          gene.col = "red4",
          regu.col = "green4")

## Given a gene
par(mfrow = c(2,2))
plotGLM(GLMoutput = SimGLM,
          gene = "ENSMUSG00000012535",
          regulator = NULL,
          plotPerOmic = FALSE,
          gene.col = "blue4")

## Given a regulator
plotGLM(GLMoutput = SimGLM,
          gene = NULL,
          regulator = "mmu-miR-15a-3p",
          plotPerOmic = FALSE,
          gene.col = "red4",
          regu.col = "green4")

## With cont.var parameter
plotGLM(GLMoutput = SimGLM,
          gene = "ENSMUSG00000026563",
          regulator = "Sp100",
          plotPerOmic = FALSE,
          gene.col = "darkblue",
          regu.col = "green4",
          cont.var = c(1,2,3,4,5,6,7,8),
          xlab = "CondA | CondB")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
