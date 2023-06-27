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
#'\code{GetGLM} fits a regression model for all the genes in the data set to identify
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
#' @param degree If cont.var is not NULL, non-linear (polynomial) relationships between the cont.var
#' and the gene expression must be studied. By default, degree = 1 and must be an integer.
#' A higher number will allow for in quadratic, cubic, etc. terms for cont.var in the model.
#' @param Res.df Number of degrees of freedom in the residuals. By default, 5. Increasing
#' Res.df will increase the power of the statistical model.
#' @param alfa Significance level. By default, 0.05.
#' @param MT.adjust Multiple testing correction method to be used within the stepwise variable
#' selection procedure. By default, "none". See the different options in ?\code{p.adjust}.
#' @param family Error distribution and link function to be used in the model (see ?\code{glm}
#' for more information). By default, negative.binomial(theta=10).
#' @param stepwise Stepwise variable selection method to be applied. It can be one of: "none", "backward"
#' (default), "forward", "two.ways.backward" or "two.ways.forward".
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
#'
#' @return
#' @export
#'
#' @examples
GetGLM = function(GeneExpression,
                  associations = NULL,
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
                  omic.type = 0){


  cont.var = NULL
  # Converting matrix to data.frame
  GeneExpression = as.data.frame(GeneExpression)
  data.omics = lapply(data.omics, as.data.frame)
  
  # If associations is NULL create a list of associations NULL
  if (is.null(associations)){
    associations=vector('list',length(data.omics))
    names(associations)=names(data.omics)
  }

  # Checking that Res.df is coherent with the number of samples
  if (Res.df >= (ncol(GeneExpression)-1)) stop("ERROR: You must decrease the Res.df so that a model can be computed.")


  # Checking that the number of samples per omic is equal to number of samples for gene expression and the number of samples for edesign
  for (i in 1:length(names(data.omics))){
    if(!length(colnames(data.omics[[i]])) == length(colnames(GeneExpression)) ) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }
  if(!is.null(edesign)){
    if(!length(colnames(GeneExpression)) == length(rownames(edesign)) ) {
      stop("ERROR: Samples in data.omics must be the same as in GeneExpression and in edesign")
    }
  }

  ## Checking that samples are in the same order in GeneExpressionDE, data.omics and edesign
  orderproblem<-FALSE
  if(is.null(edesign)){
    orderproblem<-all(sapply(data.omics, function(x) isTRUE(all.equal(sort(names(x)),sort(names(GeneExpression))))))
    if(orderproblem){
      data.omics<-lapply(data.omics, function(x) x[,order(colnames(GeneExpression))])
    }
    else{
      cat('Warning. GeneExpression and data.omics samples have not same names. We assume that they are ordered.\n')
    }
  } else{
    orderproblem<-all(c(sapply(data.omics, function(x) isTRUE(all.equal(sort(names(x)),sort(names(GeneExpression))))), isTRUE(all.equal(sort(rownames(edesign)),sort(names(GeneExpression))))))
    if(orderproblem){
      data.omics<-lapply(data.omics, function(x) x[,sort(colnames(GeneExpression))])
      edesign<-edesign[sort(colnames(GeneExpression)), , drop=FALSE]
    } else{
      cat('Warning. GeneExpression, edesign and data.omics samples have not same names. We assume that they are ordered.\n')
    }
  }

  ## Checking if there are regulators with "_R", "_P" or "_N" or with ":"
  message = FALSE
  for (i in 1:length(names(data.omics))){

    problemas = c(rownames(data.omics[[i]])[grep("_R$", rownames(data.omics[[i]]))],
                  rownames(data.omics[[i]])[grep("_P$", rownames(data.omics[[i]]))],
                  rownames(data.omics[[i]])[grep("_N$", rownames(data.omics[[i]]))])

    problema = c(grep(":", rownames(data.omics[[i]]), value = TRUE))
    rownames(data.omics[[i]]) = gsub(':', '-', rownames(data.omics[[i]]))
    rownames(data.omics[[i]]) = gsub('_R$', '-R', rownames(data.omics[[i]]))
    rownames(data.omics[[i]]) = gsub('_P$', '-P', rownames(data.omics[[i]]))
    rownames(data.omics[[i]]) = gsub('_N$', '-N', rownames(data.omics[[i]]))

    #Change the name in the association matrix only if associations is not NULL
    if(!is.null(associations[[i]])){
      associations[[i]][[2]]=gsub(':', '-', associations[[i]][[2]])
    }

    if(length(problemas) > 0) {
      cat("In",names(data.omics)[i], ',', problemas ,"regulators have names that may cause conflict with the algorithm by ending in _R, _P or _N", "\n")
      cat("Endings changed with -R, -P or -N, respectively", "\n")
    }

    if(length(problema) > 0) {
      cat("Some regulators in the omic", names(data.omics)[i],  "have names with \":\" that could cause conflict, replaced with \"-\" ", "\n")
      cat("Changed identifiers: ", problema, "\n")
    }
  }


  ##Checking that there are no replicates in the identifiers and changing identifiers in case of need
  
  if(length(names(data.omics))>1){
    for (i in 1:(length(names(data.omics))-1)){
      for(j in (i+1):(length(names(data.omics)))){
        repeated = intersect(rownames(data.omics[[i]]), rownames(data.omics[[j]]))
        
        
        if(length(repeated) > 0) {
          cat(names(data.omics)[i], "and", names(data.omics)[j], "omics have shared identifiers in regulators:", repeated, "\n")
          #Change the name in the association matrix only if is not NULL
          if(!is.null(associations[[i]])){
            associations[[i]][[2]][which(associations[[i]][[2]]==repeated)] =  paste(names(data.omics)[i],'-', repeated, sep='')
          }
          if(!is.null(associations[[j]])){
            associations[[j]][[2]][which(associations[[i]][[2]]==repeated)] =  paste(names(data.omics)[j],'-', repeated, sep='')
          }
          #Change the name in data.omics
          rownames(data.omics[[i]])[which(rownames(data.omics[[i]])==repeated)] =  paste(names(data.omics)[i],'-', repeated,sep='')
          rownames(data.omics[[j]])[which(rownames(data.omics[[j]])==repeated)] =  paste(names(data.omics)[j],'-', repeated,sep = '')
          
        }
      }
    }
  }

  

  # Preparing family for ElasticNet variable selection
  family2 = family$family
  family2 = strsplit(family2, "(", fixed = TRUE)[[1]][1]

  if (family2 %in% c("poisson", "quasipoisson", "Negative Binomial")) {
    family2 = "poisson"
  } else if (family2 %in% c("gaussian", "binomial")) {
    family2 = family2
  } else {
    family2 = NULL
    message(sprintf("Warning message:"))
    message(sprintf("Elasticnet variable selection cannot be applied for family %s", family2))
  }

  #Checking there are not -Inf/Inf values and eliminate genes/regulator that contain them
  infproblemgene<-is.infinite(rowSums(GeneExpression))
  infproblemreg<-lapply(data.omics, function(x) is.infinite(rowSums(x)))
  if(any(infproblemgene)){
    genesInf<-rownames(GeneExpression)[infproblemgene]
    GeneExpression<-GeneExpression[!infproblemgene,]
  }else{genesInf <-NULL}
  for (i in 1:length(names(data.omics))){
    if(any(infproblemreg[[i]])){
      cat(rownames(data.omics[[i]])[infproblemreg[[i]]], 'regulators of the omic', names(data.omics)[i] ,'have been deleted due to -Inf/Inf values. \n')
      data.omics[[i]]<-data.omics[[i]][!infproblemreg[[i]],]
    }
  }

  ## Removing genes with too many NAs and keeping track
  genesNotNA = apply(GeneExpression, 1, function (x) sum(!is.na(x)))
  genesNotNA = names(which(genesNotNA >= min.obs))
  genesNA = setdiff(rownames(GeneExpression), genesNotNA)
  GeneExpression = GeneExpression[genesNotNA,]

  ## Removing genes with no regulators only if associations does not have an associations = NULL in any omic
  genesNOreg = NULL
  genesNOreg = lapply(associations, function(x) if(!is.null(x)) {setdiff( rownames(GeneExpression),x[,1])})
  genesNOreg = Reduce(intersect, genesNOreg)
  GeneExpression = GeneExpression[!(rownames(GeneExpression) %in% genesNOreg),]
  if (length(genesNOreg) > 0){
    cat(length(genesNOreg), "genes had no initial regulators. Models will be computed for", length(rownames(GeneExpression)), 'genes.\n')
  }

  ## Removing constant genes
  constantGenes = apply(GeneExpression, 1, sd, na.rm = TRUE)
  notConstant = names(constantGenes)[constantGenes > 0]
  constantGenes = names(constantGenes)[constantGenes == 0]
  GeneExpression = GeneExpression[notConstant,]

  Allgenes=rownames(GeneExpression)
  nGenes = length(Allgenes)

  # Experimental groups
  if (is.null(edesign)) {
    cat("No experimental covariates were provided.\n")
    Group = 1:ncol(GeneExpression)
    names(Group) = colnames(GeneExpression)
    des.mat = NULL
  } else {
    Group = apply(edesign, 1, paste, collapse = "_")
    ## Experimental design matrix (with polynomial terms for cont.var, dummies for factors and interactions)
    # des.mat = GenerateDesignMatrix(interactions.exp, degree, edesign, cont.var)
    des.mat = model.matrix(~Group)[, -1, drop = FALSE]
    rownames(des.mat) = colnames(GeneExpression)
  }


  ## Omic types
  if (length(omic.type) == 1) omic.type = rep(omic.type, length(data.omics))
  names(omic.type) = names(data.omics)


  ## Remove regulators with NA
  cat("Removing regulators with missing values...\n")

  myregNA = lapply(data.omics, rownames)
  data.omics = lapply(data.omics, na.omit)
  myregNA = lapply(1:length(data.omics), function (i) setdiff(myregNA[[i]], rownames(data.omics[[i]])))
  names(myregNA)=names(data.omics)

  cat("Number of regulators with missing values:\n")
  print(sapply(myregNA, length))
  cat("\n")


  ## Remove regulators with Low Variability
  cat("Removing regulators with low variation...\n")

  tmp = LowVariationRegu(min.variation, data.omics, Group, associations, Allgenes, omic.type)
  data.omics = tmp[["data.omics"]]
  associations = tmp[["associations"]]
  myregLV = tmp[["myregLV"]]
  rm("tmp"); gc()


  ## Centering/Scaling quantitative predictors
  for (i in 1:length(omic.type)){
    if (omic.type[[i]] == 0){
      data.omics[[i]] = t(scale(t(data.omics[[i]]), center = center, scale = scale))
    }
  }


  ### Results objects

  ## Global summary for all genes
  GlobalSummary = vector("list", length = 4)
  names(GlobalSummary) = c("GoodnessOfFit", "ReguPerGene", "GenesNOmodel", "GenesNOregulators")

  GlobalSummary$GenesNOmodel = NULL
  if (length(genesNA) > 0) {
    GlobalSummary$GenesNOmodel = data.frame("gene" = genesNA,
                                            "problem" = rep("Too many missing values", length(genesNA)))
  }
  if (length(constantGenes) > 0) {
    GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                       data.frame("gene" = constantGenes,
                                                  "problem" = rep("Response values are constant", length(constantGenes))))
  }
  if (length(genesInf) > 0){
    GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                       data.frame("gene" = genesInf,
                                                  "problem" = rep("-Inf/Inf values", length(genesInf))))

  }

  GlobalSummary$GenesNOregulators = NULL
  if (length(genesNOreg) > 0){
    GlobalSummary$GenesNoregulators = data.frame("gene" = genesNOreg, "problem" = rep("Gene had no initial regulators", length(genesNOreg)))
  }

  GlobalSummary$GoodnessOfFit = matrix(NA, ncol = 4, nrow = nGenes)
  rownames(GlobalSummary$GoodnessOfFit) = Allgenes
  colnames(GlobalSummary$GoodnessOfFit) = c("Rsquared", "RMSE","CV(RMSE)", "relReg")

  GlobalSummary$ReguPerGene = matrix(0, ncol = 3*length(data.omics), nrow = nGenes)
  rownames(GlobalSummary$ReguPerGene) = Allgenes
  colnames(GlobalSummary$ReguPerGene) = c(paste(names(data.omics), "Ini", sep = "-"),
                                          paste(names(data.omics), "Mod", sep = "-"),
                                          paste(names(data.omics), "Rel", sep = "-"))

  ## Specific results for each gene
  ResultsPerGene=vector("list", length=length(Allgenes))
  names(ResultsPerGene) = Allgenes


  ### Computing model for each gene
  cat("Checking multicollinearity, selecting predictors and fitting model for ...\n")

  pap = c(1, 1:round(nGenes/100) * 100, nGenes)

  for (i in 1:nGenes) {

    gene=Allgenes[i]

    ResultsPerGene[[i]] = vector("list", length = 5)
    names(ResultsPerGene[[i]]) = c("Y", "X", "coefficients", "allRegulators", "relevantRegulators")

    if (is.element(i, pap)) cat(paste("Fitting model for gene", i, "out of", nGenes, "\n"))

    RetRegul = GetAllReg(gene=gene, associations=associations, data.omics = data.omics)
    RetRegul.gene = RetRegul$Results  ## RetRegul$TableGene: nr reg per omic
    ## Some of these reg will be removed, because they are not in data.omics


    # RetRegul.gene--> gene/regulator/omic/area
    RetRegul.gene=RetRegul.gene[RetRegul.gene[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators


    ### NO INITIAL REGULATORS
    if(length(RetRegul.gene)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experiment

      if (is.null(edesign)) {
        ResultsPerGene[[i]]$X = NULL
        ResultsPerGene[[i]]$relevantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
        isModel = NULL

      } else {
        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
        colnames(des.mat2)[1] = "response"
        des.mat2 = na.omit(des.mat2)

        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        des.mat2 = des.mat2[,sdNo0]

        myGLM = ComputeGLM(matrix.temp = data.frame(des.mat2, check.names = FALSE),
                           alfa = alfa, stepwise = 'none', Res.df = Res.df,
                           family = family, epsilon = epsilon, MT.adjust = MT.adjust)

        ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
        ResultsPerGene[[i]]$relevantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
      }


      # GlobalSummary$ReguPerGene  # this is initially set to 0 so no need to modify it


      ### WITH INITIAL REGULATORS
    }
    else { ## There are regulators for this gene at the beginning

      ResultsPerGene[[i]]$allRegulators = data.frame(RetRegul.gene, rep("Model",nrow(RetRegul.gene)), stringsAsFactors = FALSE)
      colnames(ResultsPerGene[[i]]$allRegulators) = c("gene","regulator","omic","area","filter")

      GlobalSummary$ReguPerGene[gene, grep("-Ini", colnames(GlobalSummary$ReguPerGene))] = as.numeric(RetRegul$TableGene[-1])
      # the rest of columns remain 0

      ## Identify which regulators where removed because of missing values or low variation
      res = RemovedRegulators(RetRegul.gene = ResultsPerGene[[i]]$allRegulators,
                              myregLV=myregLV, myregNA=myregNA, data.omics=data.omics)

      if(length(res$RegulatorMatrix)==0){ ## No regulators left after the filtering to compute the model

        if (is.null(edesign)) {
          ResultsPerGene[[i]]$X = NULL
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene
          myGLM = NULL

        } else {
          des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
          colnames(des.mat2)[1] = "response"
          des.mat2 = na.omit(des.mat2)

          # Removing predictors with constant values
          sdNo0 = apply(des.mat2, 2, sd)
          sdNo0 = names(sdNo0)[sdNo0 > 0]
          des.mat2 = des.mat2[,sdNo0]

          myGLM = ComputeGLM(matrix.temp = data.frame(des.mat2, check.names = FALSE),
                             alfa = alfa, stepwise = 'none', Res.df = Res.df,
                             family = family, epsilon = epsilon, MT.adjust = MT.adjust)

          ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene

        }

      }
      else {  ## Regulators for the model!!

        ## Apply multicollinearity filter only if there is more than one regulator for a gene
        
        if (ncol(res$RegulatorMatrix)>1){
          
          res = CollinearityFilter(data = res$RegulatorMatrix, reg.table = res$SummaryPerGene,
                                 correlation = correlation, omic.type = omic.type)
        }

        ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene


        ## Creating data matrix with regulators and with/without interactions

        ## Scaling predictors for ElasticNet only in case they were not already scaled
        if (!is.null(elasticnet)) {
          if (scale) {
            des.mat2EN = RegulatorsInteractions(interactions.reg, reguValues = res$RegulatorMatrix,
                                                des.mat, cont.var, GeneExpression, gene)
          } else {
            ScaleMatrix = res$RegulatorMatrix
            for(k in 1:ncol(ScaleMatrix)){
              if(any(res$RegulatorMatrix[,k] != 1 & res$RegulatorMatrix[,k] != 0)){
                ScaleMatrix[,k] = scale(ScaleMatrix[,k])
              }
            }

            des.mat2EN = RegulatorsInteractions(interactions.reg,
                                                reguValues = ScaleMatrix,
                                                des.mat, cont.var, GeneExpression, gene)
          }
        } else {
          des.mat2EN = RegulatorsInteractions(interactions.reg, reguValues = res$RegulatorMatrix,
                                              des.mat, cont.var, GeneExpression, gene)
        }


        # Removing observations with missing values
        des.mat2EN = na.omit(des.mat2EN)

        ###  Variable selection --> Elasticnet
        tmp = ElasticNet(family2, des.mat2EN, epsilon, elasticnet, Res.df)
        regulatorcoef = tmp[['coefficients']]
        isModel = tmp[['isModel']]
        m = tmp[['m']]
        des.mat2 = as.data.frame(des.mat2EN[,colnames(tmp[["des.mat2"]]),drop = FALSE])
        ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
        rm(des.mat2EN); gc()
        
        if (ncol(des.mat2) == 1 || is.null(isModel)) {
          
          ## Extracting significant regulators
          ResultsPerGene[[i]]$relevantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
          
          ## Counting original regulators in the model per omic
          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
        } else{
          
          isModel = TRUE
          mycoef = colnames(des.mat2[,-1,drop = FALSE])
          myvariables = unlist(strsplit(mycoef, ":", fixed = TRUE))
          mycondi = intersect(myvariables, colnames(des.mat))
          myvariables = intersect(myvariables, rownames(ResultsPerGene[[i]]$allRegulators))
          
          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Rel" = 0, stringsAsFactors = FALSE)
          ResultsPerGene[[i]]$allRegulators[myvariables, "Rel"] = 1
          ResultsPerGene[[i]]$coefficients = regulatorcoef
          #ResultsPerGene[[i]]$coefficients = regulatorcoef[myvariables,, drop =FALSE]
          colnames(ResultsPerGene[[i]]$coefficients) = c('coefficient')
          
          ## A las variables significativas le quito "_R", solo quedara omica_mc"num". Luego creo un objeto que contenga a mi tabla de "allRegulators"
          ## para poder modificar los nombres de la misma forma.
          myvariables = sub("_R", "", myvariables)
          ResultsPerGene[[i]]$relevantRegulators = myvariables # significant regulators including "new" correlated regulators without _R
          
          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
          
          mytable = ResultsPerGene[[i]]$allRegulators
          mytable[,"filter"] = sub("_P", "", mytable[,"filter"])
          mytable[,"filter"] = sub("_N", "", mytable[,"filter"])
          mytable[,"filter"] = sub("_R", "", mytable[,"filter"])
          
          collin.regulators = intersect(myvariables, mytable[,"filter"])
          
          if (length(collin.regulators) > 0) {  # there were correlated regulators
            original.regulators = mytable[mytable[,"filter"] %in% collin.regulators, "regulator"]
            
            ResultsPerGene[[i]]$allRegulators[original.regulators, "Rel"] = 1
            ResultsPerGene[[i]]$relevantRegulators = c(ResultsPerGene[[i]]$relevantRegulators, as.character(original.regulators))
            ResultsPerGene[[i]]$relevantRegulators = setdiff(ResultsPerGene[[i]]$relevantRegulators, collin.regulators)
            
            ## Counting original regulators in the model per omic
            contando = ResultsPerGene[[i]]$allRegulators
            quitar = which(contando[,"filter"] == "MissingValue")
            if (length(quitar) > 0) contando = contando[-quitar,]
            quitar = which(contando[,"filter"] == "LowVariation")
            if (length(quitar) > 0) contando = contando[-quitar,]
            contando = contando[-grep("_R", rownames(contando)),]
          } else {
            contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          }
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
          
          ## TO DO: Corregir GblobalSummary$ReguPerGene
          ## Counting significant regulators per omic
          if (length(ResultsPerGene[[i]]$relevantRegulators) > 0) {
            contando = ResultsPerGene[[i]]$allRegulators[ResultsPerGene[[i]]$relevantRegulators,, drop=FALSE]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(data.omics)])
            contando[is.na(contando)] = 0
            GlobalSummary$ReguPerGene[gene, grep("-Rel", colnames(GlobalSummary$ReguPerGene))] = contando
          }
          
        }
      }

      } ## Close "else" --> None regulators from begining

    if (is.null(isModel)) {
      
      ResultsPerGene[[i]]$Y = GeneExpression[i,]
      ResultsPerGene[[i]]$coefficients = NULL
      
      GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != gene,]
      
      
    } else {
      ResultsPerGene[[i]]$Y = data.frame("y" = des.mat2[,1], "fitted.y" = tmp[['fitted.values']],
                                         "residuals" = des.mat2[,1] - tmp[['fitted.values']], check.names = FALSE)
      colnames(ResultsPerGene[[i]]$Y) <- c("y", "fitted.y", "residuals")
      GlobalSummary$GoodnessOfFit[gene,] = c(m$R.squared, m$RMSE, m$cvRMSE,as.character(length(ResultsPerGene[[gene]]$relevantRegulators)))
      
      
    }  

  }  ## At this point the loop for all genes is finished
  
  # Remove from GoodnessOfFit genes with no significant regulators
  
  genesNosig = names(which(GlobalSummary$GoodnessOfFit[,4]==0))
  genessig = setdiff(rownames(GlobalSummary$GoodnessOfFit), genesNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[genessig,]

  myarguments = list(edesign = edesign, finaldesign = des.mat, Res.df = Res.df, groups = Group, alfa = alfa, family = family,
                     stepwise = 'none', center = center, scale = scale, elasticnet = elasticnet,
                     min.variation = min.variation, correlation = correlation,
                     MT.adjust = MT.adjust, min.obs = min.obs, epsilon = epsilon,
                     GeneExpression = GeneExpression, dataOmics = data.omics, omic.type = omic.type)

  return(list("ResultsPerGene" = ResultsPerGene, "GlobalSummary" = GlobalSummary, "arguments" = myarguments))

}






# Get all regulators ------------------------------------------------------


## Auxiliar function to paste different areas

uniquePaste = function (x) {
  x = unique(x)

  if (length(x) > 1) {
    x = paste(x, collapse = ";")
  }

  return (x)
}




## INPUT
## gene: we are seeking information
## associations: List containing as many elements as original association files between genes and each omic (association$miRNA, association$DNase, etc.)
## data.omics: For removing regulators with no-information

GetAllReg=function(gene, associations, data.omics){

  Reg.matrix=NULL
  NrReg=NULL
  myomic=names(data.omics)
  
  for(ov in myomic){
    if (is.null(associations[[ov]])){
      myregulators=rownames(data.omics[[ov]])
      ov.nr=length(myregulators)
      
      if(ov.nr==0){
        myregulators=c("No-regulator")
        myregulators=t(as.matrix(myregulators))
      }
      
      Reg.matrix.temp=cbind(rep(gene,ov.nr),myregulators,ov,"")  ## Lo tengo que dejar igual que las otras omicas
      NrReg.temp=ov.nr
      colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
      Reg.matrix.temp = t(apply(Reg.matrix.temp,1, as.character))  ## adding this here to avoid factors
      Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
      colnames(Reg.matrix)=c("gene","regulator","omic","area")
      NrReg=c(NrReg, NrReg.temp)
      
    } else{

    colnames(associations[[ov]])[1]="gene"

    ## Regulator with area

    if(ncol(associations[[ov]])>2){

      myregulators=associations[[ov]][associations[[ov]]$gene==gene, ,drop=FALSE] ## "regulator" with Area--> Matrix

      if(length(myregulators)==0){
        ov.nr=length(myregulators) ## regulators nr per omic --> It could be 0
        myregulators=c("No-regulator","")
        myregulators=t(as.matrix(myregulators))

      } else {
        myregulators=aggregate(myregulators,by=list(myregulators[,2]),uniquePaste) ## Agrupado por "region", me devuelve Area separado por comas
        myregulators=myregulators[,-c(1:2)]
        myregulators[,2]=sapply(myregulators[,2], paste, collapse = ";")
        ov.nr=nrow(myregulators)
      }

      Reg.matrix.temp = cbind(gene,myregulators,ov)
      Reg.matrix.temp = Reg.matrix.temp[,c(1,2,4,3),drop=FALSE] ## Dejo el mismo orden que tenia
      colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
      NrReg.temp=ov.nr
      Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
      colnames(Reg.matrix)=c("gene","regulator","omic","area")
      NrReg=c(NrReg, NrReg.temp)

      ## Regulators without area

    } else {

      myregulators=associations[[ov]][associations[[ov]]$gene==gene,2]
      myregulators=unique(myregulators) ## Could it be repeated??
      ov.nr=length(myregulators) ## regulators nr per omic --> It could be 0

      if(ov.nr==0){
        myregulators=c("No-regulator")
        myregulators=t(as.matrix(myregulators))
      }

      Reg.matrix.temp=cbind(gene,myregulators,ov,"")  ## Lo tengo que dejar igual que las otras omicas
      NrReg.temp=ov.nr
      colnames(Reg.matrix.temp)=c("gene","regulator","omic","area")
      Reg.matrix.temp = t(apply(Reg.matrix.temp,1, as.character))  ## adding this here to avoid factors
      Reg.matrix=rbind(Reg.matrix,Reg.matrix.temp)
      colnames(Reg.matrix)=c("gene","regulator","omic","area")
      NrReg=c(NrReg, NrReg.temp)
      }
    }
  }
  Results2=c(gene,NrReg)
  myomic=paste(myomic,"Ini",sep="-")
  names(Results2)=c("gene",myomic)

  Results=vector("list", length=2)
  Results[[1]]=Reg.matrix
  Results[[2]]=Results2
  names(Results)=c("Results","TableGene")

  return(Results)

}






# Obtain which regulators have been removed and why -----------------------
## Input:
## RetRegul.gene: Initial regulator matrix
## myregLV: regulators removed for low variability. List by omics
## myregNA: regulators removed for NA. List by omics

RemovedRegulators = function(RetRegul.gene, myregLV, myregNA, data.omics){

  RegulatorsValue=NULL
  RetRegul.geneNEW = NULL

  rownames(RetRegul.gene) = RetRegul.gene[,"regulator"]
  myregini = RetRegul.gene[,"regulator"]  ## In our case, one regulator cannot belong to 2 omics
  mygene = RetRegul.gene[1,"gene"]

  for(ov in unique(RetRegul.gene[,"omic"])){

    ## remove regulators not in data.omics
    regmodel = intersect(RetRegul.gene[,"regulator"], rownames(data.omics[[ov]]))
    RegulatorsValue = cbind(RegulatorsValue, t(data.omics[[ov]][regmodel,,drop=FALSE]))
    RetRegul.geneNEW = rbind(RetRegul.geneNEW, RetRegul.gene[regmodel, , drop=FALSE])

    ## NA

    if(length(intersect(myregini, myregNA[[ov]]))>0){
      RetRegul.geneNEW = rbind(RetRegul.geneNEW, data.frame(gene = mygene, regulator = intersect(myregini, myregNA[[ov]]),
                                                            omic = ov, area=RetRegul.gene[intersect(myregini, myregNA[[ov]]),"area"],
                                                            filter = "MissingValue", stringsAsFactors = FALSE))
    }

    ## LV

    if (length(intersect(myregini, myregLV[[ov]]))>0){
      RetRegul.geneNEW = rbind(RetRegul.geneNEW, data.frame(gene = mygene, regulator = intersect(myregini, myregLV[[ov]]),
                                                            omic = ov, area=RetRegul.gene[intersect(myregini, myregLV[[ov]]),"area"],
                                                            filter = "LowVariation", stringsAsFactors = FALSE))
    }


    ## RetRegul.geneNEW=as.data.frame(RetRegul.geneNEW) ## Si es una matriz con una fila,al hacer el apply me lo convierte en vector
    RetRegul.geneNEW=apply(RetRegul.geneNEW,c(1,2),as.character)
    ##RetRegul.geneNEW[,"filter"] = as.character(RetRegul.geneNEW[,"filter"])


  }

  return(list("SummaryPerGene" = RetRegul.geneNEW, "RegulatorMatrix" = RegulatorsValue))

}




# Checking multi-collinearity ---------------------------------------------

CollinearityFilter = function(data, reg.table, correlation = 0.8, omic.type) {

  ## data = Regulator data matrix for all omics where missing values and regulators with low variation have been filtered out
  #         (regulators must be in columns)
  ## reg.table = Table with "gene", "regulator", "omic", "area", filter" where omics with no regulators have been removed

  row.names(reg.table) = reg.table[,"regulator"]
  #resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)

  for (j in 1:length(omic.type)){

    if(omic.type[[j]] == 0){  # continuous variable

      # Initial regulators
      myreg = reg.table[which(reg.table[,"omic"] == names(omic.type)[[j]]), ,drop=FALSE]
      myreg = as.character(myreg[which(myreg[,"filter"] == "Model"),"regulator"])

      if (length(myreg) > 1) {  # if there is more than one regulator for this omic:

        # Dejo la matriz original para poder tomar luego las correlaciones de la red.
        mycorrelations = data.frame(t(combn(myreg,2)), as.numeric(as.dist(cor(data[,myreg]))), stringsAsFactors = FALSE)
        mycor = mycorrelations[abs(mycorrelations[,3]) >= correlation,]

        if (nrow(mycor) == 1) {  ### only 2 regulators are correlated in this omic

          correlacionados = unlist(mycor[,1:2])
          regulators = colnames(data)
          keep = sample(correlacionados, 1) # Regulador al azar de la pareja

          ## Lo siguiente elimina el no representante de la matriz de reguladores. Al regulador escogido como representante,
          ## le cambia el nombre por "mc_1_R" para que despues pase la seleccion de variables y asi, en reg.table se conserva
          ## la info de que fue escogido como representante.

          remove = setdiff(correlacionados, keep)
          regulators = setdiff(regulators, remove)
          data = as.matrix(data[ ,regulators])
          colnames(data) = regulators
          index.reg = which(colnames(data) == as.character(keep))
          colnames(data)[index.reg] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "R", sep = "_")

          # Cambio en reg.table. Asignacion de los nombres segun sea representante,
          # correlacion positiva o negativa. Creacion de una nueva fila con el representante
          # para la seleccion de variables y asi, no perder la info del representante.

          reg.table = rbind(reg.table, reg.table[keep,])
          reg.table[nrow(reg.table), "regulator"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "R", sep = "_")
          reg.table[keep, "filter"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "R", sep = "_")
          rownames(reg.table) = reg.table[ ,"regulator"]

          if(mycor[,3] > 0){
            index = mycor[1, which(with(mycor, mycor[1 ,c(1,2)] != keep))]
            reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "P", sep = "_")
          } else{
            index = mycor[1, which(with(mycor, mycor[1, c(1,2)] != keep))]
            reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "N", sep = "_")
          }
        }

        if (nrow(mycor) >= 2) {   ### more than 2 regulators might be correlated in this omic

          mygraph = graph.data.frame(mycor, directed=F)
          mycomponents = clusters(mygraph)

          for (i in 1:mycomponents$no) {
            correlacionados = names(mycomponents$membership[mycomponents$membership == i])
            regulators = colnames(data)

            ## Escoge un regulador al azar como representante de cada componente conexa. Para cada componente conexa elimina aquellos reguladores
            ## que no han sido escogidos como representante.
            keep = sample(correlacionados, 1)  # mantiene uno al azar
            reg.remove = setdiff(correlacionados, keep) # correlated regulator to remove
            regulators = setdiff(regulators, reg.remove)  # all regulators to keep
            data = as.matrix(data[ ,regulators])
            colnames(data) = regulators
            index.reg = which(colnames(data) == as.character(keep))
            colnames(data)[index.reg] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "R", sep = "_")

            # Asignacion de nombre al representante y nueva fila para el filtro de
            # seleccion de variables (asi no se pierde la info del representante).
            reg.table = rbind(reg.table, reg.table[keep,])
            reg.table[nrow(reg.table), "regulator"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "R", sep = "_")
            reg.table[keep, "filter"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "R", sep = "_")
            rownames(reg.table) = reg.table[ ,"regulator"]

            # Matriz que recoge los reguladores correlacionados con el representante: actual.correlation. Asi se puede ver si la correlacion es
            # positiva o negativa y asignar el nombre. Intente hacer merge(), expand.grid(), pero no daba las mismas combinaciones que combn(), por lo que
            # vi necesario hacer un bucle para quedarme con aquellas parejas que interesan (representante - resto de reguladores).
            actual.couple = data.frame(t(combn(correlacionados,2)), stringsAsFactors = FALSE)
            colnames(actual.couple) = colnames(mycorrelations[,c(1,2)])

            actual.correlation = NULL
            for(k in 1:nrow(actual.couple)){
              if (any(actual.couple[k,c(1,2)] == keep)){
                actual.correlation = rbind(actual.correlation, actual.couple[k,])
              }
            }

            actual.correlation = merge(actual.correlation[,c(1,2)],mycorrelations)

            # Uso la matriz anterior para recorrer las correlaciones y segun sea positiva
            # o negativa, asigno un nombre.
            for(k in 1:nrow(actual.correlation)){
              if(actual.correlation[k,3] > 0){
                index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
                reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "P", sep = "_")
              } else{
                index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
                reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "N", sep = "_")
              }
            }
          }
        }
      }
    } else {

      # Variables categoricas. El procedimiento es analogo, solo que las correlaciones
      # se obtienen mediante tablas de contingencia 2 a 2. Se usa el coeficiente de correlacion phi.
      myreg = reg.table[which(reg.table[,"omic"] == names(omic.type)[[j]]), ,drop=FALSE]
      myreg = as.character(myreg[which(myreg[,"filter"] == "Model"),"regulator"])

      if (length(myreg) > 1){

        # Phi. Contingency table.
        couple = t(combn(myreg,2))
        categorical.correlation = NULL

        # Aqui hago las tablas de contingencia y almaceno en la matriz categorical.correlation las correlaciones por
        # parejas de reguladores.
        for (k in 1:nrow(couple)){
          contingency.table = table(data[,couple[k,1]], data[,couple[k,2]])
          categorical.correlation = rbind(categorical.correlation, phi(contingency.table))
        }

        # He creido conveniente tomar una correlacion de 0.6 por defecto para las variables
        # cualitativas, ya que con ejemplos he visto que 0.6 ya presenta un alto numero de
        # valores iguales en mismas posiciones.
        mycorrelations = data.frame(couple, categorical.correlation)
        correlations = mycorrelations[abs(mycorrelations[,3]) >= correlation,]


        if (nrow(correlations) == 1) { ### only 2 regulators are correlated in this omic

          correlacionados = as.character(unlist(correlations[,1:2]))
          regulators = colnames(data)
          keep = sample(correlacionados, 1)
          remove = setdiff(correlacionados, keep)
          regulators = setdiff(regulators, remove)  # all regulators to keep
          data = as.matrix(data[ ,regulators])
          colnames(data) = regulators
          index.reg = which(colnames(data) == as.character(keep))
          colnames(data)[index.reg] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "R", sep = "_")

          reg.table = rbind(reg.table, reg.table[keep,])
          reg.table[nrow(reg.table), "regulator"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "R", sep = "_")
          reg.table[keep, "filter"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "R", sep = "_")
          rownames(reg.table) = reg.table[ ,"regulator"]

          if(correlations[,3] > 0){
            index = as.character(correlations[1, which(with(correlations, correlations[1,c(1,2)] != keep))])
            reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "P", sep = "_")
          } else{
            index = as.character(correlations[1, which(with(correlations, correlations[1,c(1,2)] != keep))])
            reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", 1, sep = ""), "N", sep = "_")
          }

        }

        if (nrow(correlations) >= 2) {   ### more than 2 regulators might be correlated in this omic

          mygraph = graph.data.frame(correlations, directed=F)
          mycomponents = clusters(mygraph)

          for (i in 1:mycomponents$no) {
            correlacionados = names(mycomponents$membership[mycomponents$membership == i])
            regulators = colnames(data)

            keep = sample(correlacionados, 1)  # mantiene uno al azar
            reg.remove = setdiff(correlacionados, keep) # correlated regulator to remove
            regulators = setdiff(regulators, reg.remove)  # all regulators to keep
            data = as.matrix(data[ ,regulators])
            colnames(data) = regulators
            index.reg = which(colnames(data) == as.character(keep))
            colnames(data)[index.reg] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "R", sep = "_")

            reg.table = rbind(reg.table, reg.table[keep,])
            reg.table[nrow(reg.table), "regulator"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "R", sep = "_")
            reg.table[keep, "filter"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "R", sep = "_")
            rownames(reg.table) = reg.table[ ,"regulator"]

            actual.couple = as.data.frame(t(combn(correlacionados,2)))
            colnames(actual.couple) = colnames(mycorrelations[,c(1,2)])

            actual.correlation = NULL
            for(k in 1:nrow(actual.couple)){
              if (any(actual.couple[k,c(1,2)] == keep)){
                actual.correlation = rbind(actual.correlation, actual.couple[k,])
              }
            }

            actual.correlation = merge(actual.correlation[,c(1,2)],mycorrelations)

            for(k in 1:nrow(actual.correlation)){
              if(actual.correlation[k,3] > 0){
                index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
                reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "P", sep = "_")
              } else{
                index = as.character(actual.correlation[k, which(with(actual.correlation, actual.correlation[k,c(1,2)] != keep))])
                reg.table[index, "filter"] = paste(names(omic.type)[[j]], paste("mc", i, sep = ""), "N", sep = "_")
              }
            }
          }
        }
      }
    }
  }
  resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)
  rownames(resultado$SummaryPerGene) = resultado$SummaryPerGene[,"regulator"]
  return(resultado)
}









# CollinearityFilterOLD = function(data, reg.table, correlation = 0.8, omic.type) {
#
#   ## data = Regulator data matrix for all omics where missing values and regulators with low variation have been filtered out
#   #         (regulators must be in columns)
#   ## reg.table = Table with "gene", "regulator", "omic", "area", filter" where omics with no regulators have been removed
#
#   row.names(reg.table) = reg.table[,"regulator"]
#   # resultado = list(RegulatorMatrix = data, SummaryPerGene = reg.table)
#
#   for (omic in unique(reg.table[,"omic"])) {
#
#     # Initial regulators
#     myreg = reg.table[which(reg.table[,"omic"] == omic), ,drop=FALSE]
#     myreg = as.character(myreg[which(myreg[,"filter"] == "Model"),"regulator"])
#
#     if (length(myreg) > 1) {  # if there is more than one regulator for this omic:
#
#       # Checking if correlation is higher than threshold value
#       mycor = data.frame(t(combn(myreg,2)), as.numeric(as.dist(cor(data[,myreg]))), stringsAsFactors = FALSE)
#       mycor = mycor[mycor[,3] >= correlation,]
#
#       if (nrow(mycor) == 1) {  ### only 2 regulators are correlated in this omic
#         correlacionados = unlist(mycor[,1:2])
#         resultado = CorreAction(correlacionados, data = resultado$RegulatorMatrix, reg.table = resultado$SummaryPerGene,
#                                 action = action, nombre = "mc1", omic = omic)
#       }
#
#       if (nrow(mycor) >= 2) {   ### more than 2 regulators might be correlated in this omic
#
#         mygraph = graph.data.frame(mycor, directed=F)
#         mycomponents = clusters(mygraph)
#
#         for (i in 1:mycomponents$no) {
#           correlacionados = names(mycomponents$membership[mycomponents$membership == i])
#           resultado = CorreAction(correlacionados, data = resultado$RegulatorMatrix, reg.table = resultado$SummaryPerGene,
#                                   action = action, nombre = paste("mc", i, sep = ""), omic = omic)
#         }
#       }
#
#     }
#
#   }
#
#   rownames(resultado$SummaryPerGene) = resultado$SummaryPerGene[,"regulator"]
#   return(resultado)
# }




# Action to perform when having correlated regulators ---------------------

# # CorreAction = function (correlacionados, data, reg.table, omic, action = c("mean", "random"), nombre = "mc1") {
#
#   regus = colnames(data)
#
#   if (action == "random") {
#     dejar = sample(correlacionados, 1)  # correlated regulator to keep
#     quitar = setdiff(correlacionados,dejar) # correlated regulator to remove
#     regus = setdiff(regus, quitar)  # all regulators to keep
#     data = data[,regus]  # excluding filtered regulator from data
#     reg.table[quitar,"filter"] = dejar    # update info in regulators table
#     reg.table = apply(reg.table, 2, as.character)
#     rownames(reg.table) = reg.table[,"regulator"]
#   }
#
#   if (action == "mean") {
#     nuevo = rowMeans(data[,correlacionados])  # mean of the correlated regulators
#
#     # Collapse the areas of all the regulators
#     myarea = unique(reg.table[correlacionados, "area"])
#     if (length(grep(";", myarea)) > 0) {
#       myarea = unique(unlist(sapply(myarea, strsplit, ";")))
#     }
#     myarea = paste(myarea, collapse = ";")
#
#     regus = setdiff(regus, correlacionados) # remove correlated regulators
#     data = data[,regus]   # new data without correlated regulators
#     data = cbind(data, nuevo)  # add new "mean" regulator
#     nombre = paste(omic, nombre, sep = "_")
#     colnames(data)[ncol(data)] = nombre
#     reg.table[correlacionados,"filter"] = rep(nombre, length(correlacionados)) # update info in regulators table
#     reg.table = rbind(reg.table,
#                       data.frame("gene" = reg.table[1,1], "regulator" = nombre, "omic" = omic, "area" = myarea, "filter" = "Model",
#                                  stringsAsFactors = FALSE))  # add a new row to regulators table with new regulators
#     reg.table = apply(reg.table, 2, as.character)
#     rownames(reg.table) = reg.table[,"regulator"]
#   }
#
#   return(list(RegulatorMatrix = data, SummaryPerGene = reg.table))
#
# }












# Plot GLM results --------------------------------------------------------

# order: Should the experimental groups be ordered for the plot? If TRUE, omic values are also ordered accordingly.
#        If FALSE, the function assumes they were provided in the right order for a meaningful plot.

plotGLM = function (GLMoutput, gene, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL,...) {

  # Colors for omics
  omic.col = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390,
                        100,200,300,400,500,10,20,30,40,50,60,70,80,90,150,250,350,450,550)]

  if (is.null(regu.col)) {
    any.col = omic.col[1:length(GLMoutput$arguments$dataOmics)]
  } else {
    if (length(regu.col) == length(GLMoutput$arguments$dataOmics)) {
      any.col = regu.col
    } else {
      any.col = rep(regu.col, length(GLMoutput$arguments$dataOmics))
    }
  }
  names(any.col) = names(GLMoutput$arguments$dataOmics)


  # Changing margin
  par(mar = c(5,3,4,3)+0.1)

  # Groups to plot
  if (is.null(cond2plot)) {
    if (!is.null(GLMoutput$arguments$edesign)) {
      cond2plot = apply(GLMoutput$arguments$edesign, 1, paste, collapse = "_")
    }
  }

  # Replicates
  if (!is.null(cont.var)) {  # we have continuous variable
    if (!is.null(cond2plot)) { # cont.var + cond2plot
      myreplicates = apply(cbind(cond2plot, cont.var), 1, paste, collapse = "_")

    } else {  # only continuous variable is provided
      myreplicates = cont.var
    }

  } else {   # no cont.var

    if (!is.null(cond2plot)) { # only cond2plot
      myreplicates = colnames(GLMoutput$arguments$GeneExpression)
    } else {  # nothing
      myreplicates = colnames(GLMoutput$arguments$GeneExpression)
    }
  }

  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = colnames(GLMoutput$arguments$GeneExpression)
  if (order) {
    myorder = order(myreplicates)
    myreplicates = sort(myreplicates)
    cond2plot = cond2plot[myorder]
  }
  myrepliUni = unique(myreplicates)


  if (max(table(myreplicates)) == 1) {
    replicates = FALSE
  } else {
    replicates = TRUE
  }


  if (is.null(cond2plot)) {
    numLines = NULL
  } else {
    condi1 = unique(cond2plot)
    num1 = 1:length(condi1); names(num1) = condi1
    num2 = num1[as.character(cond2plot)]
    if (replicates) {
      num2 = aggregate(num2, by = list("rep" = myreplicates), unique)$x
    }
    numLines = which(diff(num2) != 0)+0.5
  }


  # Error values
  getErrorValues = function(realValues, repsInfo) {
    # Disable it
    if (! replicates)
      return(NULL)

    out_values = tapply(realValues, repsInfo, function(reps) sd(reps)/sqrt(length(reps)))

    return(out_values)
  }


  myrepliUni = unique(myreplicates)

  ## REGULATOR = NULL

  if (is.null(regulator)) {  ### Plot all regulators for the given gene

    GLMgene = GLMoutput$ResultsPerGene[[gene]]

    if (is.null(GLMgene)) {
      stop(paste("No GLM was obtained for gene", gene))
    }

    if (is.null(GLMgene$significantRegulators)) { ## No significant regulators

      cat("No significant regulators were found for this gene.\n")


    } else
    {  ## Significant regulators:

      # Considering multicollinearity
      SigReg = GLMgene$allRegulators
      SigReg = SigReg[SigReg$Sig == 1, c("regulator", "omic", "area", "filter")]

      SigReg = SigReg[GLMgene$significantRegulators,,drop = FALSE]

      cat(paste(nrow(SigReg), "significant regulators are to be plotted for gene", gene)); cat("\n")

      # Gene values
      geneValues = GLMgene$Y$y
      if (order) geneValues = geneValues[myorder]
      errorValues = getErrorValues(geneValues, myreplicates)
      geneValues = tapply(geneValues, myreplicates, mean)
      geneValues = geneValues[myrepliUni]
      names(geneValues) = myrepliUni

      # X values
      x.points = 1:length(myrepliUni)
      eje = myrepliUni

      if (plotPerOmic) { ## All regulators from the same omic in the same plot

        myomics = unique(SigReg$omic)

        for (oo in myomics) {

          SigRegOmic = SigReg[SigReg$omic == oo,]

          omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
          omicValues = omicValues[rownames(GLMgene$X),]
          reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[1]]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)

          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni

          mycol = any.col[oo]

          if (nrow(SigRegOmic) == 1) {
            leftlab = SigRegOmic$regulator[1]
          } else { leftlab = oo }

          yleftlim = range(omicValues[,SigRegOmic$regulator], na.rm = TRUE)

          if (! is.null(errorValuesRegu)) {
            yleftlim = range(
              apply(omicValues[, SigRegOmic$regulator, drop = FALSE], 2, function(x) {
                if (order) x = x[myorder]
                errorInd = getErrorValues(x, myreplicates)
                meanValues = tapply(x, myreplicates, mean)
                meanValues = meanValues[myrepliUni]

                return(c(meanValues - errorInd, meanValues + errorInd))
              }))
          }

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), yleftlim = yleftlim,
                       xlab = xlab, yylab = c(gene, leftlab), pch = c(16,16),
                       main = oo, numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)

          if (nrow(SigRegOmic) > 1) {
            for (i in 2:nrow(SigRegOmic)) {

              reguValues = omicValues[, colnames(omicValues) == SigRegOmic$regulator[i]]
              if (order) reguValues = reguValues[myorder]
              errorValuesRegu = getErrorValues(reguValues, myreplicates)

              reguValues = tapply(reguValues, myreplicates, mean)
              reguValues = reguValues[myrepliUni]
              names(reguValues) = myrepliUni

              lines(x.points, reguValues, type = "o", lwd = 2, pch = i, col = mycol, lty = i)

              if (! is.null(errorValuesRegu)) {
                arrows(x.points, reguValues - errorValuesRegu, x.points, reguValues + errorValuesRegu,
                       code = 3, length = 0.02, angle = 90, col = mycol)
              }
            }
          }

        }

      } else {  ## Each regulator in a separate plot

        for (rr in SigReg$regulator) {

          oo = GLMoutput$ResultsPerGene[[gene]]$allRegulators[rr,"omic"]
          omicValues = t(GLMoutput$arguments$dataOmics[[oo]])
          reguValues = omicValues[, colnames(omicValues) == rr]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)

          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni

          mycol = any.col[oo]

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, mycol), xlab = xlab,
                       yylab = c(gene, rr), pch = c(16,16),
                       main = paste(as.character(SigReg[rr, c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)

        }

      }

      return(GLMgene$allRegulators[GLMgene$significantRegulators, -6])
    }

  }


  ## GENE = NULL

  if (is.null(gene)) {  ### Plot all genes regulated by the regulator

    SigniReguGene = GetPairsGeneRegulator(genes = NULL, getGLMoutput = GLMoutput)
    SigniReguGene = SigniReguGene[SigniReguGene[,"regulator"] == regulator,]
    myomics = SigniReguGene[,"omic"]
    myomic = unique(myomics)

    if (nrow(SigniReguGene) > 0) {  # When there are genes regulated by this regulator

      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(GLMoutput$arguments$dataOmics[[myomic]][regulator,])
      }

      numGenes = length(SigniReguGene$gene)
      cat(paste(numGenes, "genes are regulated by", regulator)); cat("\n")

      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)

        if (order) reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni

        lapply(1:numGenes, function (i) {

          geneValues = GLMoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$Y$y
          if (order) geneValues = geneValues[myorder]
          errorValues = getErrorValues(geneValues, myreplicates)
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[myrepliUni]
          names(geneValues) = myrepliUni

          x.points = 1:length(myrepliUni)
          eje = myrepliUni

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col[myomics[i]]), xlab = xlab,
                       yylab = c(SigniReguGene[i,"gene"], regulator), pch = c(16,16),
                       main = paste(as.character(SigniReguGene[1,c("omic", "area")]), collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        })
      } else { cat("Regulator values could not be recovered from GLMoutput. Please provide them in reguValues argument to generate the plot.\n") }

      return(SigniReguGene$gene)

    } else { cat(paste("There are no genes significantly regulated by", regulator)); cat("\n") }

  }


  ## GENE + REGULATOR

  if (!is.null(gene) && !is.null(regulator)) {  ### Plot only the given gene and the given regulator

    geneResults = GLMoutput$ResultsPerGene[[gene]]

    if (is.null(geneResults)) {
      stop(paste("No GLM was obtained for gene", gene))
    } else
    {
      myomic = geneResults$allRegulators[regulator, "omic"]

      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(GLMoutput$arguments$dataOmics[[myomic]][regulator,]) # regulator values
      }

      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)

        if (order)  reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni

        geneValues = GLMoutput$ResultsPerGene[[gene]]$Y$y
        if (order) geneValues = geneValues[myorder]
        errorValues = getErrorValues(geneValues, myreplicates)
        geneValues = tapply(geneValues, myreplicates, mean)
        geneValues = geneValues[myrepliUni]
        names(geneValues) = myrepliUni

        x.points = 1:length(myrepliUni)
        eje = myrepliUni

        plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                     col = c(gene.col, any.col[myomic]), xlab = xlab,
                     yylab = c(gene, regulator), pch = c(16,16),
                     main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                  collapse = " "),
                     numLines = numLines, x.names = eje,
                     geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)

      } else {

        cat("Regulator values could not be recovered from GLMoutput.\n")

        regulator = geneResults$allRegulators[regulator,"filter"]

        if (regulator %in% rownames(geneResults$allRegulators)) {

          cat(paste(regulator, "values will be plotted instead.")); cat("\n")
          cat(paste(regulator, "summarizes information from the following correlated regulators:")); cat("\n")
          cat(geneResults$allRegulators[geneResults$allRegulators[,"filter"] == regulator,"regulator"]); cat("\n")

          reguValues = geneResults$X
          reguValues = reguValues[, colnames(reguValues) == regulator]
          if (order) reguValues = reguValues[myorder]
          errorValuesRegu = getErrorValues(reguValues, myreplicates)

          reguValues = tapply(reguValues, myreplicates, mean)
          reguValues = reguValues[myrepliUni]
          names(reguValues) = myrepliUni

          geneValues = GLMoutput$ResultsPerGene[[gene]]$Y$y
          if (order) geneValues = geneValues[myorder]
          errorValues = getErrorValues(geneValues, myreplicates)
          geneValues = tapply(geneValues, myreplicates, mean)
          geneValues = geneValues[myrepliUni]
          names(geneValues) = myrepliUni

          x.points = 1:length(myrepliUni)
          eje = myrepliUni

          plotGeneRegu(x.points = x.points, geneValues = geneValues, reguValues = reguValues,
                       col = c(gene.col, any.col[myomic]), xlab = xlab,
                       yylab = c(gene, regulator), pch = c(16,16),
                       main = paste(as.character(geneResults$allRegulators[regulator, c("omic", "area")]),
                                    collapse = " "),
                       numLines = numLines, x.names = eje,
                       geneErrorValues = errorValues, reguErrorValues = errorValuesRegu)
        } else {
          cat("The selected regulator was not declared as significant by the GLM.\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }

      }

    }

  }

}





# Function to obtain all significant pairs gene-regulator per omic --------

# For all genes
GetPairsGeneRegulator = function (genes = NULL, getGLMoutput) {

  if (is.null(genes)) genes = rownames(getGLMoutput$GlobalSummary$ReguPerGene)

  myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, getGLMoutput))

  #   colnames(myresults) = c("gene", "regulator", "omic", "area")
  return(myresults)
}


RegulationPerCondition = function(getGLMoutput){
  # getGLMoutput: results of the getGLM function.
  design = getGLMoutput$arguments$finaldesign
  Group = getGLMoutput$arguments$groups

  # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
  genes = rownames(getGLMoutput$GlobalSummary$ReguPerGene)
  myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, getGLMoutput))
  colnames(myresults) = c(colnames(myresults)[1:4], "representative")
  myresults[myresults[, "representative"] == "Model", "representative"] = ""

  if (is.null(design)){

    # Anyado la columna de coeficientes.
    coeffs = matrix(1, nrow(myresults), 1)
    colnames(coeffs) = "coefficients"
    rownames(coeffs) = rownames(myresults)
    myresults = cbind(myresults, coeffs)
    myresults[grep("_N", myresults[, "representative"]), "coefficients"] = -1  # Para cambiar el signo si pertenece al grupo de correlacionados negativamente

    for(k in unique(myresults[,"gene"])){

      # Posicion y reguladores que son representantes.
      counts = grep("_R", myresults[myresults[,"gene"] == k, "representative"]) # positions of representatives of mc
      representatives = myresults[myresults[,"gene"] == k, "regulator"][counts]      # Devuelve el nombre real de los reguladores representantes
      omic.representative = myresults[myresults[,"gene"] == k, c("regulator", "representative")][counts,]   # Columna Regulator y Representative

      # Necesito el if, si no da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante
      if(length(representatives) != 0){
        norow.nulls = which(myresults[myresults[,"gene"] == k, "representative"] != "")
        myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
        myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
        myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])

        for(i in 1:length(representatives)){
          # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante
          reg.rep = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == representatives[i], "representative"]
          myresults[myresults[,"gene"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
        }

        # Reguladores significativos del GLM. Pongo gsub() porque haciendo pruebas he visto que hay reguladores que se nombran `nombre regulador`.
        # Las comitas haran que no pueda encontrar el regulador
        # en la tabla. Sin embargo, creo sign.glm para manterner las comitas y poder acceder a la tabla de coeficientes
        significatives = gsub("`", "", names(getGLMoutput$ResultsPerGene[[k]]$coefficients[2:nrow(getGLMoutput$ResultsPerGene[[k]]$coefficients), 1]))
        sign.glm = names(getGLMoutput$ResultsPerGene[[k]]$coefficients[2:nrow(getGLMoutput$ResultsPerGene[[k]]$coefficients), 1])

        for(i in 1:length(significatives)){
          if(any(significatives[i] == omic.representative[,2])){
            # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo GLM.
            index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
            PN = myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*getGLMoutput$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
          } else {
            # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
            myresults[myresults[,"gene"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = getGLMoutput$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]
          }
        }

      } else {
        # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
        myresults[myresults[,"gene"] == k, "coefficients"] = getGLMoutput$ResultsPerGene[[k]]$coefficients[2:nrow(getGLMoutput$ResultsPerGene[[k]]$coefficients), 1]
      }
    }


  } else {

    # Anyado las columnas de las condiciones experimentales. Pongo "Group" porque al hacer model.matrix() siempre coloca "Group" y lo que se almacena en el objeto Group
    index = unique(Group)
    names.groups = paste("Group", index, sep = "")
    conditions = matrix(0, nrow(myresults), length(names.groups))
    colnames(conditions) = names.groups
    rownames(conditions) = rownames(myresults)
    myresults = cbind(myresults, conditions)

    for(k in unique(myresults[,"gene"])){
      significant.regulators = getGLMoutput$ResultsPerGene[[k]]$significantRegulators                    # Reguladores significativos.
      model.variables = gsub("`", "", rownames(getGLMoutput$ResultsPerGene[[k]]$coefficients))[-1]       # Reguladores e interacciones en el modelo.

      # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
      interactions.model = gsub("`", "", rownames(getGLMoutput$ResultsPerGene[[k]]$coefficients)[grep(":", rownames(getGLMoutput$ResultsPerGene[[k]]$coefficients))])

      inter.variables = unlist(strsplit(interactions.model, ":", fixed = TRUE))
      if(is.null(inter.variables)){
        inter.variables = NULL                                                                            # No hay interacciones.
      } else {
        inter.variables = inter.variables[seq(2, length(inter.variables), by = 2)]                        # Reguladores que presentan interseccion con algun grupo.
      }

      variables.only = setdiff(setdiff(model.variables, interactions.model), inter.variables)             # Reguladores solos en el modelo, sin interacciones.
      if(length(grep("Group", variables.only)) != 0){                                                     # No puedo hacer la interseccion con las variables significativas porque me cargo tambien omic_mc: hay que eliminar Group si esta.
        variables.only = variables.only[-grep("Group", variables.only)]
      }

      variables.inter.only = intersect(inter.variables, model.variables)                                  # Reguladores con interaccion y solas.
      variables.inter = setdiff(inter.variables, model.variables)                                         # Reguladores con solo interaccion (no aparecen solas en el modelo).

      for(j in 2:nrow(getGLMoutput$ResultsPerGene[[k]]$coefficients)){
        regul = unlist(strsplit(gsub("`", "", rownames(getGLMoutput$ResultsPerGene[[k]]$coefficients)[j]), ":"))

        # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
        if(any(regul %in% variables.only)){
          if(any(regul %in% significant.regulators)){
            myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
          } else {
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] = getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
          }
        }

        if(any(regul %in% variables.inter)){
          if(any(regul %in% significant.regulators)){
            myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] + getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
          } else {
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] + getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
          }
        }

        if(any(regul %in% variables.inter.only)){
          if(any(regul %in% significant.regulators)){
            if(length(regul) == 1){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] + getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] + getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
            }
          } else {
            if(length(regul) == 1){
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] + getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] + getGLMoutput$ResultsPerGene[[k]]$coefficients[j,]
            }
          }
        }
      }

      # Veo si hay representantes, en caso de haberlos asignara la misma fila del representante a los reguladores que acaben en "_P" y el opuesto a los que acaban en "_N".
      countsR = grep("_R", myresults[myresults[,"gene"] == k, "representative"])

      if(length(countsR) != 0){
        countsR = myresults[myresults[,"gene"] == k, 5:ncol(myresults)][countsR,]

        # Para los correlacionados positivamente: mete la misma fila de coeficientes del representante.
        countsP = countsR
        countsP[,"representative"] = sub("_R", "", countsP[,"representative"])
        countsP[,"representative"] = paste(countsP[,"representative"], "_P", sep = "")

        for(l in 1:nrow(countsP)){
          myresults[myresults[,"gene"] == k & myresults[,"representative"] == countsP[l,"representative"], 6:ncol(myresults)] = countsP[l,2:ncol(countsP)]
        }

        # Para los correlacionados negativamente: mete la fila opuesta de coeficientes del representante.
        countsN = countsR
        countsN[,"representative"] = sub("_R", "", countsN[,"representative"])
        countsN[,"representative"] = paste(countsN[,"representative"], "_N", sep = "")

        for(l in 1:nrow(countsN)){
          myresults[myresults[,"gene"] == k & myresults[,"representative"] == countsN[l,"representative"], 6:ncol(myresults)] = -countsN[l,2:ncol(countsN)]
        }

        counts = grep("_R", myresults[myresults[,"gene"] == k, "representative"])
        representatives = myresults[myresults[,"gene"] == k, "regulator"][counts]                                    # Devuelve el nombre real de los reguladores representantes
        omic.representative = myresults[myresults[,"gene"] == k, c("regulator", "representative")][counts,]          # Columna Regulator y Representative

        # Necesito el if, sino da error. En caso de entrar, elimino las coletillas para que sea mas sencillo buscar y asignar el representante.
        if(length(representatives) != 0){
          norow.nulls = which(myresults[myresults[,"gene"] == k, "representative"] != "")
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_P", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_N", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])
          myresults[myresults[,"gene"] == k, "representative"][norow.nulls] = sub("_R", "", myresults[myresults[,"gene"] == k, "representative"][norow.nulls])

          for(i in 1:length(representatives)){
            # Aquellos que se llamen igual "omica_mc(numero)", se les asignara el representante.
            reg.rep = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == representatives[i], "representative"]
            myresults[myresults[,"gene"] == k & myresults[,"representative"] == reg.rep, "representative"] = representatives[i]
          }
        }
      }
    }
  }
  myresults[,6:ncol(myresults)] = signif(myresults[,6:ncol(myresults)], digits = 4) # Para que no salgan los numeros en diferentes notaciones
  return(myresults)
}






# For only 1 gene
GetPairs1GeneAllReg = function (gene, getGLMoutput) {

  reguSignif = getGLMoutput$ResultsPerGene[[gene]]$significantRegulators

  if (is.null(reguSignif)) {  # NO significant regulators
    return (NULL)

  } else {  # Significant regulators

    reguSignif = getGLMoutput$ResultsPerGene[[gene]]$allRegulators[reguSignif,]
    reguSignif = reguSignif[,c("gene", "regulator", "omic", "area", "filter")]
    return (reguSignif)
  }
}



BetaTest = function(coeffs, myGene, MOREresults){

  family = MOREresults$arguments$family
  alfa = MOREresults$arguments$alfa

  # Coeficientes distintos de 0 en todas las columnas: significa que tendran una suma y debera hacerse el test.
  mycoeffs = coeffs[coeffs[,"gene"] == myGene,]
  #row.null =  apply(mycoeffs[, 6:ncol(mycoeffs)], 1, function(x) all(x != 0))   ###  no seria en la columna 6?????

  # Tabla que solo recoge aquellos reguladores que poseen una suma en sus coeficientes.
  mytable = mycoeffs[which(mycoeffs[,6] != 0), c(1,2,5:ncol(mycoeffs))]

  # Si todos los reguladores tienen interseccion (no aparecen solos), no habra suma de coeficientes y dara fallo.
  # Pongo el siguiente if para controlar esto.
  if(dim(mytable)[1] != 0){

    # Comparo que reguladores tienen coeficientes distintos en cada grupo: eso significara que ha habido interaccion y aparece solo en el modelo,
    # entonces habra una suma de coeficientes.
    # El primero es el referencial, comparo con el. Lo he intentado con apply, devuelve cada elemento si es TRUE o FALSE y yo quiero la fila.
    # También con which() pero no devuelve
    # la posicion de la fila.
    betas = NULL
    for(j in 1:nrow(mytable)){
      if(any(mytable[j,4] != mytable[j, 5:ncol(mytable)])){betas = rbind(betas, mytable[j,])}
    }

    # Hay suma de coeficientes.
    if(!is.null(betas)){
      mySignificatives = rownames(MOREresults$ResultsPerGene[[myGene]]$coefficients)[-1]
      mySigni = gsub("`", "", mySignificatives)
      myY = MOREresults$ResultsPerGene[[myGene]]$Y[,1]
      myX = MOREresults$ResultsPerGene[[myGene]]$X
      myX = myX[, mySigni]

      for(i in 1:nrow(betas)){

        # Para tener en cuenta si es omica_mc_R.
        if(betas[i,"representative"] == ""){
          myRegulator = betas[i, "regulator"]
        } else {
          myRegulator = betas[i, "representative"]
        }

        aquitar = mySignificatives[grep(myRegulator, mySignificatives)]
        group = unlist(strsplit(aquitar, ":", fixed = TRUE))
        group = gsub("`",  "",group[grep("Group", group)])

        # Modelo GLM
        mymodel = glm(myY~., data = myX, family = family)
        myHyp = paste0(paste(aquitar, collapse = "+"), " = 0")
        pvalue = try(linearHypothesis(mymodel, c(myHyp), test = "Chisq")$"Pr(>Chisq)"[2], silent = TRUE)

        if(class(pvalue) == "try-error"){pvalue = NA}

        if(pvalue > alfa){
          # Para escoger la fila correcta segun el nombre del regulador.
          if(betas[i,"representative"] == ""){
            coeffs[coeffs[, "gene"] == myGene & coeffs[, "regulator"] == myRegulator, group] = 0
          } else {
            coeffs[coeffs[, "gene"] == myGene & coeffs[, "representative"] == myRegulator, group] = 0
          }
        }
      }
    }
  }
  return(coeffs)
}








# Plot 1 gene versus 1 regulator ------------------------------------------

plotGeneRegu = function (x.points, geneValues, reguValues, geneErrorValues, reguErrorValues, col = c(1,2),
                         xlab = "", yylab = c("right", "left"), pch = c(16,17), main = "",
                         numLines = NULL, x.names = NULL, yleftlim, yrightlim) {

  # Adjust the axis to include the error value
  if (missing(yrightlim)) {
    if (! missing(geneErrorValues) && ! is.null(geneErrorValues)) {
      yrightlim = range(c(geneValues - geneErrorValues, geneValues + geneErrorValues), na.rm = TRUE)
    } else {
      yrightlim = range(geneValues, na.rm = TRUE)
    }
  }

  if (missing(yleftlim)) {
    if (! missing(reguErrorValues) && ! is.null(reguErrorValues)) {
      yleftlim = range(c(reguValues - reguErrorValues, reguValues + reguErrorValues), na.rm = TRUE)
    } else {
      yleftlim = range(reguValues, na.rm = TRUE)
    }
  }

  plot.y2(x = x.points, yright = geneValues, yleft = reguValues, yleftlim = yleftlim,
          col = col, xlab = xlab, yylab = yylab, pch = pch, main = main, yrightlim = yrightlim,
          yrightErrorValues = geneErrorValues, yleftErrorValues = reguErrorValues)

  if (!is.null(numLines)) {
    for (aa in numLines) {
      abline(v = aa, lty = 2, col = 1)
    }
  }

  if (!is.null(x.names)) {
    axis(side=1, at = x.points, labels = x.names, cex.axis = 0.8, las=2)
  }
}








# Plot Y2 -----------------------------------------------------------------

# By Ajay Shah (taken from [R] Plot 2 time series with different y axes (left and right),
# in https://stat.ethz.ch/pipermail/r-help/2004-March/047775.html)

# Modified by: Sonia Tarazona

### PARAMETERS (default):
# x: data to be drawn on X-axis
# yright: data to be drawn on Y right axis
# yleft: data to be drawn on Y left axis
# yrightlim (range(yright, na.rm = TRUE)): ylim for rigth Y-axis
# yleftlim (range(yleft, na.rm = TRUE)): ylim for left Y-axis
# xlab (NULL): Label for X-axis
# yylab (c("","")): Labels for right and left Y-axis
# pch (c(1,2)): Type of symbol for rigth and left data
# col (c(1,2)): Color for rigth and left data
# linky (TRUE): If TRUE, points are connected by lines.
# smooth (0): Friedman's super smoothing
# lwds (1): Line width for smoothed line
# length (10): Number of tick-marks to be drawn on axis
# ...: Other graphical parameters to be added by user (such as main, font, etc.)
###


plot.y2 <- function(x, yright, yleft, yrightlim = range(yright, na.rm = TRUE),
                    yleftlim = range(yleft, na.rm = TRUE),
                    xlim = range(x, na.rm = TRUE),
                    xlab = NULL, yylab = c("",""), lwd = c(2,2),
                    pch = c(1,2), col = c(1,2), type = c("o","o"),
                    linky = TRUE, smooth = 0, bg = c("white","white"),
                    lwds = 1, length = 10, ...,
                    x2 = NULL, yright2 = NULL, yleft2 = NULL, col2 = c(3,4),
                    yrightErrorValues, yleftErrorValues
)
{
  #par(mar = c(5,2,4,2), oma = c(0,3,0,3))

  ## Plotting RIGHT axis data

  plot(x, yright, axes = FALSE, ylab = "", xlab = xlab, ylim = yrightlim,
       xlim = xlim, pch = pch[1], type = type[1], lwd = lwd[1],
       col = col[1], ...)

  axis(4, pretty(yrightlim, length), col = 1, col.axis = 1)

  if (is.null(yright2) == FALSE) {
    points(x2, yright2, type = type[1], pch = pch[1], lwd = lwd[1], col = col2[1], ...)
  }

  #if (linky) lines(x, yright, col = col[1], ...)

  if (smooth != 0) lines(supsmu(x, yright, span = smooth), col = col[1], lwd = lwds, ...)

  if(yylab[1]=="") {
    mtext(deparse(substitute(yright)), side = 4, outer = FALSE, line = 2,
          col = col[1], cex = 0.9,...)
  } else {
    mtext(yylab[1], side = 4, outer = FALSE, line = 2, col = col[1], cex = 0.9,...)
  }

  # Plot arrows showing standard error
  if (!missing(yrightErrorValues) && ! is.null(yrightErrorValues)) {
    arrows(x, yright - yrightErrorValues, x, yright + yrightErrorValues,
           code = 3, length = 0.02, angle = 90, col = col[1])
  }


  par(new = T)

  ## Plotting LEFT axis data
  plot(x, yleft, axes = FALSE, ylab = "" , xlab = xlab, ylim = yleftlim,
       xlim = xlim, bg = bg[1],
       pch = pch[2], type = type[2], lwd = lwd[2], col = col[2], ...)

  box()

  axis(2, pretty(yleftlim, length), col = 1, col.axis = 1)

  if (is.null(yleft2) == FALSE) {
    points(x2, yleft2, type = type[2], pch = pch[2], bg = bg[2],
           lwd = lwd[2], col = col2[2], ...)
  }


  #if (linky) lines(x, yleft, col = col[2], ...)

  if (smooth != 0) lines(supsmu(x, yleft, span = smooth), col = col[2], lwd=lwds, ...)

  if(yylab[2] == "") {
    mtext(deparse(substitute(yleft)), side = 2, outer = FALSE, line = 2, col = col[2], cex = 0.9, ...)
  } else {
    mtext(yylab[2], side = 2, outer = FALSE, line = 2, col = col[2], cex = 0.9, ...)
  }

  if (!missing(yleftErrorValues) && ! is.null(yleftErrorValues)) {
    arrows(x, yleft - yleftErrorValues, x, yleft + yleftErrorValues,
           code = 3, length = 0.02, angle = 90, col = col[2])
  }

  ## X-axis
  ##  axis(1, at = pretty(xlim, length))  ## Comment last line


}
