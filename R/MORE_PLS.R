#########################################################################################
######           Functions to integrate omics data using PLSs                      ######
#########################################################################################


## By Maider
## 20-July-2023
## Last modified: July 2023


options(stringsAsFactors = FALSE)

# library(maSigPro)  ### ver que pasa si no cargamos maSigPro
library(igraph)
library(MASS)
library(psych)
library(car)
library(Hmisc)
library(ropls)
library(fastDummies)
# Partial Least Squares  -----------------------------------------------
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
#'\code{GetPLS} fits a PLS model for all the genes in the data set to identify
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
#' @param alfa Significance level. By default, 0.05.
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
#' @return
#' @export
#'
#' @examples
GetPLS = function(GeneExpression,
                  associations,
                  data.omics,
                  clinic = NULL,
                  edesign = NULL,
                  center = TRUE, scale = TRUE,
                  epsilon = 0.00001,
                  alfa = 0.05, 
                  interactions.reg = TRUE,
                  min.variation = 0,
                  min.obs = 10,
                  omic.type = 0,
                  clinic.type = 0,
                  scaletype = 'auto',
                  p.method ='jack',
                  vip = 0.8){
  
  cont.var = NULL
  # Converting matrix to data.frame
  GeneExpression = as.data.frame(GeneExpression)
  data.omics = lapply(data.omics, as.data.frame)
  if(!is.null(clinic)){
    clinic = as.data.frame(clinic)
    data.omics = c(list(clinic = clinic),data.omics)
    
    ## clinic.type
    if (length(clinic.type) == 1) clinic.type = rep(clinic.type, ncol(clinic))
    if (length(clinic.type) != nrow(clinic))  stop("ERROR: clinic.type must be a vector with length equal to the number of clinical variables or a number (0 or 1) if all the clinical variables are of the same type. ")
    
    #Add in associations clinic to consider all the clinical variables in all genes
    if(!is.null(associations)){associations = c(list(clinic = NULL),associations)}
    
    #Add information to omic.type even if it is not relevant
    omic.type = c(0,omic.type)
    
  }else{clinic.type=NULL}
  
  # If associations is NULL create a list of associations NULL
  if (is.null(associations)){
    associations=vector('list',length(data.omics))
    names(associations)=names(data.omics)
  }
  
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
  
  ##Verify that the scaletype selected is a possible selection
  
  if(!scaletype %in% c('auto','pareto','block')){stop('ERROR: The selected type of scaling is not one of auto, pareto or block')}
  if(!p.method %in% c('pval','jack')){stop('ERROR: The selected method for p value calculation is not one of pval or jack')}
  
  ## Checking that samples are in the same order in GeneExpressionDE, data.omics and edesign
  orderproblem<-FALSE
  if(is.null(edesign)){
    nameproblem<-!all(sapply(data.omics, function(x) length(intersect(colnames(x),colnames(GeneExpression))==ncol(GeneExpression))))
    if(nameproblem){
      cat('Warning. GeneExpression and data.omics samples have not same names. We assume that they are ordered.\n')
    }else{
      orderproblem<-!all(sapply(data.omics, function(x) identical(colnames(x),colnames(GeneExpression))))
      if(orderproblem){
        data.omics<-lapply(data.omics, function(x) x[,colnames(GeneExpression)])
      }
    }
    
  } else{
    nameproblem<-!all(c(sapply(data.omics, function(x) length(intersect(colnames(x),colnames(GeneExpression)))==ncol(GeneExpression)), length(intersect(rownames(edesign),colnames(GeneExpression)))==ncol(GeneExpression)))
    if(nameproblem){
      cat('Warning. GeneExpression, edesign and data.omics samples have not same names. We assume that they are ordered.\n')
    } else{
      orderproblem<-!all(c(sapply(data.omics, function(x) identical(colnames(x),colnames(GeneExpression))), identical(colnames(GeneExpression),rownames(edesign))))
      if(orderproblem){
        data.omics<-lapply(data.omics, function(x) x[,colnames(GeneExpression)])
        edesign<-edesign[colnames(GeneExpression), , drop=FALSE]
      }
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
    
    #Change the name in the association matrix
    if(!is.null(associations[[i]])){
      associations[[i]][[2]]=gsub(':', '-', associations[[i]][[2]])
    }
    
    if(length(problemas) > 0) {
      cat("In",names(data.omics)[i], problemas ,"regulators have names that may cause conflict with the algorithm by ending in _R, _P or _N", "\n")
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
  }
  des.mat = ScalePLSdesmat(edesign, scaletype, center, scale)
  
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
  
  tmp = LowVariationRegu(min.variation, data.omics, Group, associations, Allgenes, omic.type, clinic.type)
  data.omics = tmp[["data.omics"]]
  associations = tmp[["associations"]]
  myregLV = tmp[["myregLV"]]
  rm("tmp"); gc()
  
  if(all(sapply(data.omics, function(x)nrow(x)==0))) stop("ERROR: No regulators left after LowVariation filter. Consider being less restrictive.")
  
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
  colnames(GlobalSummary$GoodnessOfFit) = c( "RsquaredY", "Qsquared",'RMSEE',"sigReg")
  
  GlobalSummary$ReguPerGene = matrix(0, ncol = 3*length(data.omics), nrow = nGenes)
  rownames(GlobalSummary$ReguPerGene) = Allgenes
  colnames(GlobalSummary$ReguPerGene) = c(paste(names(data.omics), "Ini", sep = "-"),
                                          paste(names(data.omics), "Mod", sep = "-"),
                                          paste(names(data.omics), "Sig", sep = "-"))
  
  ## Specific results for each gene
  ResultsPerGene=vector("list", length=length(Allgenes))
  names(ResultsPerGene) = Allgenes
  
  ### Computing model for each gene
  cat("Checking multicollinearity, selecting predictors and fitting model for ...\n")
  
  pap = c(1, 1:round(nGenes/100) * 100, nGenes)
  
  for (i in 1:nGenes) {
    
    gene=Allgenes[i]
    print(gene)
    cat(i)
    ResultsPerGene[[i]] = vector("list", length = 5)
    names(ResultsPerGene[[i]]) = c("Y", "X", "coefficients", "allRegulators", "significantRegulators")
    
    if (is.element(i, pap)) cat(paste("Fitting model for gene", i, "out of", nGenes, "\n"))
    
    RetRegul = GetAllReg(gene=gene, associations=associations, data.omics = data.omics)
    RetRegul.gene = RetRegul$Results  ## RetRegul$TableGene: nr reg per omic
    ## Some of these reg will be removed, because they are not in data.omics
    
    
    # RetRegul.gene--> gene/regulator/omic/area
    RetRegul.gene=RetRegul.gene[RetRegul.gene[,"regulator"]!= "No-regulator", ,drop=FALSE] ## Remove rows with no-regulators
    
    
    ### NO INITIAL REGULATORS
    if(length(RetRegul.gene)==0){ ## En el caso de que no hayan INICIALMENTE reguladores -> Calcular modelo con variables experimentales
      
      if (is.null(edesign)) {
        ResultsPerGene[[i]]$X = NULL
        ResultsPerGene[[i]]$significantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
        myPLS = NULL
        
      } else {
        des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
        colnames(des.mat2)[1] = "response"
        des.mat2 = na.omit(des.mat2)
        
        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        des.mat2 = des.mat2[,sdNo0]
        
        myPLS = NULL
        
        ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
        ResultsPerGene[[i]]$significantRegulators = NULL
        ResultsPerGene[[i]]$allRegulators = NULL
      }
      
      # GlobalSummary$ReguPerGene  # this is initially set to 0 so no need to modify it
      
      ### WITH INITIAL REGULATORS
    } else { ## There are regulators for this gene at the beginning
      
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
          ResultsPerGene[[i]]$significantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene
          myPLS = NULL
          
        } else {
          des.mat2 = cbind(t(GeneExpression[gene,]), des.mat)
          colnames(des.mat2)[1] = "response"
          des.mat2 = na.omit(des.mat2)
          
          # Removing predictors with constant values
          sdNo0 = apply(des.mat2, 2, sd)
          sdNo0 = names(sdNo0)[sdNo0 > 0]
          des.mat2 = des.mat2[,sdNo0]
          
          myPLS = NULL
          
          GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                             data.frame("gene" = gene,
                                                        "problem" = 'No regulators left after NA/LowVar filtering'))
          
          ResultsPerGene[[i]]$X = des.mat2[,-1, drop = FALSE]
          ResultsPerGene[[i]]$significantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene
        }
        
      } else {  ## Regulators for the model!!
        ResultsPerGene[[i]]$allRegulators = res$SummaryPerGene
        ## Make the groups of omics
        regupero = lapply(unique(res$SummaryPerGene[,'omic']), function(x) rownames(res$SummaryPerGene)[res$SummaryPerGene[,'omic'] == x & res$SummaryPerGene[,'filter'] == "Model"])
        names(regupero) = unique(res$SummaryPerGene[,'omic'])
        
        ## Create the interactions between regulators and edesign 
        
        des.mat2 = RegulatorsInteractionsPLS(interactions.reg, reguValues = res$RegulatorMatrix,
                                             edesign, clinic.type, cont.var, GeneExpression, gene, regupero, omic.type)
        
        #Add the ones related to the interactions
        regupero = lapply(names(regupero), function(x) if(length(regupero[[x]])!=0){colnames(des.mat2)[grepl(paste(regupero[[x]],collapse = "|"), colnames(des.mat2))]})
        names(regupero) = unique(res$SummaryPerGene[,'omic'])
        
        res$RegulatorMatrix = ScalePLS(des.mat2[,-1,drop=FALSE], regupero, omic.type, scaletype, center, scale)
        
        #Use them jointly
        des.mat2 = data.frame(des.mat2[,1,drop=FALSE], des.mat, res$RegulatorMatrix,check.names = FALSE)
        
        # Removing predictors with constant values
        sdNo0 = apply(des.mat2, 2, sd)
        sdNo0 = names(sdNo0)[sdNo0 > 0]
        quito = setdiff(colnames(des.mat2), sdNo0)
        des.mat2 = des.mat2[,sdNo0]
        sdNo0 = unique(gsub(".*[:?](.*)$", "\\1", sdNo0))
        quito = unique(gsub(".*[:?](.*)$", "\\1", quito))
        quito = setdiff(quito, sdNo0)
        ResultsPerGene[[i]]$allRegulators[quito,"filter"] = "Constant"
        
        #Include the gene expression information. Necesario cuando no considerabamos RegulatorsInteractions. Pero lo hace internamente no es neceario en este caso
        #des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
        #colnames(des.mat2)[1] = 'response'
        
        ## Computing GLM model
        if (nrow(des.mat2)<7){cross = nrow(des.mat)-2}else{cross =7}
        myPLS = try(suppressWarnings( opls(des.mat2[,-1], scale(des.mat2[,1], center = center,scale=scale), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0)),silent = TRUE)
        
        if(class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
          myPLS = try(suppressWarnings( opls(des.mat2[,-1], scale(des.mat2[,1], center =center, scale = scale), info.txtC = 'none', fig.pdfC='none', scaleC = 'none', crossvalI = cross, permI=0, predI=1)), silent = TRUE)
        }
        
        if (class(myPLS)=='try-error' || length(myPLS@modelDF)==0){
          
          myPLS = NULL
          
          GlobalSummary$GenesNOmodel = rbind(GlobalSummary$GenesNOmodel,
                                             data.frame("gene" = gene, "problem" = "No significant components on PLS"))
          
          ## Extracting significant regulators
          ResultsPerGene[[i]]$significantRegulators = NULL
          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Sig" = NA, stringsAsFactors = FALSE)
          
          ## Counting original regulators in the model per omic
          contando = ResultsPerGene[[i]]$allRegulators[which(ResultsPerGene[[i]]$allRegulators[,"filter"] == "Model"),]
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
          
          ## Counting significant regulators per omic
          GlobalSummary$ReguPerGene[gene, grep("-Sig", colnames(GlobalSummary$ReguPerGene))] = NA
          
        }else{
          
          if (p.method == 'jack'){
            pval = p.valuejack(myPLS, des.mat2, alfa)
          } else {
            pval = p.coef(myPLS, 100, des.mat2)
          }
          
          
          #ResultsPerGene[[i]]$significantRegulators  = intersect(colnames(res$RegulatorMatrix), colnames(des.mat2[,which(pvalor<alfa & myPLS@vipVn>0.5)]) )
          
          #Tratar como significativas tan solo las que cumplan ambas condiciones
          sigvariables = intersect(names(myPLS@vipVn[which(myPLS@vipVn>vip)]), rownames(pval)[which(pval<alfa)])
          ResultsPerGene[[i]]$coefficients = data.frame('coefficient' = myPLS@coefficientMN[sigvariables,], 'p-value' = pval[sigvariables,])
          rows_to_remove = rownames(ResultsPerGene[[i]]$coefficients)[grepl("_0$", rownames(ResultsPerGene[[i]]$coefficients)) & !rownames(ResultsPerGene[[i]]$coefficients) %in% colnames(des.mat)]
          ResultsPerGene[[i]]$coefficients = ResultsPerGene[[i]]$coefficients[!rownames(ResultsPerGene[[i]]$coefficients) %in% rows_to_remove, ]
          
          # Obtain the indices of the rows to modify
          rows_to_modify = grepl("_1$", rownames(ResultsPerGene[[i]]$coefficients)) & !(rownames(ResultsPerGene[[i]]$coefficients) %in% colnames(des.mat))
          rownames(ResultsPerGene[[i]]$coefficients)[rows_to_modify] =  gsub("_1$", "", rownames(ResultsPerGene[[i]]$coefficients)[rows_to_modify])
          
          ## Extracting significant regulators and recovering correlated regulators
          myvariables = unlist(strsplit(sigvariables, ":", fixed = TRUE))
          #Eliminar _0 y _1 correspondientes a los dummy
          myvariables = gsub('_0$','',myvariables)
          myvariables = gsub('_1$','',myvariables)
          myvariables = intersect(myvariables, rownames(ResultsPerGene[[i]]$allRegulators))
          
          ResultsPerGene[[i]]$allRegulators = data.frame(ResultsPerGene[[i]]$allRegulators, "Sig" = 0, stringsAsFactors = FALSE)
          ResultsPerGene[[i]]$allRegulators[myvariables, "Sig"] = 1
          
          ResultsPerGene[[i]]$significantRegulators = myvariables
          
          ## Counting original regulators in the model per omic              
          contando = ResultsPerGene[[i]]$allRegulators
          quitar = which(contando[,"filter"] == "MissingValue")
          if (length(quitar) > 0) contando = contando[-quitar,]
          quitar = which(contando[,"filter"] == "LowVariation")
          if (length(quitar) > 0) contando = contando[-quitar,]
          
          contando = table(contando[,"omic"])
          contando = as.numeric(contando[names(data.omics)])
          contando[is.na(contando)] = 0
          GlobalSummary$ReguPerGene[gene, grep("-Mod", colnames(GlobalSummary$ReguPerGene))] = contando
          
          ## Counting significant regulators per omic
          if (length(ResultsPerGene[[i]]$significantRegulators) > 0) {
            contando = ResultsPerGene[[i]]$allRegulators[ResultsPerGene[[i]]$significantRegulators,]
            contando = table(contando[,"omic"])
            contando = as.numeric(contando[names(data.omics)])
            contando[is.na(contando)] = 0
            GlobalSummary$ReguPerGene[gene, grep("-Sig", colnames(GlobalSummary$ReguPerGene))] = contando
          }
          
          
        }
        
        
        
      } ## Close "else" --> None regulators from begining
      
      if (is.null(myPLS)) {
        
        ResultsPerGene[[i]]$Y = GeneExpression[i,]
        ResultsPerGene[[i]]$coefficients = NULL
        
        GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[rownames(GlobalSummary$GoodnessOfFit) != gene,]
        
        
      } else {
        ResultsPerGene[[i]]$Y = data.frame("y" = myPLS@suppLs$y, "fitted.y" = myPLS@suppLs$yPreMN,
                                           "residuals" = residuals(myPLS))
        
        
        GlobalSummary$GoodnessOfFit[gene,] = c(myPLS@modelDF[,'R2Y(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@modelDF[,'Q2(cum)'][myPLS@summaryDF[,'pre']],
                                               myPLS@summaryDF[,'RMSEE'][myPLS@summaryDF[,'pre']],
                                               as.integer(length(ResultsPerGene[[i]]$significantRegulators)))
        
      }
      
    }
    
  }  ## At this point the loop for all genes is finished
  end_t =Sys.time()
  
  genesNosig = names(which(GlobalSummary$GoodnessOfFit[,4]==0))
  genessig = setdiff(rownames(GlobalSummary$GoodnessOfFit), genesNosig)
  GlobalSummary$GoodnessOfFit = GlobalSummary$GoodnessOfFit[genessig,]
  
  myarguments = list(edesign = edesign, finaldesign = des.mat, groups = Group, alfa = alfa, 
                     center = center, scale = scale, clinic.type = clinic.type,
                     min.variation = min.variation, 
                     min.obs = min.obs, epsilon = epsilon, vip = vip,
                     GeneExpression = GeneExpression, dataOmics = data.omics, omic.type = omic.type, method ='pls')
  
  # Create the results for the scale filter check
  
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
        
        if(nrow(myregulators)==0){
          ov.nr=nrow(myregulators) ## regulators nr per omic --> It could be 0
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



