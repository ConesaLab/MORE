#########################################
#### Auxiliar functions
#########################################

## By Sonia
## 6-Jun-2017



# Generating design matrix with interactions between experimental variables ------------------------------
## We only consider interactions of order 2

# GenerateDesignMatrix = function (interactions.exp, degree, edesign, cont.var) {
# 
#   if ((cont.var %in% colnames(edesign)) && (degree > 1)) {  # polynomial terms
#     politerms = paste0(cont.var, 2:degree)
#     poliEquat = NULL
#     for (i in 2:degree) {
#       poliEquat = cbind( poliEquat, edesign[,cont.var]^i)
#     }
#     colnames(poliEquat) = politerms
# 
#   } else { politerms = NULL }
# 
# 
#   if (interactions.exp) {   # we compute all possible interactions between exp variables
# 
#     if (NCOL(edesign) > 1) {
#       pares = combn(colnames(edesign), 2) # all possible pairs
#       pares = apply(pares, 2, function (x) paste(sprintf("`%s`", x), collapse="*")) # interactions of order 2
#       fff = paste0("~ ", paste(pares, collapse = "+")) # complete formula
#       fff = as.formula(fff)
#       desmat = model.matrix(fff, edesign)[,-1]  # design matrix with all predictors and all interactions of order 2
# 
#       # polynomial terms?
#       if (!is.null(politerms)) {
#         expcond = setdiff(colnames(desmat), cont.var) # remove cont.var
#         # remove interactions to avoid having int > 2
#         if (length(grep(":", expcond, fixed = TRUE)) > 0) expcond = expcond[-grep(":", expcond, fixed = TRUE)]
#         desmat = as.data.frame(cbind(desmat, poliEquat))
#         fff = paste0("~ ", paste(sapply(politerms, function (x) paste(expcond, x, sep = ":")), collapse = "+"))
#         fff = as.formula(fff)
#         desmat = as.data.frame(cbind(desmat, model.matrix(fff, desmat)[,-1, drop = FALSE]))
#       }
# 
# 
#     } else {  # only 1 experimental covariate
#       desmat = model.matrix(~., data = edesign)[, -1, drop = FALSE]
#       if (!is.null(politerms)) {
#         expcond = setdiff(colnames(desmat), cont.var) # remove cont.var
#         desmat = as.data.frame(cbind(desmat, poliEquat))
#         if (length(expcond) > 0) {
#           fff = paste0("~ ", paste(sapply(politerms, function (x) paste(expcond, x, sep = ":")), collapse = "+"))
#           fff = as.formula(fff)
#           desmat = as.data.frame(cbind(desmat, model.matrix(fff, desmat)[,-1]))
#         }
#       }
#     }
# 
#   } else {  # no interactions
#     desmat = model.matrix(~., data = edesign)[, -1, drop = FALSE]
#     if (!is.null(politerms)) desmat = as.data.frame(cbind(desmat, poliEquat))
#   }
# 
#   return(desmat)
# 
# }
# 


# Removing regulators with low variation ----------------------------------

LowVariationRegu = function(min.variation, data.omics, ExpGroups, associations, Allgenes, omic.type) {

  if (all(is.na(min.variation))){
    
    for (ov in names(associations)){
      if(!is.null(associations[[ov]])){
        myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2] # removing regulators not associated to our genes
        data.omics[[ov]]=data.omics[[ov]][intersect(myreg, rownames(data.omics[[ov]])),] ## Reduced data.omics
        rm("myreg")
      }
    }
    #### Low variation cutoff is computed automatically
    data.omicsMean = vector("list", length=length(data.omics))
    for(i in 1:length(data.omics)){
      data.omicsMean[[i]]=t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
    }
    names(data.omicsMean)=names(data.omics)

    percVar = c(10, 0.1)
    names(percVar) = 0:1
    percVar = percVar[as.character(omic.type)]
    names(percVar)=names(data.omicsMean)

    # Applying Low Variation filter
    LowVar=LowVariatFilter(data=data.omicsMean, method="sd", percVar=percVar, omic.type = omic.type)

    # data.omicsMean reduced: without NA and LV
    data.omicsMean=LowVar$data
    
    ## data.omics reduced: only mygenes and without NA and LV
    for (ov in names(data.omics)){
      data.omics[[ov]] = data.omics[[ov]][rownames(data.omicsMean[[ov]]),]
      # Remove regulators from associations that have been removed due to LowVariation
      associations[[ov]] = associations[[ov]][associations[[ov]][,2] %in% rownames(data.omics[[ov]]),,drop = FALSE]
      
    }

  }
  else {

    #### Low variation cutoff is set by the user

    # removing regulators not associated to our genes only when there is associations matrix
    for (ov in names(associations)){
      if(!is.null(associations[[ov]])){
        myreg=associations[[ov]][associations[[ov]][,1] %in% Allgenes,2]
        data.omics[[ov]]=data.omics[[ov]][intersect(myreg, rownames(data.omics[[ov]])),]
        rm("myreg")
      }
    }
    
    # Creating vector for min.variation
    if (length(min.variation) == 1) {  ## Including min.variation = 0. I need a vector with omics names
      min.variation=rep(min.variation,length(data.omics))
      names(min.variation)=names(data.omics)
    } 

    # computing mean per condition in data.omics
    data.omicsMean=vector("list", length = length(data.omics))
    for(i in 1:length(data.omics)){
      data.omicsMean[[i]] = t(apply(data.omics[[i]], 1, tapply, ExpGroups, mean))
    }
    names(data.omicsMean) = names(data.omics)

    # Applying Low Variation filter
    LowVar=LowVariatFilter(data = data.omicsMean, method = "user", percVar = min.variation, omic.type = omic.type)

    data.omicsMean=LowVar$data  ## data.omicsMean reduced

    ## data.omics reduced: only mygenes and without NA and LV
    for (ov in names(data.omics)){
      data.omics[[ov]] = data.omics[[ov]][rownames(data.omicsMean[[ov]]),]
      # Remove regulators from associations that have been removed due to LowVariation
      associations[[ov]] = associations[[ov]][associations[[ov]][,2] %in% rownames(data.omics[[ov]]),,drop = FALSE]
      
    }

  }

  rm("data.omicsMean"); gc()

  # Regulators removed due to low variation filter
  myregLV=LowVar$LV.reg
  rm("LowVar"); gc()

  cat("Number of regulators with low variation:\n")
  print(sapply(myregLV, length))
  cat("\n")

  return(list("myregLV" = myregLV, "data.omics" = data.omics, "associations" = associations))
}





# Low variation filtering -------------------------------------------------
## data: For computing variability. It could be mean values between replicates
## method: One of "sd","range", "IQrange" or "user"
## percVar: percentage of variation defined by the user

LowVariatFilter=function(data, method, percVar, omic.type){
  
  SummaryRes = LV.reg = vector("list", length=length(data))
  names(SummaryRes) = names(LV.reg) = names(data)
  
  for (ov in names(data)) {
    if(is.na(percVar[ov])){
      method.low="sd"
      if(omic.type[ov]==0){
        percVar[ov]=10
      }else{
        percVar[ov]==0.1
      }
    }else{
      method.low =method
    }
    if (omic.type[ov] == 0) {  # numerical regulators
      
      if (method.low=="sd") {
        met=apply(data[[ov]], 1, sd, na.rm=TRUE)  ## Compute standard deviation between conditions
        maxMet=max(met)*(percVar[ov]/100)  ## Compute minimum variation allowed
        myreg=met[met>maxMet]  # Regulators to be kept
        LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
        data[[ov]]=data[[ov]][names(myreg), ,drop=FALSE]
      }
      
      if(method.low=="user") {
        if (min(dim(data[[ov]])) > 0) {
          met = apply(data[[ov]], 1, function(x) max(x, na.rm=TRUE)-min(x, na.rm=TRUE) )
          maxMet = percVar[ov] ## We don't consider the max, just the percentage defined by the user
          myreg=met[met>maxMet]
          LV.reg[[ov]]=names(met[met<=maxMet]) ## Keep names of removed regulators
          data[[ov]]=data[[ov]][names(myreg), , drop=FALSE]
        }
      }
    }
    
    if (omic.type[ov] == 1) {  # binary categorical regulators
      
      if(method.low=='sd'){
        met = apply(data[[ov]], 1, function (x) { max(x, na.rm = TRUE)-min(x, na.rm = TRUE) })  ## Compute maximun variation between groups
        maxMet=max(met)/10 ## Compute the minimum variation allowed
        myreg = met[met > maxMet]  # Regulators to be kept
        LV.reg[[ov]] = names(met[met <= maxMet]) ## Keep names of removed regulators
        data[[ov]] = data[[ov]][names(myreg), ,drop=FALSE]
      }
      
      if(method.low == 'user'){ 
        met = apply(data[[ov]], 1, function (x) { max(x, na.rm = TRUE)-min(x, na.rm = TRUE) })  ## Compute maximun variation between groups
        myreg = met[met > percVar[ov]]  # Regulators to be kept
        LV.reg[[ov]] = names(met[met <= percVar[ov]]) ## Keep names of removed regulators
        data[[ov]] = data[[ov]][names(myreg), ,drop=FALSE]
      }
    }
  }
  
  results=vector("list", length=2)
  results[[1]]=data
  results[[2]]=LV.reg
  names(results)=c("data", "LV.reg")
  
  return(results)
}


# Adding interations with regulators -------------------------------------

RegulatorsInteractions = function (interactions.reg, reguValues, des.mat, cont.var, GeneExpression, gene) {

  # Adding regulators
  if (is.null(des.mat)) {
    des.mat2 = data.frame(reguValues, check.names = FALSE)
    des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
    colnames(des.mat2)[1] = "response"
  } else {
    des.mat2 = data.frame(des.mat, reguValues, check.names = FALSE)
    ### WITH INTERACTIONS with regulators
    if (interactions.reg > 0) {
      expcond = colnames(des.mat)

      if (interactions.reg == 1) {  # Max order of interaction = 2 & Interactions with cont.var not allowed
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
        if (!is.null(cont.var)) {
          if (length(grep(cont.var, expcond)))  expcond = expcond[-grep(cont.var, expcond)]
        }
      }

      if (interactions.reg == 2) { # Max order of interaction = 2
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      }

      fff = paste0("~ ",
                   paste(sapply(colnames(reguValues),
                                function (x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                         collapse = "+"))

      fff = as.formula(fff)
      inter.var = model.matrix(fff, des.mat2)[,-1, drop = FALSE]
      colnames(inter.var) = strsplit(as.character(fff),"\\s*\\+\\s*")[-1][[1]]
      des.mat2 = cbind(des.mat2, inter.var)

      des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
      colnames(des.mat2)[1] = "response"

      colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))


      ## PUEDE OCURRIR QUE AL METER LAS INTERACCIONES TODA LA COLUMNA SEA 0. HAGO UN FILTRO PREVIO
      sd.regulators = apply(des.mat2[,-1, drop = FALSE], 2, sd, na.rm=TRUE)
      regulators0 = names(sd.regulators[sd.regulators==0])
      if (length(regulators0)>0) des.mat2 = des.mat2[, setdiff(colnames(des.mat2), regulators0), drop=FALSE]


    } else  {    ### WITHOUT INTERACTIONS

      des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
      colnames(des.mat2)[1] = "response"

    }
  }

  return(des.mat2)

}


# ElasticNet variable selection -------------------------------------------

ElasticNet = function (family2, des.mat2, epsilon, elasticnet) {

  if (!is.null(elasticnet) && (ncol(des.mat2) > 2)) {  # Variable selection with ElasticNet

    if (!is.null(family2)) {

      if (nrow(des.mat2) > 200) { mynfolds =  ceiling(nrow(des.mat2)/20) } else { mynfolds = nrow(des.mat2) }  # for CV

      cvEN = cv.glmnet(x = as.matrix(des.mat2[,-1]),
                       y = des.mat2[,1], 
                       nfolds = mynfolds, alpha = elasticnet, standardize = FALSE, thres = epsilon,
                       family = family2, grouped = FALSE)

      myS = cvEN$lambda.min  # optimum penalization parameter
      y.fitted = predict(cvEN, s = myS, newx = as.matrix(des.mat2[,-1]))

      if (length(myS) > 1) {  # more than 1 values for lambda
        myCVerror = cvEN$cvm[sapply(myS, function (i) which(cvEN$lambda == i))]
        myS = myS[which.min(myCVerror)]
      }

      if (length(myS) == 1) {
        x = des.mat2[,-1]
        sel = colnames(x)[which(coef(cvEN, s = myS)[-1,1] != 0)] # selected coefficients without intercept
        selcoef = as.data.frame(as.matrix(coef(cvEN, s = myS)[which(coef(cvEN, s = myS)[,1] != 0),,drop=FALSE]))
        mycoef = c(colnames(des.mat2)[1], sel)
        
        removedCoefs = setdiff(colnames(des.mat2), mycoef)
        des.mat2 = des.mat2[, mycoef, drop = FALSE]
      }
      
      m = modelcharac(cvEN, myS, des.mat2[,1], y.fitted)
      isModel = TRUE

    } else { removedCoefs = NULL; selcoef=NULL; y.fitted = NULL;m = list('AICc'=NULL,'R.squared'=NULL,'cvRMSE'=NULL ); isModel = NULL }

  } else { removedCoefs = NULL; selcoef=NULL;y.fitted = NULL; m = list('AICc'=NULL,'R.squared'=NULL,'cvRMSE'=NULL ); isModel = NULL }

  return(list("des.mat2" = des.mat2, "removedCoefs" = removedCoefs, "coefficients" = selcoef, 'fitted.values' = y.fitted, 'm' = m, 'isModel'=isModel))
}

modelcharac = function(fitted.glm,s, y, y.fitted){
  n = fitted.glm$glmnet.fit$nobs
  
  R2 = fitted.glm$glmnet.fit$dev.ratio[which(fitted.glm$glmnet.fit$lambda == s)]
  RMSE = sqrt(sum((y-y.fitted)^2)/n)
  cvRMSE= sqrt(sum((y-y.fitted)^2)/n)/mean(y)
  
  return(list('R.squared'=R2, 'RMSE' = RMSE ,'cvRMSE'=cvRMSE))
}

# Functions needed for PLS


# Adding interations with regulators -------------------------------------

RegulatorsInteractionsPLS = function (interactions.reg, reguValues, des.mat, clinic.type, cont.var, GeneExpression, gene, regu, omic.type) {
  
  # Adding regulators
  if (is.null(des.mat)) {
    des.mat2 = data.frame(reguValues, check.names = FALSE)
    des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
    colnames(des.mat2)[1] = "response"
  } else {
    res.mat = NULL
    for (i in 1:ncol(des.mat)) {
      res.mat = cbind(res.mat,as.matrix(dummy_cols(des.mat[,i, drop = FALSE])[,-1]))
    }
    colnames(res.mat) = sub('.*_','Group_',colnames(res.mat))
    rownames(res.mat) = rownames(des.mat)
    des.mat = res.mat
    
    res.mat = NULL
    for (ov in names(regu)){
      if (ov =='clinic'){
        des.mat2 = as.data.frame(reguValues[, regu[[ov]],drop =FALSE])
        for (i in 1:length(clinic.type)) {
          if(clinic.type[i]==1){
            res.mat = cbind(res.mat,as.matrix(dummy_cols(des.mat2[,i,drop=FALSE])[,-1, drop=FALSE]))
          }else{
            res.mat = cbind(res.mat,as.matrix(des.mat2[,i,drop=FALSE]))
          }
        }
      }else{
        des.mat2 = as.data.frame(reguValues[, regu[[ov]],drop =FALSE])
        if (ncol(des.mat2!=0)){
          if (omic.type[ov] == 1){
            des.mat2 = apply(des.mat2, 2, factor)
            res.mat = cbind(res.mat,as.matrix(dummy_cols(des.mat2)[,-c(1:ncol(des.mat2)), drop=FALSE]))
          } else{
            res.mat = cbind(res.mat,as.matrix(des.mat2))
          } 
        }
      }
    }
    
    rownames(res.mat) = rownames(reguValues)
    reguValues = res.mat
    rm(res.mat);rm(des.mat2);gc();
    des.mat2 = data.frame(des.mat, reguValues, check.names = FALSE)
    ### WITH INTERACTIONS with regulators
    if (interactions.reg > 0) {
      expcond = colnames(des.mat)
      
      if (interactions.reg == 1) {  # Max order of interaction = 2 & Interactions with cont.var not allowed
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
        if (!is.null(cont.var)) {
          if (length(grep(cont.var, expcond)))  expcond = expcond[-grep(cont.var, expcond)]
        }
      }
      
      if (interactions.reg == 2) { # Max order of interaction = 2
        if (length(grep(":", expcond, fixed = TRUE)))  expcond = expcond[-grep(":", expcond, fixed = TRUE)]
      }
      
      fff = paste0("~ ",
                   paste(sapply(colnames(reguValues),
                                function (x) paste(expcond, sprintf("`%s`", x), sep = ":")),
                         collapse = "+"))
      
      fff = as.formula(fff)
      inter.var = model.matrix(fff, des.mat2)[,-1, drop = FALSE]
      colnames(inter.var) = strsplit(as.character(fff),"\\s*\\+\\s*")[-1][[1]]
      des.mat2 = cbind(des.mat2, inter.var)
      
      des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
      colnames(des.mat2)[1] = "response"
      
      colnames(des.mat2) = gsub("\`", "", colnames(des.mat2))
      
      
      ## PUEDE OCURRIR QUE AL METER LAS INTERACCIONES TODA LA COLUMNA SEA 0. HAGO UN FILTRO PREVIO
      sd.regulators = apply(des.mat2[,-1, drop = FALSE], 2, sd, na.rm=TRUE)
      regulators0 = names(sd.regulators[sd.regulators==0])
      if (length(regulators0)>0) des.mat2 = des.mat2[, setdiff(colnames(des.mat2), regulators0), drop=FALSE]
      
      
    } else  {    ### WITHOUT INTERACTIONS
      
      des.mat2 = cbind(t(GeneExpression[gene,]), des.mat2)
      colnames(des.mat2)[1] = "response"
      
    }
  }
  
  return(des.mat2)
  
}

# Scale for PLS models -------------------------------

ScalePLSdesmat = function(des.mat2, scaletype, center, scale){
  
  res.mat = NULL
  for (i in 1:ncol(des.mat2)) {
    res.mat = cbind(res.mat,as.matrix(dummy_cols(des.mat2[i])[,-1]))
  }
  #Change the name to avoid conflicts with RegulationPerCondition
  colnames(res.mat) = sub('.*_','Group_',colnames(res.mat))
  if(scaletype=='auto'){
    #Autoescalado: todas las variables tienen el mismo peso en el modelo, indiferentemente del modelo
    res.mat = scale(res.mat, center = center, scale = scale)
  }
  if(scaletype=='pareto'){
    #Escalado tipo pareto
    res.mat = scale(res.mat, center = center, scale = scale) / (ncol(des.mat2)^(1/4))
  }
  if(scaletype=='block'){
    #Cada bloque de variables tiene el mismo peso total en el modelo
    res.mat = scale(res.mat, center = center, scale = scale) / (ncol(des.mat2)^(1/2))
  }
  rownames(res.mat) = rownames(des.mat2)
  return(res.mat)
}

ScalePLS= function(reguVal, regu, omic.type, scaletype, center, scale){
  res.mat = NULL
  for (ov in names(regu)){
    des.mat2 = as.data.frame(reguVal[, regu[[ov]],drop =FALSE])
    
    if(scaletype=='auto'){
      #Autoescalado: todas las variables tienen el mismo peso en el modelo, indiferentemente del modelo
      des.mat2 = scale(des.mat2, center = center, scale = scale)
    }
    if(scaletype=='pareto'){
      #Escalado tipo pareto
      des.mat2 = scale(des.mat2, center = center, scale = scale) / (ncol(des.mat2)^(1/4))
    }
    if(scaletype=='block'){
      #Cada bloque de variables tiene el mismo peso total en el modelo
      des.mat2 = scale(des.mat2, center = center, scale = scale) / (ncol(des.mat2)^(1/2))
    }
    
    res.mat = cbind(res.mat,des.mat2)
    
  }
  
  rownames(res.mat) = rownames(reguVal)
  
  return(res.mat)
}

# p-values for PLS models ------------------------------------------------

p.coef<-function(pls,R, datospls){
  #pls: modelo PLS generado con la libreria ropls
  #R: número de veces a repetir la prueba
  #datospls: matriz de datos utilizada para crear el modelo
  k=pls@summaryDF$pre
  coefmod<-pls@coefficientMN
  Y=datospls[,1]
  a<-NULL
  for (i in 1:R){
    Yperm=sample(Y, replace=FALSE)
    plsda.opls<-opls(datospls[,-1], scale(Yperm), scaleC='none', predI=k,
                     info.txtC='none', fig.pdfC='none', crossvalI=1,permI = 0)
    a<-cbind(a,plsda.opls@coefficientMN)
  }
  p.coefs<-matrix(2, nrow=nrow(a), ncol=1)
  for(i in 1:nrow(a)){
    coefs<-a[i,]
    pvalor = ifelse(coefmod[i] > 0, (sum(coefs > coefmod[i]) + sum(coefs < -coefmod[i]) )/ R, (sum(coefs < coefmod[i]) + sum(coefs > -coefmod[i])) / R)
    p.coefs[i,1]<-pvalor
  }
  rownames(p.coefs)<-rownames(coefmod)
  return(p.coefs)
}


p.valuejack<-function(pls, datospls,alfa){
  
  #pls: modelo PLS generado con la libreria ropls
  #datospls: matriz de datos utilizada para crear el modelo
  
  #Caso de Leave-one-out
  k=pls@summaryDF$pre
  coefmod=pls@coefficientMN
  a=NULL
  
  for (i in 1: nrow(datospls)) {
    
    pls.opls=suppressWarnings(opls(datospls[-i,-1], scale(datospls[-i,1]), scaleC='none', predI=k,
                                   info.txtC='none', fig.pdfC='none', crossvalI=1, permI = 0))
    
    if(length(pls.opls)<length(coefmod)){
      exclvar=setdiff(rownames(coefmod),rownames(pls.opls@coefficientMN))
      b<-matrix(0,ncol=1,nrow = length(exclvar))
      rownames(b)<-exclvar
      pls.opls@coefficientMN=rbind(pls.opls@coefficientMN,b)
      
    }
    order <- match(rownames(coefmod), rownames(pls.opls@coefficientMN))
    plscoefficientMN <- pls.opls@coefficientMN[order, ]
    a=cbind(a,plscoefficientMN)
  }
  
  est = (a-matrix(rep(coefmod, ncol(a)), ncol = ncol(a)))^2
  estjack = sqrt(((nrow(datospls)-1)/nrow(datospls)) * rowSums(est))
  
  pvalor = 2*pt(abs(coefmod/estjack), df = nrow(datospls)-1, lower.tail = FALSE)
  #calculo de los intervalos de confianza para ver que todo va bien, aún y todo podemos ahorrar esto no es necesario en nuestro caso
  #res =cbind(pvalor, paste('(',coefmod- qt(1 - (alfa)/ 2, df = nrow(datospls)-1)* estjack,'-', coefmod+ qt(1 - (alfa)/ 2, df = nrow(datospls)-1)* estjack,')' ))
  #colnames(res) = c('pvalue', 'CI')
  colnames(pvalor) = c('pvalue')
  return(pvalor)
}
