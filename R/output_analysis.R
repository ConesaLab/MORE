
#########################################
#### Functions for output analysis
#########################################

## By Sonia, Monica, Maider
## 15-Oct-2023

# Function to obtain all significant pairs gene-regulator per omic --------

# For only 1 gene
GetPairs1GeneAllReg = function (gene, output) {
  
  if(output$arguments$method=='glm'){
    
    reguSignif = output$ResultsPerGene[[gene]]$relevantRegulators
    
    if (is.null(reguSignif)) {  # NO significant regulators
      return (NULL)
      
    } else {  # Significant regulators
      
      reguSignif = output$ResultsPerGene[[gene]]$allRegulators[reguSignif,]
      reguSignif = reguSignif[,c("gene", "regulator", "omic", "area", "filter")]
      return (reguSignif)
    }
    
  }
  if(output$arguments$method=='pls1' || output$arguments$method=='pls2'){
    
    reguSignif = output$ResultsPerGene[[gene]]$significantRegulators
    
    if (is.null(reguSignif)) {  # NO significant regulators
      return (NULL)
      
    } else {  # Significant regulators
      
      reguSignif = output$ResultsPerGene[[gene]]$allRegulators[reguSignif,]
      reguSignif = reguSignif[,c("gene", "regulator", "omic", "area", "filter")]
      return (reguSignif)
    }
    
  }
  
}


# For all genes
GetPairsGeneRegulator = function (genes = NULL, output) {
  
  if (is.null(genes)) genes = rownames(output$GlobalSummary$ReguPerGene)
  
  myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, output))
  
  #   colnames(myresults) = c("gene", "regulator", "omic", "area")
  return(myresults)
}

#' RegulationPerCondition
#'
#' \code{RegulationPerCondition} Function to be applied to more main function output.
#' 
#' @param output Output object of MORE main function.
#' 
#' @return Summary table containing all the relevant/significant regulators. Moreover, it provides the regression coefficient that relates the gene and the regulator for each experimental condition after testing if this coefficient is relevant/significant or not.
#'


RegulationPerCondition = function(output){
  # output: results of the getGLM/getPLS function.
  method = output$arguments$method
  
  if(method =='glm'){
    design = output$arguments$finaldesign
    Group = output$arguments$groups
    
    # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
    genes = rownames(output$GlobalSummary$ReguPerGene)
    myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, output))
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
          significatives = gsub("`", "", names(output$ResultsPerGene[[k]]$coefficients[2:nrow(output$ResultsPerGene[[k]]$coefficients), 1]))
          sign.glm = names(output$ResultsPerGene[[k]]$coefficients[2:nrow(output$ResultsPerGene[[k]]$coefficients), 1])
          
          for(i in 1:length(significatives)){
            if(any(significatives[i] == omic.representative[,2])){
              # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo GLM.
              index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
              PN = myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
            } else {
              # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]
            }
          }
          
        } else {
          # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
          myresults[myresults[,"gene"] == k, "coefficients"] = output$ResultsPerGene[[k]]$coefficients[2:nrow(output$ResultsPerGene[[k]]$coefficients), 1]
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
        significant.regulators = output$ResultsPerGene[[k]]$relevantRegulators                    # Reguladores significativos.
        model.variables = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients))[-1]       # Reguladores e interacciones en el modelo.
        
        # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
        interactions.model = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[grep(":", rownames(output$ResultsPerGene[[k]]$coefficients))])
        
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
        
        for(j in 2:nrow(output$ResultsPerGene[[k]]$coefficients)){
          regul = unlist(strsplit(gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[j]), ":"))
          
          # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
          if(any(regul %in% variables.only)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = output$ResultsPerGene[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] = output$ResultsPerGene[[k]]$coefficients[j,]
            }
          }
          
          if(any(regul %in% variables.inter)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
            } else {
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
            }
          }
          
          if(any(regul %in% variables.inter.only)){
            if(any(regul %in% significant.regulators)){
              if(length(regul) == 1){
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] + output$ResultsPerGene[[k]]$coefficients[j,]
              } else {
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
              }
            } else {
              if(length(regul) == 1){
                myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul, c(names.groups)] + output$ResultsPerGene[[k]]$coefficients[j,]
              } else {
                myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] = myresults[myresults[,"gene"] == k & myresults[,"representative"] == regul[2], regul[1]] + output$ResultsPerGene[[k]]$coefficients[j,]
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
    
  }
  if(method=='pls1' || method=='pls2'){
    design = output$arguments$finaldesign
    Group = output$arguments$groups
    
    # Creo a partir de la funcion que ya estaba hecha (linea 1657) la tabla y le anyado los huecos en blanco y cambio el nombre a "representative".
    genes = rownames(output$GlobalSummary$ReguPerGene)
    myresults = do.call("rbind", lapply(genes, GetPairs1GeneAllReg, output))
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
          significatives = gsub("`", "", names(output$ResultsPerGene[[k]]$coefficients[1:nrow(output$ResultsPerGene[[k]]$coefficients), 1]))
          sign.glm = names(output$ResultsPerGene[[k]]$coefficients[1:nrow(output$ResultsPerGene[[k]]$coefficients), 1])
          
          for(i in 1:length(significatives)){
            if(any(significatives[i] == omic.representative[,2])){
              # index.regul: para saber que regulador es el representante y asi todos los que tengan su nombre en la columna "representative" tendran su coeficiente del modelo GLM.
              index.regul = rownames(omic.representative)[which(omic.representative[,2] == significatives[i])]
              PN = myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"]                        # Sera 1 o -1, segun tenga "_P" o "_N"
              myresults[myresults[,"gene"] == k & myresults[,"representative"] == index.regul, "coefficients"] = PN*output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]                                                    # Tendra signo de la tabla si es "_P" y signo opuesto si es "_N".
            } else {
              # En caso de no pertenecer a un grupo de reguladores correlacionados, cogera su coeficiente de la tabla y lo asignara a la posicion correspondiente
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == significatives[i], "coefficients"] = output$ResultsPerGene[[k]]$coefficients[sign.glm[i], 1]
            }
          }
          
        } else {
          # Si no presenta grupo de reguladores correlacionados, simplemente sacara los coeficientes de la tabla "coefficients"
          myresults[myresults[,"gene"] == k, "coefficients"] = output$ResultsPerGene[[k]]$coefficients[1:nrow(output$ResultsPerGene[[k]]$coefficients), 1]
        }
      }
      
      
    } else {
      
      # Añado las columnas de las condiciones experimentales. Pongo "Group" porque al hacer model.matrix() siempre coloca "Group" y lo que se almacena en el objeto Group
      index = unique(Group)
      names.groups = paste("Group", index, sep = "_")
      conditions = matrix(0, nrow(myresults), length(names.groups))
      colnames(conditions) = names.groups
      rownames(conditions) = rownames(myresults)
      myresults = cbind(myresults, conditions)
      
      for(k in unique(myresults[,"gene"])){
        significant.regulators = output$ResultsPerGene[[k]]$significantRegulators                    # Reguladores significativos.
        model.variables = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients))           # Reguladores e interacciones en el modelo.
        
        # Cojo las interacciones y creo objetos que contengan los reguladores que aparecen con interaccion, solas o ambas.
        interactions.model = gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[grep(":", rownames(output$ResultsPerGene[[k]]$coefficients))])
        
        inter.variables = unlist(strsplit(interactions.model, ":", fixed = TRUE))
        if(is.null(inter.variables)){
          inter.variables = NULL                                                                            # No hay interacciones.
        } else {
          inter.variables = inter.variables[seq(2, length(inter.variables), by = 2)]                        # Reguladores que presentan interseccion con algun grupo.
        }
        
        variables.only = setdiff(setdiff(model.variables, interactions.model), inter.variables)             # Reguladores solos en el modelo, sin interacciones.
        if(length(grep("Group_", variables.only)) != 0){                                                     # No puedo hacer la interseccion con las variables significativas porque me cargo tambien omic_mc: hay que eliminar Group si esta.
          variables.only = variables.only[-grep("Group_", variables.only)]
        }
        
        variables.inter.only = intersect(inter.variables, model.variables)                                  # Reguladores con interaccion y solas.
        variables.inter = setdiff(inter.variables, model.variables)                                         # Reguladores con solo interaccion (no aparecen solas en el modelo).
        
        for(j in 1:nrow(output$ResultsPerGene[[k]]$coefficients)){
          regul = unlist(strsplit(gsub("`", "", rownames(output$ResultsPerGene[[k]]$coefficients)[j]), ":"))
          
          groups = regul[grepl(paste(names(table(Group)),collapse='|'), regul)]
          regula = setdiff(regul,groups)
          # Evaluo en que conjunto se encuentra el regulador correspondiente y segun eso asigno el coeficiente o sumo el nuevo coeficiente a lo que ya habia en esa posicion.
          if(any(regul %in% variables.only)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = output$ResultsPerGene[[k]]$coefficients[j,1]
            } 
          }
          #TO DO: regul[1] se supone que tendría que ser algo sobre los grupos y no un regulador
          if(any(regul %in% variables.inter)){
            if(any(regul %in% significant.regulators)){
              myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] + output$ResultsPerGene[[k]]$coefficients[j,1]
            } 
          }
          
          if(any(regul %in% variables.inter.only)){
            if(any(regul %in% significant.regulators)){
              if(length(regul) == 1){
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regul, c(names.groups)] + output$ResultsPerGene[[k]]$coefficients[j,1]
              } else {
                myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] = myresults[myresults[,"gene"] == k & myresults[,"regulator"] == regula, groups] + output$ResultsPerGene[[k]]$coefficients[j,1]
              }
            } 
          }
        }
      }
    }
    myresults[,6:ncol(myresults)] = signif(myresults[,6:ncol(myresults)], digits = 4) # Para que no salgan los numeros en diferentes notaciones
    myresults = myresults[,-5,drop=FALSE]
    
  }
  
  return(myresults)
}



# Plot GLM results --------------------------------------------------------

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

#' plotmore
#'
#' \code{plotmore} Graphical representation of the relationship between genes and regulators.
#' 
#' @param output Output object of MORE main function.
#' 
#' @param gene ID of the gene to be plotted.
#' @param regulator ID of the regulator to be plotted. If NULL (default), all regulators of the gene are plotted.
#' @param reguValues Vector containing the values of a regulator. If NULL (default), these values are taken from the output object as long as they are available. 
#' @param plotPerOmic If TRUE, all the relevant/significant regulators of the given gene and the same omic are plotted in the same graph. If FALSE (default), each regulator is plotted in a separate plot.
#' @param gene.col Color to plot the gene. By default, 1 (black). 
#' @param regul.col  Color to plot the regulator. If NULL (default), a color will be assigned by the function, that will be different for each regulatory omic.
#' @param order If TRUE (default), the values in X-axis are ordered.
#' @param xlab Label for the X-axis.
#' @param cont.var  Vector with length equal to the number of observations in data, which optionally may contain the values of the numerical variable (e.g. time) to be plotted on the X-axis. By default, NULL.
#' @param cond2plot Vector or factor indicating the experimental group of each value to represent. If NULL (default), the labels are taken from the experimental design matrix. 
#' 
#' @return Graphical representation of the relationship between genes and regulators.
#'


plotmore = function(output, gene, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL,...) {
  
  if(output$arguments$method=='glm'){
    
    return(plotGLM(output, gene, regulator = regulator, reguValues = reguValues, plotPerOmic = plotPerOmic,
                   gene.col = gene.col, regu.col = regu.col, order = order,
                   xlab = xlab, cont.var = cont.var, cond2plot = cond2plot,...))
  }
  
  if(output$arguments$method=='pls1'||output$arguments$method=='pls2'){
    
    return(plotPLS(output, gene, regulator = regulator, reguValues = reguValues, plotPerOmic = plotPerOmic,
                   gene.col = gene.col, regu.col = regu.col, order = order,
                   xlab = xlab, cont.var = cont.var, cond2plot = cond2plot,...))
  }
}

# order: Should the experimental groups be ordered for the plot? If TRUE, omic values are also ordered accordingly.
#        If FALSE, the function assumes they were provided in the right order for a meaningful plot.

plotGLM = function (GLMoutput, gene, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL, verbose =TRUE, ...) {
  
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
    
    if (is.null(GLMgene$relevantRegulators)) { ## No significant regulators
      
      cat("No relevant regulators were found for this gene.\n")
      
      
    } else
    {  ## Significant regulators:
      
      # Considering multicollinearity
      SigReg = GLMgene$allRegulators
      SigReg = SigReg[SigReg$Rel == 1, c("regulator", "omic", "area", "filter")]
      
      SigReg = SigReg[GLMgene$relevantRegulators,,drop = FALSE]
      
      cat(paste(nrow(SigReg), "relevant regulators are to be plotted for gene", gene)); cat("\n")
      
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
      
      return(GLMgene$allRegulators[GLMgene$relevantRegulators, -6])
    }
    
  }
  
  
  ## GENE = NULL
  
  if (is.null(gene)) {  ### Plot all genes regulated by the regulator
    
    SigniReguGene = GetPairsGeneRegulator(genes = NULL, output = GLMoutput)
    SigniReguGene = SigniReguGene[SigniReguGene[,"regulator"] == regulator,]
    myomics = SigniReguGene[,"omic"]
    myomic = unique(myomics)
    
    if (nrow(SigniReguGene) > 0) {  # When there are genes regulated by this regulator
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(GLMoutput$arguments$dataOmics[[myomic]][regulator,])
      }
      
      numGenes = length(SigniReguGene$gene)
      if(verbose) {cat(paste(numGenes, "genes are regulated by", regulator)); cat("\n")}
      
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
      
    } else { cat(paste("There are no genes relevantly regulated by", regulator)); cat("\n") }
    
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
          cat("The selected regulator was not declared as relevant by the ElasticNet\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }
        
      }
      
    }
    
  }
  
}

plotPLS = function (PLSoutput, gene, regulator = NULL, reguValues = NULL, plotPerOmic = FALSE,
                    gene.col = 1, regu.col = NULL, order = TRUE,
                    xlab = "", cont.var = NULL, cond2plot = NULL, verbose = TRUE,...) {
  
  # Colors for omics
  omic.col = colors()[c(554,89,111,512,17,586,132,428,601,568,86,390,
                        100,200,300,400,500,10,20,30,40,50,60,70,80,90,150,250,350,450,550)]
  
  if (is.null(regu.col)) {
    any.col = omic.col[1:length(PLSoutput$arguments$dataOmics)]
  } else {
    if (length(regu.col) == length(PLSoutput$arguments$dataOmics)) {
      any.col = regu.col
    } else {
      any.col = rep(regu.col, length(PLSoutput$arguments$dataOmics))
    }
  }
  names(any.col) = names(PLSoutput$arguments$dataOmics)
  
  
  # Changing margin
  par(mar = c(5,3,4,3)+0.1)
  
  # Groups to plot
  if (is.null(cond2plot)) {
    if (!is.null(PLSoutput$arguments$edesign)) {
      cond2plot = apply(PLSoutput$arguments$edesign, 1, paste, collapse = "_")
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
      myreplicates = colnames(PLSoutput$arguments$GeneExpression)
    } else {  # nothing
      myreplicates = colnames(PLSoutput$arguments$GeneExpression)
    }
  }
  
  # Cast myreplicates to character
  myreplicates = as.character(myreplicates)
  names(myreplicates) = colnames(PLSoutput$arguments$GeneExpression)
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
    
    PLSgene = PLSoutput$ResultsPerGene[[gene]]
    
    if (is.null(PLSgene)) {
      stop(paste("No GLM was obtained for gene", gene))
    }
    
    if (is.null(PLSgene$significantRegulators)) { ## No significant regulators
      
      cat("No significant regulators were found for this gene.\n")
      
      
    } else
    {  ## Significant regulators:
      
      # Considering multicollinearity
      SigReg = PLSgene$allRegulators
      SigReg = SigReg[SigReg$Sig == 1, c("regulator", "omic", "area", "filter")]
      
      SigReg = SigReg[PLSgene$significantRegulators,,drop = FALSE]
      
      cat(paste(nrow(SigReg), "significant regulators are to be plotted for gene", gene)); cat("\n")
      
      # Gene values
      geneValues = PLSgene$Y$y
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
          
          omicValues = t(PLSoutput$arguments$dataOmics[[oo]])
          omicValues = omicValues[rownames(PLSgene$X),]
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
          
          oo = PLSoutput$ResultsPerGene[[gene]]$allRegulators[rr,"omic"]
          omicValues = t(PLSoutput$arguments$dataOmics[[oo]])
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
      
      return(PLSgene$allRegulators[PLSgene$significantRegulators, -6])
    }
    
  }
  
  
  ## GENE = NULL
  
  if (is.null(gene)) {  ### Plot all genes regulated by the regulator
    
    SigniReguGene = GetPairsGeneRegulator(genes = NULL, output = PLSoutput)
    SigniReguGene = SigniReguGene[SigniReguGene[,"regulator"] == regulator,]
    myomics = SigniReguGene[,"omic"]
    myomic = unique(myomics)
    
    if (nrow(SigniReguGene) > 0) {  # When there are genes regulated by this regulator
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(PLSoutput$arguments$dataOmics[[myomic]][regulator,])
      }
      
      numGenes = length(SigniReguGene$gene)
      if(verbose) {cat(paste(numGenes, "genes are regulated by", regulator)); cat("\n")}
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        if (order) reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        lapply(1:numGenes, function (i) {
          
          geneValues = PLSoutput$ResultsPerGene[[SigniReguGene[i,"gene"]]]$Y$y
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
      } else { cat("Regulator values could not be recovered from output. Please provide them in reguValues argument to generate the plot.\n") }
      
      return(SigniReguGene$gene)
      
    } else { cat(paste("There are no genes significantly regulated by", regulator)); cat("\n") }
    
  }
  
  
  ## GENE + REGULATOR
  
  if (!is.null(gene) && !is.null(regulator)) {  ### Plot only the given gene and the given regulator
    
    geneResults = PLSoutput$ResultsPerGene[[gene]]
    
    if (is.null(geneResults)) {
      stop(paste("No GLM was obtained for gene", gene))
    } else
    {
      myomic = geneResults$allRegulators[regulator, "omic"]
      
      if (is.null(reguValues)) {  # User does not provide reguValues
        reguValues = as.numeric(PLSoutput$arguments$dataOmics[[myomic]][regulator,]) # regulator values
      }
      
      if (length(reguValues) > 0) {  # reguValues are available (recovered or given by user)
        
        if (order)  reguValues = reguValues[myorder]
        errorValuesRegu = getErrorValues(reguValues, myreplicates)
        reguValues = tapply(reguValues, myreplicates, mean, na.rm = TRUE)
        reguValues = reguValues[myrepliUni]
        names(reguValues) = myrepliUni
        
        geneValues = PLSoutput$ResultsPerGene[[gene]]$Y$y
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
        
        cat("Regulator values could not be recovered from output.\n")
        
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
          
          geneValues = PLSoutput$ResultsPerGene[[gene]]$Y$y
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
          cat("The selected regulator was not declared as significant by the ElasticNet\n")
          cat("Please, either select another regulator or provide the regulator values.\n")
        }
        
      }
      
    }
    
  }
  
}

## Summary ------------

library(ggplot2)

summary.MORE <-function(object, plot.more=FALSE){
  
  cat('A model was computed for',length(object$ResultsPerGene), 'genes.' ,'\n')
  cat(nrow(object$GlobalSummary$GenesNoregulators), 'genes had no intial regulators.' ,'\n')
  
  if(object$arguments$method == 'glm'){
    cat('For', ifelse(is.null(object$GlobalSummary$GenesNOmodel),0,object$GlobalSummary$GenesNOmodel), 'genes, the final GLM model could not be obtained.','\n')
    cat('Genes presented a mean of ',mean(object$GlobalSummary$GoodnessOfFit[,'relReg']),'relevant regulators.','\n')
    
    #Top hub genes
    relevant_regulators<-object$GlobalSummary$ReguPerGene[,c(grep('-Rel$',colnames(object$GlobalSummary$ReguPerGene)))]
    #globally
    
    s_rel_reg<-apply(relevant_regulators, 1, sum)
    cat('These are the top 10 hub genes and the number of relevant regulators for each:\n')
    print(s_rel_reg[rev(tail(order(s_rel_reg),10))])
    
    cat('The top 10 hub genes by omics:\n')
    
    for (i in 1:ncol(relevant_regulators)){
      cat('These are the top 10 hub genes for ',gsub('-Rel$','',colnames(relevant_regulators)[i]) ,'and the number of relevant regulators for each in this omic:\n')
      print(relevant_regulators[rev(tail(order(relevant_regulators[,i]),10)),i])
    }
    
    #Master regulators
    
    m_rel_reg<-lapply(object$ResultsPerGene, function(x) x$relevantRegulators)
    m_rel_reg <- unlist(m_rel_reg)
    
    ## Count occurrences
    mrel_vector <- table(m_rel_reg)
    cat('These are the top 10 master regulators and the number of genes that they regulate:\n')
    mreg<-mrel_vector[rev(tail(order(mrel_vector),10))]
    for (i in 1:10) {
      cat(names(mreg)[i], 'relevantly regulates', mreg[i][[1]], 'genes. \n')
      if(plot.more){
        par(mfrow=c(2,4))
        plotGLM(object, gene = NULL, regulator = names(msig)[i], plotPerOmic = FALSE ,order = FALSE, gene.col = 'skyblue', regu.col = 'tan1', verbose = FALSE)
      }
    }
    
  }
  else{
    cat('Genes presented a mean of ',mean(object$GlobalSummary$GoodnessOfFit[,'sigReg']),'significant regulators.','\n')
    
    #Top hub genes
    significant_regulators<-object$GlobalSummary$ReguPerGene[,c(grep('-Sig$',colnames(object$GlobalSummary$ReguPerGene)))]
    #globally
    
    s_sig_reg<-apply(significant_regulators, 1, sum)
    cat('These are the top 10 hub genes and the number of significant regulators for each:\n')
    print(s_sig_reg[tail(order(s_sig_reg),10)])
    
    cat('The top 10 hub genes by omics:\n')
    
    for (i in 1:ncol(significant_regulators)){
      cat('These are the top 10 hub genes for ',gsub('-Rel$','',colnames(significant_regulators)[i]) ,'and the number of relevant regulators for each in this omic:\n')
      print(significant_regulators[rev(tail(order(significant_regulators[,i]),10)),i])
    }
    
    #Master regulators
    
    m_sig_reg<-lapply(object$ResultsPerGene, function(x) x$significantRegulators)
    m_sig_reg <- unlist(m_sig_reg)
    
    ## Count occurrences
    msig_vector <- table(m_sig_reg)
    cat('These are the top 10 master regulators and the genes that they significantly regulate:\n')
    msig<-msig_vector[rev(tail(order(msig_vector),10))]
    for (i in 1:10) {
      cat(names(msig)[i], 'significantly regulates', msig[i][[1]], 'genes. \n')
      if(plot.more){
        par(mfrow=c(2,4))
        plotPLS(object, gene = NULL, regulator = names(msig)[i], plotPerOmic = FALSE ,order = FALSE, gene.col = 'skyblue', regu.col = 'tan1', verbose = FALSE)
        
      }
    }
    
  }
  
}


summary_plot<-function(output, output_regpcond, by_genes =TRUE){
  
  #output should by a MORE object
  #output_regpcond should by the output of calculating RegulationPerCondition to a MORE object
  #by_genes by default TRUE calculates the percentage of genes with significant regulators globally and per omic. FALSE to calculate the percentage of significant regulations globally and per omic.
  
  if(by_genes){
    #Calculate the vector with % of significant regulators by condition and globally
    
    ngroups = length(unique(output$arguments$groups))
    omics = names(output$arguments$dataOmics)
    totalgenes = length(output$ResultsPerGene)
    
    #Create all the counts needed globally and per groups
    
    cts = matrix(NA, nrow=(ngroups)+1,ncol=length(omics)+1)
    
    for (i in 1:((ngroups)+1)){
      #Create the global values
      for (j in 1:((length(omics))+1)){
        if (i==1 && j==1){
          cts[i,j]= length(unique(output_regpcond$gene))
        }else if(i==1){
          cts[i,j]= length(unique(output_regpcond[output_regpcond$omic==omics[j-1],]$gene))
        }else if(j==1){
          cts[i,j] = length(unique(output_regpcond[output_regpcond[,5+(i-1)]!=0,]$gene))
        }else{
          cts[i,j] = length(unique(output_regpcond[intersect(which(output_regpcond[,5+(i-1)]!=0),which(output_regpcond$omic==omics[j-1])),]$gene))
        }
        
      }
    }
    
    #Convert to percentage
    
    cts<-cts/totalgenes
    
    #Create a df with the percentage of genes with significant regulators by omic and condition
    df <- data.frame(Group=rep(c('Global',unique(output$arguments$groups)), times=length(omics) +1),
                     omic=rep(c('Any',names(output$arguments$dataOmics)),each = ngroups+1),
                     genes=as.vector(cts))
    
    
    ggplto2::ggplot(data=df, aes(x=omic, y=genes, fill=Group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme_minimal()+
      scale_fill_conesa(palette = "complete")+  
      labs(x="Omic", y = "% genes with significant regulators") +
      theme(legend.text = element_text(size = 12)) 
    
  } else{
    #Calculate the vector with % of significant regulations by condition in each omic
    
    ngroups = length(unique(output$arguments$groups))
    omics = names(output$arguments$dataOmics)
    
    #Create all the counts needed globally and per groups
    
    cts = matrix(NA, nrow=(ngroups),ncol=length(omics))
    
    for (i in 1:ngroups){
      #Create the global values
      for (j in 1:length(omics)){
        cts[i,j] = length(output_regpcond[intersect(which(output_regpcond[,5+i]!=0),which(output_regpcond$omic==omics[j])),]$regulator)/ nrow(output$arguments$associations[[omics[j]]])
        
      }
    }
    
    
    #Create a df with the percentage of genes with significant regulators by omic and condition
    df <- data.frame(Group=rep(unique(output$arguments$groups), times=length(omics)),
                     omic=rep(names(output$arguments$dataOmics),each = ngroups),
                     genes=as.vector(cts))
    
    
    ggplot2::ggplot(data=df, aes(x=omic, y=genes, fill=Group)) +
      geom_bar(stat="identity", position=position_dodge()) +
      theme_minimal()+
      scale_fill_conesa(palette = "complete")+
      labs(x="Omic", y = "% significant regulations") +
      theme(legend.text = element_text(size = 12)) 
    
  }
  
}