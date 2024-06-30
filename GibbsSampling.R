
#####################  Función auxiliar para la previa SBSB  ################ 

Gdelta <- function (r, variables, groups) {
  
  p1 <- length(groups) # número de grupos originales
  p2 <- length(groups[as.logical(variables)[1:p1]]) # número de grupos activos
  grouped <- as.vector(lengths(groups)[as.logical(variables)[1:p1]]) 
  # num de variables agrupadas originales de los grupos activos
  
  indices <- 1:(min(c(grouped[1], r))) 
  if (p2 > 1) {
    for (i in 2:p2) {
      # creo matriz con todos los índices posibles, de 1 a l_j para cada j
      indices <- merge(indices, 1:(min(c(grouped[i], r))), by = NULL)
    }
  } else {indices <- matrix(indices, ncol = 1)} # si no hay grupos activos
  
  indices_to_remove <- c()
  for (i in 1:nrow(indices)) {
    if (sum(indices[i,]) != r) { # me quedo con los que son solución
      indices_to_remove <- c(indices_to_remove, i) # qué filas quito
    }
  }
  if (!is.null(indices_to_remove)) {
    indices <- as.matrix(indices[-indices_to_remove, ])
  }
  
  suma <- 0 # sumatorio G
  for (i in 1:nrow(indices)) {
    suma <- suma + prod(choose(grouped, indices[i,]))
  }
  return(suma)
}

###########  Funciones para calcular las previas de los modelos  #############

SBSB <- function (gamma, positions, groups) {
  
  p1 <- length(groups) # número de grupos originales
  k1 <- nrow(positions) - p1 # número de vars sing originales
  
  variables <- positions %*% gamma # qué variables tiene gamma
  
  m2 <- sum(as.logical(variables[1:p1])) # num de grupos activos
  if (k1 > 0) {
    # número de variables singulares activas
    num_vars_sing <- sum(variables[(p1 + 1):nrow(positions),])
  } else { num_vars_sing <- 0 }
  
  # numero de variables agrupadas originales de cada grupo activo
  num_grouped <- sum(lengths(groups)[as.logical(variables[1:p1])])
  
  r <- sum(variables) # dimensión del modelo
  SB1 <- 1/((k1 + p1 + 1) * choose(k1 + p1, m2 + num_vars_sing))
  if (m2 == 0) {SB2 <- 1} else{
    SB2 <- 1/(Gdelta(r - num_vars_sing, variables, groups) * 
                (num_grouped - m2 + 1))
  }
  return(SB1 * SB2)
}


ConstConst <- function (gamma, positions, groups) {
  
  p1 <- length(groups) # numero de grupos originales
  k1 <- nrow(positions) - p1 # num de vars singulares originales
  
  variables <- positions %*% gamma # variables en cada grupo(fila)
  
  num_models <- 1
  for (j in 1:p1) {
    if (variables[j] > 0) {
      grouped <- as.numeric(lengths(groups)[j])
      num_models <- num_models * (2^grouped - 1)
      # calcula el número de modelos teniendo en cuenta todas las vars 
      # originales de los grupos que están activos en gamma
    }
  }
  
  Const1 <- 1/(2^(k1 + p1))
  Const2 <- 1/num_models
  return(Const1 * Const2)
}


SBConst <- function (gamma, positions, groups) {
  
  p1 <- length(groups) # numero de grupos originales
  k1 <- nrow(positions) - p1 # num de vars singulares originales
  
  variables <- positions %*% gamma
  
  num_models <- 1
  for (j in 1:p1) {
    if (variables[j] > 0) {
      grouped <- as.numeric(lengths(groups)[j])
      num_models <- num_models * (2^grouped - 1)
      # calcula el número de modelos teniendo en cuenta todas las vars
      # originales de los grupos que están activos en gamma
    } 
  }
  
  m2 <- sum(as.logical(variables[1:p1])) # num de grupos activos
  if (k1 > 0) {num_vars_sing <- sum(variables[(p1 + 1):nrow(positions),])
  } else {
    num_vars_sing <- 0 # si no hay variables singulares
  }
  
  SB1 <- 1/((k1 + p1 + 1) * choose(k1 + p1, m2 + num_vars_sing))
  Const2 <- 1/num_models
  return(SB1 * Const2)
}


auxSB <- function (positions, groups) {
  
  previasSB <- list()
  p1 <- length(groups) # número de grupos originales
  k1 <- nrow(positions) - p1 # número de vars sing originales
  grouped <- lengths(groups) # numero de variables originales de cada grupo
  
  rf <- sum(rep(1, sum(positions))) # num de variables máximo
  for(i in 1:rf) {
    matriz <- 0:min(c(k1, i))
    for (j in 1:p1) {
      matriz <- merge(matriz, 0:(min(c(grouped[j], i))), by = NULL)
    }
    
    indices_to_remove <- c()
    for (j in 1:nrow(matriz)) {
      if (sum(matriz[j,]) != i) {
        indices_to_remove <- c(indices_to_remove, j) 
        # quito las filas que no lo verifican
      }
    }
    if (!is.null(indices_to_remove)) {
      indices <- as.matrix(matriz[-indices_to_remove, ])
    } else {
      indices <- as.matrix(matriz)
    }
    
    suma <- 0
    for (j in 1:nrow(indices)) {
      suma <- suma + choose(k1, indices[j, 1]) * 
        prod(choose(grouped, indices[j, -1]))
    }
    
    previasSB[[i]] <- 1/(suma * (rf + 1))
  }
  return(previasSB)
}

#########  Función para escribir la fórmula del modelo gamma  #########

formula.gamma <- function (gamma, formula_original, orig.names) {
  
  formula_gamma <- paste0(formula_original[2], " ~ 1 ") # intercepto
  
  for(i in 1:length(gamma)) {
    if(gamma[i] == 1) {
      formula_gamma <- paste0(formula_gamma, " + ", orig.names[i])
    }
  }
  return(as.formula(formula_gamma)) 
}

### Función para calcular el factor bayes de cada modelo respecto al nulo  ###

BayesFactor <- function (formula_gamma, data, log_marginal_nulo) {
  
  bas_object <- bas.glm(formula_gamma,
                        data = data, 
                        n.models = 1, # solo el que queremos
                        method = "deterministic",
                        betaprior = robust(), 
                        family = binomial(),
                        modelprior = beta.binomial(1,1))
  
  # BF = e^logmar(i)/e^logmar(0) = e^(logmar(i) - logmar(0))
  BF <- exp(bas_object$logmarg - log_marginal_nulo)
  return(BF)
}

###############          CÓDIGO GIBBS SAMPLING           ##################


GS <- function (formula, 
                data,
                null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                groups = NULL,
                prior.betas = "Robust",
                prior.models = "SBSB",
                n.iter = 10000,
                init.model = "Full",
                n.burnin = 1000,
                n.thin = 10, 
                seed = runif(1, 0, 16091956)) 
{
  # Calculamos el tiempo de computación
  t <- proc.time()
  
  # Comprobaciones previas
  if(!requireNamespace("BAS", quietly = TRUE)) {
    install.packages("BAS")
  }
  library(BAS)
  
  formula <- as.formula(formula)
  null.model <- as.formula(null.model)
  
  # groups debe ser una lista
  if (!is.list(groups)) {
    stop("Argument groups must be a list.\n")
  }
  warning("Be sure that the names in the groups list coincide with some of the 
  explanatory vars\n and these cannot be in the null model\n")
  
  # Respuesta en el nulo y en el modelo completo debe ser la misma
  if (formula[[2]] != null.model[[2]]) {
    stop("The response in the full and null model does not coincide.\n")
  }
  # tempdir como working directory
  wd <- tempdir()
  # Elimina elementos previos del working directory
  unlink(paste(wd, "*", sep = "/"))
  
  # Evaluamos el modelo nulo
  glm_null <- glm(null.model, 
                  data = data, 
                  family = binomial(link = "logit"),
                  x = TRUE,
                  y = TRUE)
  
  fixed.cov <- dimnames(glm_null$x)[[2]] # covariables fijas en todos los modelos
  
  # Matriz de diseño si hay covariables fijas. Modelo completo:
  glm_full <- glm(formula, 
                  data = data, 
                  family = binomial(link = "logit"),
                  x = TRUE,
                  y = TRUE)
  
  # Variables que son función lineal de otras:
  if (glm_full$rank != dim(glm_full$x)[2]) {
    stop("Some of the explanatory variables are a linear function of others\n")
  }
  
  # Groups:
  X.full <- glm_full$x # Matriz de diseño del modelo completo
  orig.names <- new.names <- colnames(X.full)
  for (i in 1:length(groups)){
    for (j in 1:length(groups[[i]])){
      este <- which(orig.names == groups[[i]][j])
      new.names[este] <- paste("group(",  names(groups)[i], ")", 
                               groups[[i]][j], sep = "")
    }
  }
  colnames(X.full) <- new.names
  namesx <- dimnames(X.full)[[2]] # Nombres actualizados
  # group + nombre del grupo cuando es una variable agrupada
  
  # Compruebo si el modelo nulo está contenido en el completo:
  namesnull <- dimnames(glm_null$x)[[2]]
  "%notin%" <- function(x, table) match(x, table, nomatch = 0) == 0
  
  # El nulo debe tener el intercepto
  if (is.null(namesnull)) stop("The null model should contain the intercept\n")
  
  if (sum(namesnull == "Intercept") == 0 & sum(namesnull == "(Intercept)") == 0) {
    stop("The null model should contain the intercept\n")
  }
  
  for (i in 1:length(namesnull)) {
    if (namesnull[i] %notin% namesx) {
      cat("Error in var: ", namesnull[i], "\n")
      stop("null model is not nested in full model\n")
    }
  }
  
  # Cuando hay variables fijas, compruebo si hay más para seleccionar
  if (length(namesx) == length(namesnull)) {
    stop("The number of fixed covariates is equal to the number of 
         covariates in the full model. No model selection can be done\n")
  }

  fixed.pos <- which(namesx %in% namesnull) # posición de las variables fijas 
  knull <- ncol(glm_null$x) # número de variables fijas en el nulo
  
  # Factores:						
  X1 <- X.full[, -fixed.pos]
  n <- nrow(data)
  if (nrow(X1) < n) stop("NA values found for some of the competing variables")
  
  p <- dim(X1)[2] # Número de covariables a seleccionar
  
  # Groups:
  
  # Creamos una matriz cuyas primeras length(groups) filas sean los grupos y las
  # siguientes sean las variables singulares. Las columnas son las p variables a
  # seleccionar, de forma que la celda (i,j) vale 1 si la variable i pertenece al
  # grupo (si 1<=j<=length(groups)) o variable singular (si j>length(groups)) j.
  
  depvars <- c(names(groups), setdiff(colnames(X.full)[new.names == orig.names],  
                                      attr(glm_null$terms, "term.labels")))
  depvars <- depvars[depvars != "(Intercept)"]
  depvars <- depvars[depvars != "Intercept"] # nombres grupos y vars sing
  
  positions <- matrix(0, ncol = p, nrow = length(depvars))
  for (i in 1:length(depvars)) {
    if (sum(depvars[i] == names(groups)) > 0) {
      positions[i,] <- grepl(paste("group(", depvars[i], sep = ""), 
                             colnames(X1), fixed = TRUE)
    }
    else {positions[i, which(depvars[i] == colnames(X1))] <- 1}
  }
  
  # positionsx tiene un 1 donde hay variables singulares
  positionsx <- as.numeric(colSums(positions %*% t(positions)) == 1)
  rownames(positions) <- depvars
  
  # Escribimos el modelo inicial introducido por el usuario:
  if (is.character(init.model) == TRUE) {
    im <- substr(tolower(init.model), 1, 1)
    if (im != "n" &&
        im != "f" && im != "r") {
      stop("Initial model not valid\n")
    }
    if (im == "n") { # null
      init.model <- rep(0, p)
    }
    if (im == "f") { # full
      init.model <- rep(1, p)
    }
    if (im == "r") { # random
      init.model <- rbinom(n = p,
                           size = 1,
                           prob = .5)
    }
  } else{ # si el usuario da el modelo exacto
    init.model <- as.numeric(init.model > 0)
    if (length(init.model) != p) {
      stop("Initial model with incorrect length\n")
    }
  }
  
  # Info:
  cat("\n Info. . . .\n")
  cat("Most complex model has a total of", nrow(positions) + knull,
      "single numerical covariates and groups.\n")
  # Si hubiera covariables fijas: 
  # (hacerlo después del TFM, en los experimentos no hay cov fijas)
  # if (!is.null(fixed.cov)) {
  # if (knull > 1) {
  #  cat("From those",
  #     knull,
  #    "are fixed and we should select from the remaining",
  #   dim(positions)[1],
  #  "\n")
  # }
  if (knull == 1) {
  cat("From those", knull,
      "is fixed (the intercept) and we should select from the remaining",
      dim(positions)[1], ".\n")
  }
  cat(" Single variables: ", depvars[positionsx == 1], "\n",
      "Group of variables:", depvars[positionsx == 0], "\n\n")
  
  # }
  cat("The problem has a total of", 2^p, "competing models.\n")
  # cat("Of these,", n.burnin + n.iter, "are sampled with replacement.\n")
  cat("Then,", n.iter, # simulo n.burnin + (n.iter * n.thin) modelos y 
      # luego quito n.burnin  y 1 cada n.thin modelos
      "are kept and used to construct the summaries.\n")
  
  # Comprobación previas del espacio de modelos
  nombre_previas <- c("SBSB", "ConstConst", "SBConst", "SB", "Const")
  if(prior.models %notin% nombre_previas) {
    stop("Prior over the model space not supported.\n")
  }
  prior_type <- as.numeric(factor(prior.models, 
                                  levels = nombre_previas, labels = 1:5))
  
  if(prior_type == 4) {
    previasSB <- auxSB(positions, groups)
    # Creo todas las posibles soluciones para resolver i + j_1 +...+ j_p = r
    # haciendolo al principio ahorro tiempo, esto no se puede con SBSB
  }
  
  if (prior.betas != "Robust") {
    stop("Dont recognize the prior for betas.\n It should be:  \"Robust\".\n")
  }
  
  # Muestreo con reemplazamiento:
  
  # Calculo el factor Bayes de cada modelo dividiendo la marginal de los
  # parámetros de dicho modelo entre la marginal del parámetro del modelo nulo,
  # calculadas con BAS.
  
  bas_null <- bas.glm(null.model,
                      data = data, 
                      n.models = 1, # solo el nulo
                      method = "MCMC",  
                      # deterministic no va para null con n.models = 1
                      betaprior = robust(), 
                      family = binomial(),
                      modelprior = beta.binomial(1,1))
  
  gamma <- list() # lista donde guardo los modelos simulados
  BF <- model_priors <- c() # guardo factores Bayes y previas de los modelos
  gamma[[1]] <- init.model # inicializo con el modelo inicial

  log_marginal_nulo <- bas_null$logmarg
  # BayesFactor da error si es null:
  if (sum(gamma[[1]]) == 0) {BF[1] <- 1} else {
    # escribo la formula de gamma[[1]] para pasársela a BAS
    formula_gamma <- formula.gamma(gamma[[1]], formula, orig.names[-1])
    BF[1] <- BayesFactor(formula_gamma, data, log_marginal_nulo)
  }
  
  # Calculo la previa del modelo gamma[[1]]
  model_priors[1] <- switch(prior_type, 
                            SBSB(gamma[[1]], positions, groups), # SB-SB
                            ConstConst(gamma[[1]], positions, groups), #Const-Const
                            SBConst(gamma[[1]], positions, groups), # sB-Const
                            if(sum(gamma[[1]]) == 0){1} else { 
                              previasSB[[sum(gamma[[1]])]]}, #SB
                            1/(2^p)) # Const
  
  for (i in 1:(n.burnin + (n.iter * n.thin) - 1)) {
    for (j in 1:p) {
      new_gamma <- gamma[[i]]
      new_gamma[j] <- 1 - new_gamma[j] # Cambia la componente j-esima
      
      # Calculo el factor Bayes de cada modelo
      if (sum(new_gamma) == 0) {new_BF <- 1} else {
        new_formula <- formula.gamma(new_gamma, formula, orig.names[-1])
        new_BF <- BayesFactor(new_formula, data, log_marginal_nulo)
      }
      
      # Calculo la previa del modelo new_gamma
      prob_new_gamma <- switch(prior_type, 
                               SBSB(new_gamma, positions, groups), 
                               ConstConst(new_gamma, positions, groups),
                               SBConst(new_gamma, positions, groups),
                               if(sum(new_gamma) == 0){1} else {
                                 previasSB[[sum(new_gamma)]]},
                               1/(2^p))
      
      prob_gamma <- model_priors[i]
      
      r <- (new_BF * prob_new_gamma) /
        (new_BF * prob_new_gamma + BF[i] * prob_gamma)
      
      if (as.logical(rbinom(1,1,r))) {
        gamma[[i]] <- new_gamma
        BF[i] <- new_BF
        model_priors[i] <- prob_new_gamma
      }
    }
    # Siguiente iteración
    gamma[[i + 1]] <- gamma[[i]]
    BF[i + 1] <- BF[i]
    model_priors[i + 1] <- model_priors[i]
  }
  
  # n.burnin
  gamma <- gamma[-c(1:n.burnin)]
  BF <- BF[-c(1:n.burnin)]
  model_priors <- model_priors[-c(1:n.burnin)]
  
  # n.thin
  gamma <- gamma[c(1:n.iter * n.thin)]
  BF <- BF[c(1:n.iter * n.thin)]
  model_priors <- model_priors[c(1:n.iter * n.thin)]
  
  if (sum(duplicated(gamma)) != 0) { # Si hay repetidos los quito
    BF <- BF[-which(duplicated(gamma))]
    model_priors <- model_priors[-which(duplicated(gamma))]
    gamma <- unique(gamma)
  }
  
  # Calculo las posteriors
  
  N <- length(gamma)
  denominador <- sum(BF * model_priors)
  posteriors <- (BF * model_priors)/denominador
  
  # Modelo de mayor posterior
  
  bestmodel <- gamma[which.max(posteriors)]
  
  # Probabilidades de inclusión para cada grupo
  
  prob_incl <- c()
  for (j in 1:length(groups)) {
    prob <- 0
    for (i in 1:N) {
      variables <- positions %*% gamma[[i]]
      if (variables[j] > 0) {
        # la prob del grupo es la suma de las probabilidades de los modelos que
        # tienen al menos una variable del grupo activa
        prob <- prob + posteriors[i]
      }
    }
    prob_incl[j] <- prob
  }
  
  # Probabilidades de inclusión para cada variable singular
  
  k <- sum(positionsx)
  if (k != 0) {
    for (j in 1:k) {
      prob <- 0
      for (i in 1:N) {
        variables <- positions %*% gamma[[i]]
        if (variables[length(groups) + j] > 0) {
          # la prob de una variable es la suma de las probs de los modelos que
          # tienen dicha variable activa
          prob <- prob + posteriors[i]
        }
      }
      prob_incl[j + length(groups)] <- prob
      # Las guardo a continuación de los grupos para mostrar luego
    }
  }
  inclusiones <- data.frame(prob_incl, row.names = depvars)
  
  # Probabilidades de inclusión para cada variable agrupada
  
  if (prior_type == 4) { ## SB
    prob_incl <- c()
    for (j in 1:ncol(positions)) {
      prob <- 0
      for (i in 1:N) {
        if (gamma[[i]][j] > 0) {
          # Si no considero la presencia de grupos, la prob de las variables
          # agrupadas se calcula como la de las singulares
          prob <- prob + posteriors[i]
        }
      }
      prob_incl[j] <- prob
    }
  } else{
    L <- c(0, lengths(groups))
    prob_incl <- c() # Las guardo por orden
    for (l in 1:length(groups)) {
      for (j in 1:L[l + 1]) {
        prob <- 0
        # posición de la variable j-ésima del grupo l
        pos <- which(as.logical(positions[l,]))[j]
        for (i in 1:N) {
          if (gamma[[i]][pos] == 1) {
            prob <- prob + posteriors[i]
          }
        }
        if (inclusiones[l, 1] != 0 ) {
          prob_incl[j + sum(L[1:l])] <- prob/inclusiones[l, 1]
        } else {prob_incl[j + sum(L[1:l])] <- 0}
      }
    }
  }
  
  # Matriz para mostrar las probabilidades de inclusión
  # primero grupos con sus respectivas variables y luego las singulares
  salida <- data.frame()
  aux <- 1
  namesord <- c()
  L <- c(0, lengths(groups))
  for (l in 1:length(groups)) {
    salida[aux, 1] <- inclusiones[l, 1]
    salida[aux, 2] <- paste0(" Variables ", names(groups)[l], ":  ")
    namesord[aux] <- depvars[l]
    aux <- aux + 1
    for (j in 1:L[l + 1]) {
      salida[aux, 1] <- "-"
      salida[aux, 2] <- prob_incl[j + sum(L[1:l])]
      namesord[aux] <- groups[[l]][j]
      aux <- aux + 1
    }
  }
  salida[aux, 1] <- "-"
  salida[aux, 2] <- "-"
  namesord[aux] <- "Variables singulares:"
  if (k != 0) {
    salida[(aux + 1):(aux + p - sum(L)), 1] <- inclusiones[(length(groups) + 1):nrow(inclusiones),1]
    salida[(aux + 1):(aux + p - sum(L)), 2] <- rep("-", nrow(inclusiones) - length(groups))
    namesord[(aux + 1):(aux + p - sum(L))] <- depvars[(length(groups) + 1):length(depvars)]
  }
  
  row.names(salida) <- namesord
  colnamessalida <- c("Prob inclusión grupos y vars singulares", "Prob inclusión vars agrupadas")
  
  # # Print de la salida
  # knitr::kable(salida, 
  #             format = "html", 
  #             digits = 6, 
  #             col.names = colnamessalida, 
  #             linesep = "",
  #             align = 'lcc', 
  #             booktabs = TRUE, 
  #             caption = "Probabilidad de inclusión de grupos y variables 
  #               singulares (primera columna) y de variables agrupadas 
  #               (segunda columna).") %>%
  #             row_spec(0, color = "black", angle = -45) %>%
  #             kableExtra::kable_styling(latex_options = c("repeat_header"), 
  #               font_size = 10, 
  #               full_width = FALSE)
  colnames(salida) <- colnamessalida 
  
  tiempo <- proc.time() - t
  
  # Resultado final
  
  result <- list()
  # result$call <- match.call()
  result$time <- tiempo # Tiempo que le llevó ejecutar todo
  # result$glmfull <- glm_full # glm del modelo completo
  # if(!is.null(fixed.cov)){
  #   result$glmnull <- glm_null # The glm object for the null model
  # }
  result$variables <- depvars # Nombre de las variables a seleccionar
  # result$n <- n # Número de observaciones
  # result$p <- length(depvars) # Número de grupos y vars sing a seleccionar
  # result$k <- knull # Número de covariables fijas

  # result$positions <- positions
  # result$positionsx <- positionsx
  
  result$models <- gamma # Modelos visitados quitando n.burnin, thin y duplicados
  result$modelsBF <- BF # Factores bayes correspondientes a cada modelo gamma
  result$modelspriors <- model_priors
  result$modelsposteriors <- posteriors
  result$inclprob <- salida 
  #result$postprobdim <- probdim # Probabilidad de cada dimensión (hacer después del TFM)
  result$bestmodel <- bestmodel
  
  result
}
