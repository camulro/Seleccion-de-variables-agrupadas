#----------------------------------------------------------------------------

# Vamos a simular los datos de una regresión logística de coeficientes conocidos
# para ver cómo de eficaz es la metodología propuesta en la detección de las 
# variables agrupadas verdaderas, modificando diferentes aspectos.

#----------------------------------------------------------------------------

source("GibbsSampling.R")

#----------------------------------------------------------------------------

# Experimento 1: Cambiar el número de vars agrupadas

## Función para simular los datos de una regresión logística a partir de una 
# normal multivariante y los coeficientes especificados en 4.1
sim_log_regrs <- function (size, num_vars, beta_1, beta_2) {
  
  # simulo las variables del predictor lineal
  A <- matrix(runif(num_vars ^ 2, ) * 2 - 1, ncol = num_vars)
  Sigma1 <- t(A) %*% A # matriz semidefinida positiva para normal multivar
  A <- matrix(runif(num_vars ^ 2) * 2 - 1, ncol = num_vars)
  Sigma2 <- t(A) %*% A
  
  # dos grupos con num_vars cada uno
  X1 <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma1)), Sigma = Sigma1)
  X2 <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma2)), Sigma = Sigma2)
  X1 <- scale(X1, center = TRUE, scale = TRUE) # estandariza las vars
  X2 <- scale(X2, center = TRUE, scale = TRUE)
  
  intercepto <- 1
  cita <- intercepto + X1 %*% beta_1 + X2 %*% beta_2 # predictor lineal
  prob <- 1 / (1 + exp(-cita))
  y <- rbinom(n = size, size = 1, prob = prob) # bernuilli con pi = prob
  
  df <- data.frame(y, X1, X2)
  names(df) <- c("Y", "X1.1", "X1.2", "X1.3", "X1.4",
                 "X2.1", "X2.2", "X2.3", "X2.4")
  return(df)
}

## Función auiliar para introducir los grupos espurios
simulador_esp <- function (p, df, size, num_vars) {
  
  if (p < 2) stop("Como mínimo debe haber 3 grupos para comparar.")
  nombres <- names(df)
  for (i in ((ncol(df)-1)/4 + 1):p) {
    A <- matrix(runif(num_vars ^ 2) * 2 - 1, ncol = num_vars) # grupo i espurio
    Sigma <- t(A) %*% A # como antes
    
    X <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
    X <- scale(X, center = TRUE, scale = TRUE) # normaliza las vars
    df <- data.frame(df, X) # añade los nuevos grupos que no generan los datos
    nombres <- c(nombres, paste0("X",i,".1"), paste0("X",i,".2"),
                 paste0("X",i,".3"), paste0("X",i,".4"))
  }
  names(df) <- nombres
  return(df)
}

## Función para realizar la selección de variables con p variables agrupadas 
# en cada grupo para cada previa
simulacion <- function (p, dfsim, num_vars, previas) {
  
  result <- list()
  items <- 1
  groups <- list()
  nombres <- c()
  for (j in 1:p) {
    groups[[j]] <- c(paste0("X",j,".1"), paste0("X",j,".2"), 
                     paste0("X",j,".3"), paste0("X",j,".4"))
    nombres <- c(nombres, paste0("G",j))
  }
  names(groups) <- nombres
  
  formula <- "Y ~ ."
  for(j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
              dfsim,
              null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
              groups = groups,
              prior.betas = "Robust",
              prior.models = previas[j],
              n.iter = 1000,
              init.model = "Full",
              n.burnin = 100,
              n.thin = 1)
    result[[items]] <- gs
    items <- items + 1
    # hace la selección de variables para cada previa
  }
  return(result)
}

## DATOS

size <- 300 # número de observaciones
num_vars <- 4 # número de variables de cada grupo
p <- c(3, 5, 7) # número de grupos
beta_1 <- c(-0.6, -0.56, -0.53, 0.49)
beta_2 <- c(-0.45, 0.4, 0.35, 0.3)
previas <- c("ConstConst", "SBConst", "SBSB")

datalist <- list()
set.seed(123)
for (i in 1:10) { # genero los conjuntos de datos verdaderos que se usarán
  datalist[[i]] <- sim_log_regrs(size = size, num_vars = num_vars,
                                 beta_1 = beta_1, beta_2 = beta_2)
}

dfsimlist <- list()
dfsimlist[[1]] <- dfsimlist[[2]] <- dfsimlist[[3]] <- list()
for (i in 1:10) { # genero el tercer grupo espurio
  dfsimlist[[1]][[i]] <- simulador_esp(p[1], datalist[[i]], size, num_vars) 
}

for (i in 1:10) { # genero el cuarto y quinto grupos espurios
  dfsimlist[[2]][[i]] <- simulador_esp(p[2], dfsimlist[[1]][[i]], size, num_vars) 
}

for (i in 1:10) { # genero el sexto y séptimo grupos espurios
  dfsimlist[[3]][[i]] <- simulador_esp(p[3], dfsimlist[[2]][[i]], size, num_vars) 
}

## Gibbs
# for (i in 1:10) {# 3 grupos
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado1 <- simulacion(p[1], dfsimlist[[1]][[i]], num_vars, previas) 
#   loc <- paste0("./E1/SB/resultado1.", i, ".Rdata")
#   save(resultado1, file = loc)
# }
# 
# for (i in 1:10) {# 5 grupos
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado2 <- simulacion(p[2], dfsimlist[[2]][[i]], num_vars, previas) 
#   loc <- paste0("./E1/SB/resultado2.", i, ".Rdata")
#   save(resultado2, file = loc)
# }
# 
# for (i in 1:10) {# 7 grupos
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado3 <- simulacion(p[3], dfsimlist[[3]][[i]], num_vars, previas) 
#   loc <- paste0("./E1/SB/resultado3.", i, ".Rdata")
#   save(resultado3, file = loc)
# }

previas <- c("SB")
# for (i in 1:10) {# 3 grupos
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado1SB <- simulacion(p[1], dfsimlist[[1]][[i]], num_vars, previas)
#   loc <- paste0("./E1/resultado1.", i, "SB.Rdata")
#   save(resultado1SB, file = loc)
# }
# 
# for (i in 1:10) {# 5 grupos
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado2SB <- simulacion(p[2], dfsimlist[[2]][[i]], num_vars, previas)
#   loc <- paste0("./E1/resultado2.", i, "SB.Rdata")
#   save(resultado2SB, file = loc)
# }
# 
# for (i in 1:10) {# 7 grupos
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado3SB <- simulacion(p[3], dfsimlist[[3]][[i]], num_vars, previas)
#   loc <- paste0("./E1/resultado3.", i, "SB.Rdata")
#   save(resultado3SB, file = loc)
# }

#----------------------------------------------------------------------------
# Experimento 2: Cambiar el número de variables dentro de cada grupo

## Función para simular los datos de una regresión logística a partir de una 
# normal multivariante y los coeficientes especificados en 4.2
sim_log_regrs2 <- function (size, num_vars, beta_1) {
  
  # simulo las variables del predictor lineal
  A <- matrix(runif(num_vars ^ 2) * 2 - 1, ncol = num_vars)
  Sigma1 <- t(A) %*% A # matriz semidefinida posit para normal multivar
  
  X1 <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma1)), Sigma = Sigma1)
  X1 <- scale(X1, center = TRUE, scale = TRUE) # normaliza las vars
  # grupo 1 y único que genera los datos
  
  intercepto <- 1
  cita <- intercepto + X1 %*% beta_1 # predictor lineal
  prob <- 1 / (1 + exp(-cita))
  y <- rbinom(n = size, size = 1, prob = prob)
  
  df <- data.frame(y, X1) # datos originales simulados
  nombres <- c("Y")
  for (j in 1:num_vars) {
    nombres <- c(nombres, paste0("X1.",j))
  }
  names(df) <- nombres
  return(df)
}

## Función auxiliar para introducir los grupos espurios
simulador_esp <- function (p, df, size, num_vars) {
  if (p < 2) stop("Como mínimo debe haber 2 grupos para comparar.")
  nombres <- names(df)
  for (i in ((ncol(df)-1)/15 + 1):p) {
    A <- matrix(runif(num_vars ^ 2) * 2 - 1, ncol = num_vars) # grupo i espurio
    Sigma <- t(A) %*% A # como antes
    
    X <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
    X <- scale(X, center = TRUE, scale = TRUE) # normaliza las vars
    df <- data.frame(df, X) # añade los nuevos grupos que no generan los datos
    
    for (j in 1:num_vars) {
      nombres <- c(nombres, paste0("X",i,".",j))
    }
  }
  names(df) <- nombres
  return(df)
}

## Función para realizar la selección de variables con p grupos de num_vars vars 
# para las dos previas a comparar: SBC y SBSB, en dos bancos de 
# datos simulados diferentes: df1 y df2
simulacion <- function (p, dfsim1, dfsim2, num_vars, previas) {
  
  result <- list()
  items <- 1
  groups <- list()
  nombres <- c()
  for (j in 1:p) {
    nombres_vars <- c()
    for (l in 1:num_vars) {
      nombres_vars <- c(nombres_vars, paste0("X",j,".",l))
    }
    groups[[j]] <- nombres_vars
    nombres <- c(nombres, paste0("G",j))
  }
  names(groups) <- nombres
  
  formula <- "Y ~ ."
  # Primera base
  for(j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
                dfsim1,
                null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
                groups = groups,
                prior.betas = "Robust",
                prior.models = previas[j],
                n.iter = 1000,
                init.model = "Full",
                n.burnin = 100,
                n.thin = 1)
    result[[items]] <- gs
    items <- items + 1
  }
  
  # Segunda base
  for(j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
             dfsim2,
             null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
             groups = groups,
             prior.betas = "Robust",
             prior.models = previas[j],
             n.iter = 1000,
             init.model = "Full",
             n.burnin = 100,
             n.thin = 1)
    result[[items]] <- gs
    items <- items + 1
  }
  
  result[[items]] <- dfsim1
  result[[items + 1]] <- dfsim2
  return(result)
}

## DATOS:
size <- 600 # número de observaciones
num_vars <- c(5, 10, 15) # número de variables de cada grupo
p <- 4 # número de grupos
previas <- c("ConstConst", "SBConst", "SBSB")

beta_11 <- 0.56; beta_12 <- 0.49
beta_21 <- -0.43; beta_22 <- -0.37

set.seed(123)
# Primer banco de datos completo con beta_1=(beta_11, beta_12, 0...0)
beta_1 <- c(beta_11, beta_12, rep(0, num_vars[3]-2))
df1 <- sim_log_regrs2(size = size, num_vars = num_vars[3], beta_1 = beta_1)
dfsim1 <- simulador_esp(p, df1, size, num_vars[3])
# Segundo banco de datos completo con beta_1=(beta_21, beta_22, 0...0)
beta_2 <- c(beta_21, beta_22, rep(0, num_vars[3]-2))
df2 <- sim_log_regrs2(size = size, num_vars = num_vars[3], beta_1 = beta_2)
dfsim2 <- simulador_esp(p, df2, size, num_vars[3])

## RESULTADOS:
## l = 5
indices1 <- c(1,2:6,17:21,32:36,47:51)
# resultado4.1 <- simulacion(p, dfsim1[,indices1], dfsim2[indices1], num_vars[1], previas)
# save(resultado4.1, file = "./E2/resultado4.1.Rdata")

## l = 10
indices2 <- c(1,2:11,17:26,32:41,47:56)
# resultado4.2 <- simulacion(p, dfsim1[,indices2], dfsim2[indices2], num_vars[2], previas)
# save(resultado4.2, file = "./E2/resultado4.2.Rdata")

## l = 15
# resultado4.3 <- simulacion(p, dfsim1, dfsim2, num_vars[3], previas)
# save(resultado4.3, file = "./E2/resultado4.3.Rdata")

### Con SB
previas <- c("SB")
# resultado4.1SB <- simulacion(p, dfsim1[,indices1], dfsim2[indices1], num_vars[1], previas)
# save(resultado4.1SB, file = "./E2/resultado4.1SB.Rdata")

# resultado4.2SB <- simulacion(p, dfsim1[,indices2], dfsim2[indices2], num_vars[2], previas)
# save(resultado4.2SB, file = "./E2/resultado4.2SB.Rdata")

# resultado4.3SB <- simulacion(p, dfsim1, dfsim2, num_vars[3], previas)
# save(resultado4.3SB, file = "./E2/resultado4.3SB.Rdata")

## DATOS 2  
size <- 600 # número de observaciones
num_vars <- c(5, 10, 15) # número de variables de cada grupo
p <- 4 # número de grupos
previas <- c("ConstConst", "SBConst", "SBSB")

betas_1 <- c(0.63, 0.58, -0.53, 0.5, -0.46)
betas_2 <- c(-0.41, 0.39, -0.35, 0.3, 0.27)

set.seed(123)
# Primer banco de datos completo con beta_1=rep(betas_1,3)
beta_1 <- rep(betas_1, 3)
df1 <- sim_log_regrs2(size = size, num_vars = num_vars[3], beta_1 = beta_1)
dfsim1 <- simulador_esp(p, df1, size, num_vars[3])
# Segundo banco de datos completo con beta_1=rep(betas_2,3)
beta_2 <- rep(betas_2, 3)
df2 <- sim_log_regrs2(size = size, num_vars = num_vars[3], beta_1 = beta_2)
dfsim2 <- simulador_esp(p, df2, size, num_vars[3])

## l = 5
indices1 <- c(1,2:6,17:21,32:36,47:51)
# resultado4.4 <- simulacion(p, dfsim1[,indices1], dfsim2[,indices1], num_vars[1], previas)
# save(resultado4.4, file = "./E2/resultado4.4.Rdata")

## l = 10
indices2 <- c(1,2:11,17:26,32:41,47:56)
# resultado4.5 <- simulacion(p, dfsim1[,indices2], dfsim2[,indices2], num_vars[2], previas)
# save(resultado4.5, file = "./E2/resultado4.5.Rdata")

## l = 15
# resultado4.6 <- simulacion(p, dfsim1, dfsim2, num_vars[3], previas)
# save(resultado4.6, file = "./E2/resultado4.6.Rdata")

### Con SB
previas <- c("SB")
# resultado4.4SB <- simulacion(p, dfsim1[,indices1], dfsim2[,indices1], num_vars[1], previas)
# save(resultado4.4SB, file = "./E2/resultado4.4SB.Rdata")

# resultado4.5SB <- simulacion(p, dfsim1[,indices2], dfsim2[,indices2], num_vars[2], previas)
# save(resultado4.5SB, file = "./E2/resultado4.5SB.Rdata")

# resultado4.6SB <- simulacion(p, dfsim1, dfsim2, num_vars[3], previas)
# save(resultado4.6SB, file = "./E2/resultado4.6SB.Rdata")

#----------------------------------------------------------------------------

# Experimento 3: Comparar distintas metodologias

## Función para simular los datos de una regresión logística a partir de una 
# normal multivariante y los coeficientes especificados en 4
sim_log_regrs3 <- function (size, num_vars, beta_1, beta_2) {
  
  A <- matrix(runif(num_vars[1] ^ 2) * 2 - 1, ncol = num_vars[1])
  Sigma1 <- t(A) %*% A
  A <- matrix(runif(num_vars[2] ^ 2) * 2 - 1, ncol = num_vars[2]) 
  Sigma2 <- t(A) %*% A
  
  # genero las variables del predictor lineal a partir de una normal multivar
  X1 <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma1)), Sigma = Sigma1)
  X2 <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma2)), Sigma = Sigma2)
  X1 <- scale(X1, center = TRUE, scale = TRUE)
  X2 <- scale(X2, center = TRUE, scale = TRUE)
  
  intercepto <- 1
  cita <- intercepto + X1 %*% beta_1 + X2 %*% beta_2 # predictor lineal
  prob <- 1 / (1 + exp(-cita))
  y <- rbinom(n = size, size = 1, prob = prob) 
  
  df <- data.frame(y, X1, X2) # datos simulados originales
  nombres <-c("Y")
  for (i in 1:2){
    for(j in 1:num_vars[i]){
      nombres <- c(nombres, paste0("X",i,".",j))
    }
  }
  names(df) <- nombres
  return(df)
}
## Función auxiliar para introducir los grupos espurios
simulador_esp <- function (p, df, size, num_vars) {
  if (p < 3) stop("Como mínimo debe haber 3 grupos para comparar.")
  nombres <- names(df)
  for (i in 3:p) {
    A <- matrix(runif(num_vars[i] ^ 2) * 2 - 1, ncol = num_vars[i])
    Sigma <- t(A) %*% A # como antes
    
    X <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
    X <- scale(X, center = TRUE, scale = TRUE)
    df <- data.frame(df, X) # banco de datos con vars espurias
    
    for (j in 1:num_vars[i]) {
      nombres <- c(nombres, paste0("X",i,".",j))
    }
  }
  names(df) <- nombres
  return(df)
}
## Función para realizar la selección de variables con p variables 
# en cada grupo para la previa SBSB, en dos bancos de datos simulados: 
# df1 y df2, y para un número diferente de observaciones
simulacion <- function (p, dfsim1, dfsim2, num_vars, previas) {
  
  result <- list()
  items <- 1
  groups <- list()
  nombres <- c()
  for (j in 1:p) {
    nombres_vars <- c()
    for (l in 1:num_vars[j]) {
      nombres_vars <- c(nombres_vars, paste0("X",j,".",l))
    }
    groups[[j]] <- nombres_vars
    nombres <- c(nombres, paste0("G",j))
  }
  names(groups) <- nombres
  
  formula <- "Y ~ ."
  for (j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
             dfsim1,
             null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
             groups = groups,
             prior.betas = "Robust",
             prior.models = previas[j],
             n.iter = 1000,
             init.model = "Full",
             n.burnin = 100,
             n.thin = 2)
    result[[items]] <- gs
    items <- items + 1
  }
  bas1 <- bas.glm(formula,
                  data = dfsim1, 
                  n.models = 1000,
                  method = "MCMC",
                  betaprior = robust(), 
                  family = binomial())
  result[[items]] <- bas1
  items <- items + 1
  
  for (j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
             dfsim2,
             null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
             groups = groups,
             prior.betas = "Robust",
             prior.models = previas[j],
             n.iter = 1000,
             init.model = "Full",
             n.burnin = 100,
             n.thin = 2)
    result[[items]] <- gs
    items <- items + 1
  }
  bas2 <- bas.glm(formula,
                  data = dfsim2, 
                  n.models = 1000,
                  method = "MCMC",
                  betaprior = robust(), 
                  family = binomial(),
                  thin = 1,
  )
  result[[items]] <- bas2
  items <- items + 1
  
  result[[items]] <- dfsim1
  result[[items + 1]] <- dfsim2
  # guardo también los bancos de datos simulados
  return(result)
}
n <- c(300, 600)
num_vars <- c(8, 4, 8, 4)
p <- length(num_vars)
previas <- c("ConstConst", "SBConst",  "SBSB")

# Genero los datos simulados
set.seed(123)
listadf1 <- listadf2 <- listasimdf1 <- listasimdf2 <- list()
listadf1[[1]] <- listadf1[[2]] <- listadf2[[1]] <- listadf2[[2]] <- list()
listasimdf1[[1]] <- listasimdf1[[2]] <- listasimdf2[[1]] <- listasimdf2[[2]] <- list()
for (i in 1:10) {
  beta_1 <- c(0, 0, 0.55, 0.55, 0.55, 0.55, -0.49, -0.49)
  beta_2 <- c(0, 0, 0.46, 0.46)
  # df1 con beta_1, beta_2 y n=300
  listadf1[[1]][[i]] <- sim_log_regrs3(size = n[1], num_vars = num_vars[1:2],
                                       beta_1 = beta_1, beta_2 = beta_2)
  listasimdf1[[1]][[i]] <- simulador_esp(p, listadf1[[1]][[i]], n[1], num_vars)
  # df1 con beta_1, beta_2 y n=600
  listadf1[[2]][[i]] <- sim_log_regrs3(size = n[2], num_vars = num_vars[1:2],
                                       beta_1 = beta_1, beta_2 = beta_2)
  listasimdf1[[2]][[i]] <- simulador_esp(p, listadf1[[2]][[i]], n[2], num_vars)

  beta_1 <- c(0.42, rep(0, 6), -0.4)
  beta_2 <- c(0.35, 0, 0, -0.3)
  # df2 con beta_1, beta_2 y n=300
  listadf2[[1]][[i]] <- sim_log_regrs3(size = n[1], num_vars = num_vars[1:2],
                                       beta_1 = beta_1, beta_2 = beta_2)
  listasimdf2[[1]][[i]] <- simulador_esp(p, listadf2[[1]][[i]], n[1], num_vars)
  # df2 con beta_1, beta_2 y n=600
  listadf2[[2]][[i]] <- sim_log_regrs3(size = n[2], num_vars = num_vars[1:2],
                                       beta_1 = beta_1, beta_2 = beta_2)
  listasimdf2[[2]][[i]] <- simulador_esp(p, listadf2[[2]][[i]], n[2], num_vars)
}

## RESULTADOS
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado5.1 <- simulacion(p, listasimdf1[[1]][[i]], listasimdf2[[1]][[i]],
#                              num_vars, previas)
#   loc <- paste0("./E3/resultado5.1.", i, ".Rdata")
#   save(resultado5.1, file = loc)
# }
# 
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado5.2 <- simulacion(p, listasimdf1[[2]][[i]], listasimdf2[[2]][[i]],
#                              num_vars, previas)
#   loc <- paste0("./E3/resultado5.2.", i, ".Rdata")
#   save(resultado5.2, file = loc)
# }

## Con SB
previas <- c("SB")
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado5.1SB <- simulacion(p, listasimdf1[[1]][[i]], listasimdf2[[1]][[i]],
#                                num_vars, previas)
#   loc <- paste0("./E3/resultado5.1.", i, "SB.Rdata")
#   save(resultado5.1SB, file = loc)
# }
# 
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   resultado5.2SB <- simulacion(p, listasimdf1[[2]][[i]], listasimdf2[[2]][[i]],
#                                num_vars, previas)
#   loc <- paste0("./E3/resultado5.2.", i, "SB.Rdata")
#   save(resultado5.2SB, file = loc)
# }

#----------------------------------------------------------------------------

# Experimento 4: Cómo afecta la correlación y otras consideraciones

## Función para simular los datos de una regresión logística a partir de una 
# normal multivariante 
sim_log_regrs4 <- function (size, num_vars, rho, mubetas, sdbetas) {
  
  # Queremos Cov(X_i, X_j) = 0.9, Var(X_i) = Var(X_j) = 1 -> rho(X_i,X_j) = 0.9
  A <- matrix(rep(rho, num_vars ^ 2), ncol = num_vars) 
  diag(A) <- rep(1, num_vars)
  
  # genero las variables del predictor lineal a partir de una normal multivar
  X1 <- MASS::mvrnorm(size, mu = rep(0, nrow(A)), Sigma = A)
  X2 <- MASS::mvrnorm(size, mu = rep(0, nrow(A)), Sigma = A)
  X3 <- MASS::mvrnorm(size, mu = rep(0, nrow(A)), Sigma = A)
  X1 <- scale(X1, center = TRUE, scale = TRUE)
  X2 <- scale(X2, center = TRUE, scale = TRUE)
  X2 <- scale(X2, center = TRUE, scale = TRUE)
  
  beta_1 <- rnorm(num_vars, mean = mubetas, sd = sdbetas[1])
  beta_2 <- rnorm(num_vars, mean = mubetas, sd = sdbetas[2])
  beta_3 <- rnorm(num_vars, mean = mubetas, sd = sdbetas[3])
  
  intercepto <- 1
  cita <- intercepto + X1 %*% beta_1 + X2 %*% beta_2 + X3 %*% beta_3 # predictor lineal
  prob <- 1 / (1 + exp(-cita))
  y <- rbinom(n = size, size = 1, prob = prob) 
  
  df <- data.frame(y, X1, X2, X3) # datos simulados originales
  nombres <-c("Y")
  for (i in 1:3){
    for(j in 1:num_vars){
      nombres <- c(nombres, paste0("X",i,".",j))
    }
  }
  names(df) <- nombres
  return(df)
}

## Función para realizar la selección de variables con p grupos y num_vars 
# variables en cada grupo
simulacion <- function (p, df, num_vars, previas) {
  
  result <- list()
  items <- 1
  groups <- list()
  nombres <- c()
  for (j in 1:p) {
    nombres_vars <- c()
    for (l in 1:num_vars) {
      nombres_vars <- c(nombres_vars, paste0("X",j,".",l))
    }
    groups[[j]] <- nombres_vars
    nombres <- c(nombres, paste0("G",j))
  }
  names(groups) <- nombres
  
  formula <- "Y ~ ."
  for (j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
             df,
             null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
             groups = groups,
             prior.betas = "Robust",
             prior.models = previas[j],
             n.iter = 1000,
             init.model = "Full",
             n.burnin = 100,
             n.thin = 1)
    result[[items]] <- gs
    items <- items + 1
  }
  
  return(result)
}

n <- 600
num_vars <- 20
p <- 3 # número de grupos reales
rho <- 0.9
mubetas <- rep(0.5, 3)
sdbetas <- c(0.1, 0.15, 0.2)
previas <- c("ConstConst", "SBConst",  "SBSB")

dflist <- list()
set.seed(123)
for (i in 1:10) {
  dflist[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars, 
                                rho = rho, mubetas = mubetas, sdbetas = sdbetas)
}

# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultado6.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado6.1.", i,".Rdata")
#   save(resultado6.1, file = loc)
#   
#   resultado6.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado6.2.", i,".Rdata")
#   save(resultado6.2, file = loc)
#   
#   resultado6.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado6.3.", i,".Rdata")
#   save(resultado6.3, file = loc)
# }

## Con SB
# previas <- c("SB")
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultado6.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado6.1.", i,"SB.Rdata")
#   save(resultado6.1SB, file = loc)
# 
#   resultado6.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado6.2.", i,"SB.Rdata")
#   save(resultado6.2SB, file = loc)
# 
#   resultado6.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado6.3.", i,"SB.Rdata")
#   save(resultado6.3SB, file = loc)
# }

## Con BAS
# library(BAS)
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   formula <- Y ~ .
#   resultado6.1.BAS <-  bas.glm(formula,
#                                data = df1,
#                                family = binomial(link = "logit"),
#                                betaprior = robust(),
#                                method = "BAS",
#                                laplace = TRUE)
#   loc <- paste0("./E4/resultado6.1.BAS.", i,".Rdata")
#   save(resultado6.1.BAS, file = loc)
# 
#   resultado6.2.BAS <- bas.glm(formula,
#                               data = df2,
#                               family = binomial(link = "logit"),
#                               betaprior = robust(),
#                               method = "BAS",
#                               laplace = TRUE)
#   loc <- paste0("./E4/resultado6.2.BAS.", i,".Rdata")
#   save(resultado6.2.BAS, file = loc)
# 
#   resultado6.3.BAS <- bas.glm(formula,
#                               data = df3,
#                               family = binomial(link = "logit"),
#                               betaprior = robust(),
#                               method = "BAS",
#                               laplace = TRUE)
#   loc <- paste0("./E4/resultado6.3.BAS.", i,".Rdata")
#   save(resultado6.3.BAS, file = loc)
# }

## Experimento 4.1: para rho = 0.5, 0.6, 0.7, 0.8

n <- 600 
num_vars <- 20 
p <- 3 # número de grupos reales
rho <- c(0.5, 0.6, 0.7, 0.8)
mubetas <- rep(0.5, 3)
sdbetas <- c(0.1, 0.15, 0.2)
previas <- c("ConstConst", "SBConst",  "SBSB")

# set.seed(123)
# dflistrho1 <- dflistrho2 <- dflistrho3 <- dflistrho4 <- list()
# for(i in 1:5) {
#   dflistrho1[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars, 
#                                     rho = rho[1], mubetas = mubetas, sdbetas = sdbetas)
#   dflistrho2[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars, 
#                                     rho = rho[2], mubetas = mubetas, sdbetas = sdbetas)
#   dflistrho3[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars, 
#                                     rho = rho[3], mubetas = mubetas, sdbetas = sdbetas)
#   dflistrho4[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars, 
#                                     rho = rho[4], mubetas = mubetas, sdbetas = sdbetas)
# }
# 
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho1[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho1.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho1.1.", i,".Rdata")
#   save(resultadorho1.1, file = loc)
#   
#   resultadorho1.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho1.2.", i,".Rdata")
#   save(resultadorho1.2, file = loc)
#   
#   resultadorho1.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho1.3.", i,".Rdata")
#   save(resultadorho1.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho2[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho1.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho2.1.", i,".Rdata")
#   save(resultadorho1.1, file = loc)
#   
#   resultadorho1.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho2.2.", i,".Rdata")
#   save(resultadorho1.2, file = loc)
#   
#   resultadorho1.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho2.3.", i,".Rdata")
#   save(resultadorho1.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho3[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho1.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho3.1.", i,".Rdata")
#   save(resultadorho1.1, file = loc)
#   
#   resultadorho1.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho3.2.", i,".Rdata")
#   save(resultadorho1.2, file = loc)
#   
#   resultadorho1.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho3.3.", i,".Rdata")
#   save(resultadorho1.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho4[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho1.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho4.1.", i,".Rdata")
#   save(resultadorho1.1, file = loc)
#   
#   resultadorho1.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho4.2.", i,".Rdata")
#   save(resultadorho1.2, file = loc)
#   
#   resultadorho1.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho4.3.", i,".Rdata")
#   save(resultadorho1.3, file = loc)
# }

## Con SB
previas <- c("SB")
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho1[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho1.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho1.1.", i,"SB.Rdata")
#   save(resultadorho1.1SB, file = loc)
# 
#   resultadorho1.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho1.2.", i,"SB.Rdata")
#   save(resultadorho1.2SB, file = loc)
# 
#   resultadorho1.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho1.3.", i,"SB.Rdata")
#   save(resultadorho1.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho2[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho1.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho2.1.", i,"SB.Rdata")
#   save(resultadorho1.1SB, file = loc)
# 
#   resultadorho1.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho2.2.", i,"SB.Rdata")
#   save(resultadorho1.2SB, file = loc)
# 
#   resultadorho1.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho2.3.", i,"SB.Rdata")
#   save(resultadorho1.3SB, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho3[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho1.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho3.1.", i,"SB.Rdata")
#   save(resultadorho1.1SB, file = loc)
# 
#   resultadorho1.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho3.2.", i,"SB.Rdata")
#   save(resultadorho1.2SB, file = loc)
# 
#   resultadorho1.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho3.3.", i,"SB.Rdata")
#   save(resultadorho1.3SB, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho4[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho1.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho4.1.", i,"SB.Rdata")
#   save(resultadorho1.1SB, file = loc)
# 
#   resultadorho1.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho4.2.", i,"SB.Rdata")
#   save(resultadorho1.2SB, file = loc)
# 
#   resultadorho1.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho4.3.", i,"SB.Rdata")
#   save(resultadorho1.3SB, file = loc)
# }

## Experimento 4.2: Cambiar los coeficientes

n <- 600 
num_vars <- 20 
p <- 3 # número de grupos reales
rho <- c(0.5, 0.6, 0.7, 0.8)
mubetas <- rep(0.15, 3)
sdbetas <- c(0.1, 0.15, 0.2)
previas <- c("ConstConst", "SBConst",  "SBSB")

# set.seed(123)
# dflistrho1 <- dflistrho2 <- dflistrho3 <- dflistrho4 <- list()
# for(i in 1:5) {
#   dflistrho1[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars,
#                                     rho = rho[1], mubetas = mubetas, sdbetas = sdbetas)
#   dflistrho2[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars,
#                                     rho = rho[2], mubetas = mubetas, sdbetas = sdbetas)
#   dflistrho3[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars,
#                                     rho = rho[3], mubetas = mubetas, sdbetas = sdbetas)
#   dflistrho4[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars,
#                                     rho = rho[4], mubetas = mubetas, sdbetas = sdbetas)
# }

# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho1[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho2.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho1.1.", i,"c2.Rdata")
#   save(resultadorho2.1, file = loc)
# 
#   resultadorho2.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho1.2.", i,"c2.Rdata")
#   save(resultadorho2.2, file = loc)
# 
#   resultadorho2.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho1.3.", i,"c2.Rdata")
#   save(resultadorho2.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho2[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho2.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho2.1.", i,"c2.Rdata")
#   save(resultadorho2.1, file = loc)
#   
#   resultadorho2.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho2.2.", i,"c2.Rdata")
#   save(resultadorho2.2, file = loc)
#   
#   resultadorho2.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho2.3.", i,"c2.Rdata")
#   save(resultadorho2.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho3[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho2.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho3.1.", i,"c2.Rdata")
#   save(resultadorho2.1, file = loc)
#   
#   resultadorho2.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho3.2.", i,"c2.Rdata")
#   save(resultadorho2.2, file = loc)
#   
#   resultadorho2.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho3.3.", i,"c2.Rdata")
#   save(resultadorho2.3, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho4[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho2.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho4.1.", i,"c2.Rdata")
#   save(resultadorho2.1, file = loc)
#   
#   resultadorho2.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho4.2.", i,"c2.Rdata")
#   save(resultadorho2.2, file = loc)
#   
#   resultadorho2.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho4.3.", i,"c2.Rdata")
#   save(resultadorho2.3, file = loc)
# }

## Con SB
# previas <- c("SB")
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho1[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultadorho2.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho1.1.", i,"c2SB.Rdata")
#   save(resultadorho2.1SB, file = loc)
#   
#   resultadorho2.2SB  <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho1.2.", i,"c2SB.Rdata")
#   save(resultadorho2.2SB, file = loc)
#   
#   resultadorho2.3SB  <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho1.3.", i,"c2SB.Rdata")
#   save(resultadorho2.3SB, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho2[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho2.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho2.1.", i,"c2SB.Rdata")
#   save(resultadorho2.1SB, file = loc)
# 
#   resultadorho2.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho2.2.", i,"c2SB.Rdata")
#   save(resultadorho2.2SB, file = loc)
# 
#   resultadorho2.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho2.3.", i,"c2SB.Rdata")
#   save(resultadorho2.3SB, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho3[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho2.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho3.1.", i,"c2SB.Rdata")
#   save(resultadorho2.1SB, file = loc)
# 
#   resultadorho2.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho3.2.", i,"c2SB.Rdata")
#   save(resultadorho2.2SB, file = loc)
# 
#   resultadorho2.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho3.3.", i,"c2SB.Rdata")
#   save(resultadorho2.3SB, file = loc)
# }
# 
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflistrho4[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultadorho2.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultadorho4.1.", i,"c2SB.Rdata")
#   save(resultadorho2.1SB, file = loc)
# 
#   resultadorho2.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultadorho4.2.", i,"c2SB.Rdata")
#   save(resultadorho2.2SB, file = loc)
# 
#   resultadorho2.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultadorho4.3.", i,"c2SB.Rdata")
#   save(resultadorho2.3SB, file = loc)
# }

## Experimento 4.3: Correlaciones pequeñas

## Función para simular los datos de una regresión logística a partir de una 
# normal multivariante 
sim_log_regrs4.3 <- function (size, num_vars, rho, mubetas, sdbetas) {
  
  # Queremos Cov(X_i, X_j) = 0.9, Var(X_i) = Var(X_j) = 1 -> rho(X_i,X_j) = 0.9
  A <- matrix(rep(rho[1], num_vars ^ 2), ncol = num_vars) 
  diag(A) <- rep(1, num_vars)
  B <- matrix(rep(rho[2], num_vars ^ 2), ncol = num_vars) 
  diag(B) <- rep(1, num_vars)
  
  # genero las variables del predictor lineal a partir de una normal multivar
  X1 <- MASS::mvrnorm(size, mu = rep(0, nrow(A)), Sigma = A)
  X2 <- MASS::mvrnorm(size, mu = rep(0, nrow(A)), Sigma = B)
  X3 <- MASS::mvrnorm(size, mu = rep(0, nrow(A)), Sigma = B)
  X1 <- scale(X1, center = TRUE, scale = TRUE)
  X2 <- scale(X2, center = TRUE, scale = TRUE)
  X2 <- scale(X2, center = TRUE, scale = TRUE)
  
  beta_1 <- rnorm(num_vars, mean = mubetas, sd = sdbetas[1])
  beta_2 <- rnorm(num_vars, mean = mubetas, sd = sdbetas[2])
  beta_3 <- rnorm(num_vars, mean = mubetas, sd = sdbetas[3])
  
  intercepto <- 1
  cita <- intercepto + X1 %*% beta_1 + X2 %*% beta_2 + X3 %*% beta_3 # predictor lineal
  prob <- 1 / (1 + exp(-cita))
  y <- rbinom(n = size, size = 1, prob = prob) 
  
  df <- data.frame(y, X1, X2, X3) # datos simulados originales
  nombres <-c("Y")
  for (i in 1:3){
    for(j in 1:num_vars){
      nombres <- c(nombres, paste0("X",i,".",j))
    }
  }
  names(df) <- nombres
  return(df)
}

n <- 600 
num_vars <- 20
p <- 3 # número de grupos reales
rho <- c(0.2, 0.8)
mubetas <- rep(0.5, 3)
sdbetas <- c(0.1, 0.15, 0.2)
previas <- c("ConstConst", "SBConst",  "SBSB")

# dflist <- list()
# set.seed(123)
# for (i in 1:10) {
#   dflist[[i]] <- sim_log_regrs4.3(size = n, num_vars = num_vars, 
#                                   rho = rho, mubetas = mubetas, sdbetas = sdbetas)
# }
# 
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultado7.1 <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado7.1.", i,".Rdata")
#   save(resultado7.1, file = loc)
#   
#   resultado7.2 <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado7.2.", i,".Rdata")
#   save(resultado7.2, file = loc)
#   
#   resultado7.3 <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado7.3.", i,".Rdata")
#   save(resultado7.3, file = loc)
# }

## Con SB
# previas <- c("SB")
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultado7.1SB <- simulacion(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado7.1.", i,"SB.Rdata")
#   save(resultado7.1SB, file = loc)
# 
#   resultado7.2SB <- simulacion(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado7.2.", i,"SB.Rdata")
#   save(resultado7.2SB, file = loc)
# 
#   resultado7.3SB <- simulacion(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado7.3.", i,"SB.Rdata")
#   save(resultado7.3SB, file = loc)
# }

## Experimento 4.4: Grupos mal especificados

simulacion4 <- function (p, df, num_vars, previas) {
  
  result <- list()
  items <- 1
  groups <- list()
  nombres <- c()
  for (j in 1:p) {
    nombres_vars <- c()
    for (l in 1:num_vars) {
      nombres_vars <- c(nombres_vars, paste0("X",j,".",l))
    }
    groups[[j]] <- nombres_vars
    nombres <- c(nombres, paste0("G",j))
  }
  names(groups) <- nombres
  
  # Identificar incorrectamente los grupos
  groups[[1]] <- c(groups[[1]], "X2.1", "X2.2", "X3.1", "X3.2", "X3.3")
  groups[[2]] <- groups[[2]][-c(1:2)]
  groups[[3]] <- groups[[3]][-c(1:3)]
  
  formula <- "Y ~ ."
  for (j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
             df,
             null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
             groups = groups,
             prior.betas = "Robust",
             prior.models = previas[j],
             n.iter = 1000,
             init.model = "Full",
             n.burnin = 100,
             n.thin = 1)
    result[[items]] <- gs
    items <- items + 1
  }
  
  return(result)
}

n <- 600
num_vars <- 20
p <- 3 # número de grupos reales
rho <- c(0.75)
mubetas <- rep(0.5, 3)
sdbetas <- c(0.1, 0.15, 0.2)
previas <- c("ConstConst", "SBConst",  "SBSB")

# dflist <- list()
# set.seed(123)
# for (i in 1:10) {
#   dflist[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars,
#                                 rho = rho, mubetas = mubetas, sdbetas = sdbetas)
# }
# 
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
# 
#   resultado8.1 <- simulacion4(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado8.1.", i,".Rdata")
#   save(resultado8.1, file = loc)
# 
#   resultado8.2 <- simulacion4(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado8.2.", i,".Rdata")
#   save(resultado8.2, file = loc)
# 
#   resultado8.3 <- simulacion4(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado8.3.", i,".Rdata")
#   save(resultado8.3, file = loc)
# }

## Con SB
previas <- c("SB")
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:10) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dflist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50)]
#   
#   resultado8.1SB <- simulacion4(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado8.1.", i,"SB.Rdata")
#   save(resultado8.1SB, file = loc)
#   
#   resultado8.2SB <- simulacion4(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado8.2.", i,"SB.Rdata")
#   save(resultado8.2SB, file = loc)
#   
#   resultado8.3SB <- simulacion4(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado8.3.", i,"SB.Rdata")
#   save(resultado8.3SB, file = loc)
# }

## Experimento 4.5: Definición errónea del modelo y los grupos

### Función auxiliar para introducir los grupos espurios
simulador_esp4 <- function (p, df, size, num_vars) {
  
  if (p < 4) stop("Como mínimo debe haber 4 grupos para comparar.")
  nombres <- names(df)
  for (i in 4:p) {
    A <- matrix(runif(num_vars ^ 2) * 2 - 1, ncol = num_vars)
    Sigma <- t(A) %*% A
    
    X <- MASS::mvrnorm(size, mu = rep(0, nrow(Sigma)), Sigma = Sigma)
    X <- scale(X, center = TRUE, scale = TRUE)
    df <- data.frame(df, X) # banco de datos con vars espurias
    
    for (j in 1:num_vars) {
      nombres <- c(nombres, paste0("X",i,".",j))
    }
  }
  names(df) <- nombres
  return(df)
}

## Función para realizar la selección de variables con p grupos y num_vars
# variables en cada grupo
simulacion5 <- function (p, dfsim, num_vars, previas) {
  
  result <- list()
  items <- 1
  groups <- list()
  nombres <- c()
  for (j in 1:p) {
    nombres_vars <- c()
    for (l in 1:num_vars) {
      nombres_vars <- c(nombres_vars, paste0("X",j,".",l))
    }
    groups[[j]] <- nombres_vars
    nombres <- c(nombres, paste0("G",j))
  }
  names(groups) <- nombres
  
  formula <- "Y ~ ."
  for (j in 1:length(previas)) {
    cat("\n Previa:", previas[j], ":\n")
    gs <- GS(formula, 
             dfsim,
             null.model = paste(as.formula(formula)[[2]], " ~ 1", sep=""),
             groups = groups,
             prior.betas = "Robust",
             prior.models = previas[j],
             n.iter = 1000,
             init.model = "Full",
             n.burnin = 100,
             n.thin = 1)
    result[[items]] <- gs
    items <- items + 1
  }
  
  return(result)
}

n <- 800 
num_vars <- 20
p <- 4 # número de grupos reales
rho <- 0.9
mubetas <- rep(0.5, 3)
sdbetas <- c(0.1, 0.15, 0.2)
previas <- c("ConstConst", "SBConst",  "SBSB")

# dflist <- list()
# set.seed(123)
# for (i in 1:5) {
#   dflist[[i]] <- sim_log_regrs4(size = n, num_vars = num_vars,
#                                 rho = rho, mubetas = mubetas, sdbetas = sdbetas)
# }
# dfsimlist <- list()
# for (i in 1:5) {
#   dfsimlist[[i]] <- simulador_esp4(p, dflist[[i]], size = n, num_vars)
# }

# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dfsimlist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46, 62:66)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48, 62:68)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50, 62:70)]
#   
#   resultado9.1 <- simulacion5(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado9.1.", i,".Rdata")
#   save(resultado9.1, file = loc)
#   
#   resultado9.2 <- simulacion5(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado9.2.", i,".Rdata")
#   save(resultado9.2, file = loc)
#   
#   resultado9.3 <- simulacion5(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado9.3.", i,".Rdata")
#   save(resultado9.3, file = loc)
# }

## Con SB
previas <- c("SB")
# num_vars1 <- 5
# num_vars2 <- 7
# num_vars3 <- 9
# for (i in 1:5) {
#   cat("\n Iteración", i, "-ésima:\n")
#   df <- dfsimlist[[i]]
#   df1 <- df[, c(1,2:6, 22:26, 42:46, 62:66)]
#   df2 <- df[, c(1,2:8, 22:28, 42:48, 62:68)]
#   df3 <- df[, c(1,2:10, 22:30, 42:50, 62:70)]
#   
#   resultado9.1SB <- simulacion5(p, df1, num_vars1, previas)
#   loc <- paste0("./E4/resultado9.1.", i,"SB.Rdata")
#   save(resultado9.1SB, file = loc)
#   
#   resultado9.2SB <- simulacion5(p, df2, num_vars2, previas)
#   loc <- paste0("./E4/resultado9.2.", i,"SB.Rdata")
#   save(resultado9.2SB, file = loc)
#   
#   resultado9.3SB <- simulacion5(p, df3, num_vars3, previas)
#   loc <- paste0("./E4/resultado9.3.", i,"SB.Rdata")
#   save(resultado9.3SB, file = loc)
# }
