## Wynand van Staden
## Model - Data simulation & Maximum Likelihood estimation
## Copyright 2019

## Please note this file makes use of the sirModel functions in the 
## 1. One parameter calibration Least Squares.r file

library(stats4)
library(dplyr)
library(gmp)


###############################################################################################
## Model Calibration
###############################################################################################

trueGamma <- 0.02
calibModelRuns <- 1000

#NB! set the sample size according the sample size in the sirModel() functions
samSize <- 1000

#Initializing variabes to store results
MLEinfPop2 <- matrix(c(0, 0), calibModelRuns, 2)
MLEinfPop4 <- matrix(replicate(4, 0), calibModelRuns, 4)
MLEinfPop64 <- matrix(replicate(64, 0), calibModelRuns, 64)

M.result2 <- list()
M.result4 <- list()
M.result64 <- list()

MLE.result2 <- matrix(c(0), calibModelRuns, 1)
MLE.result4 <- matrix(c(0), calibModelRuns, 1)
MLE.result64 <- matrix(c(0), calibModelRuns, 1)


print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {

  MLEinfPop2[i,] <- sirModel2(trueGamma)
  MLEinfPop4[i,] <- sirModel4(trueGamma)
  MLEinfPop64[i,] <- sirModel64(trueGamma)
  
  neg.logL2 <- function(p1)  {
    p <- sirModel2(p1)
    p <- ifelse(p==0, 0.0001, p)
    x2 <- MLEinfPop2[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x2[i])) + x2[i]*log(p[i]) + (samSize-x2[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    }
    
    as.numeric(-log.lik)
  }
  
  neg.logL4 <- function(p1)  {
    p <- sirModel4(p1)
    p <- ifelse(p==0, 0.0001, p)
    x4 <- MLEinfPop4[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x4[i])) + x4[i]*log(p[i]) + (samSize-x4[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    
    }
    
    as.numeric(-log.lik)
  }
  
  neg.logL64 <- function(p1)  {
    p <- sirModel64(p1)
    p <- ifelse(p==0, 0.0001, p)
    x64 <- MLEinfPop64[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x64[i])) + x64[i]*log(p[i]) + (samSize-x64[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    }
    
    as.numeric(-log.lik)
  }
  
  M.result2[[i]] <- mle(neg.logL2, start=list(p1 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.result2[i,] <- M.result2[[i]]@details$par
  
  M.result4[[i]] <- mle(neg.logL4, start=list(p1 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.result4[i,] <- M.result4[[i]]@details$par
  
  M.result64[[i]] <- mle(neg.logL64, start=list(p1 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.result64[i,] <- M.result64[[i]]@details$par
  
  print(i) #to check number of runs
  
}

#####################################################################
## Measuring performance
#####################################################################

# Calculating the average calibrated parameteres
MLE.avgPar2 <- round(mean(MLE.result2), 3)
MLE.avgPar4 <- round(mean(MLE.result4), 3)
MLE.avgPar64 <- round(mean(MLE.result64), 3)

print("Average parameter estimate of gamma:")
print(paste0('Model with 2 target features: gamma = ', MLE.avgPar2))
print(paste0('Model with 4 target features: gamma = ', MLE.avgPar4))
print(paste0('Model with 64 target features: gamma = ', MLE.avgPar64))


#Calculating the Bias
MLE.gBias2 <- MLE.avgPar2 - trueGamma
MLE.gBias4 <- MLE.avgPar4 - trueGamma
MLE.gBias64 <- MLE.avgPar64 - trueGamma


print("The Percentage Bias of the parameter estimate gamma:")
print(paste0('Bias for Model with 2 target features: gamma Bias = ', MLE.gBias2/trueGamma*100, '%'))
print(paste0('Bias for Model with 4 target features: gamma Bias = ', MLE.gBias4/trueGamma*100, '%'))
print(paste0('Bias for Model with 64 target features: gamma Bias = ', MLE.gBias64/trueGamma*100, '%'))


#Calculating the accuracy using the Root Mean Square Error
MLE.gAccu2 <- sqrt((sum((MLE.result2 - trueGamma[1])^2)/calibModelRuns)) 
MLE.gAccu4 <- sqrt((sum((MLE.result4 - trueGamma[1])^2)/calibModelRuns)) 
MLE.gAccu64 <- sqrt((sum((MLE.result64 - trueGamma[1])^2)/calibModelRuns)) 

print("The accuracy of the parameter estimate gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: gamma = ',  round(MLE.gAccu2, 3)))
print(paste0('RMSE for Model with 4 target features: gamma = ',  round(MLE.gAccu4, 3)))
print(paste0('RMSE for Model with 64 target features: gamma = ',  round(MLE.gAccu64, 3)))


#Calculating the coverage using confidence intervals

ML.standError2 <- matrix(c(0), calibModelRuns, 1)
ML.standError4 <- matrix(c(0), calibModelRuns, 1)
ML.standError64 <- matrix(c(0), calibModelRuns, 1)

for(i in 1:calibModelRuns){
  ML.standError2[i,] <- sqrt(abs(diag(solve(M.result2[[i]]@details$hessian))))
  ML.standError4[i,] <- sqrt(abs(diag(solve(M.result4[[i]]@details$hessian))))
  ML.standError64[i,] <- sqrt(abs(diag(solve(M.result64[[i]]@details$hessian))))
}

# Now confidence intervals of each of the parameter estimates
M.CI_2_g <- matrix(c(0, 0), calibModelRuns, 2)
M.CI_4_g <- matrix(c(0, 0), calibModelRuns, 2)
M.CI_64_g <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){
  M.CI_2_g[i,1] <- MLE.result2[i] - 1.96*ML.standError2[i]
  M.CI_2_g[i,2] <- MLE.result2[i] + 1.96*ML.standError2[i]
  
  M.CI_4_g[i,1] <- MLE.result4[i] - 1.96*ML.standError4[i]
  M.CI_4_g[i,2] <- MLE.result4[i] + 1.96*ML.standError4[i]
  
  M.CI_64_g[i,1] <- MLE.result64[i] - 1.96*ML.standError64[i]
  M.CI_64_g[i,2] <- MLE.result64[i] + 1.96*ML.standError64[i]
  
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
MLE2_gcov <- sum((trueGamma >= M.CI_2_g[,1]) == TRUE & (trueGamma <= M.CI_2_g[,2]) == TRUE)/calibModelRuns * 100
MLE4_gcov <- sum((trueGamma >= M.CI_4_g[,1]) == TRUE & (trueGamma <= M.CI_4_g[,2]) == TRUE)/calibModelRuns * 100
MLE64_gcov <- sum((trueGamma >= M.CI_64_g[,1]) == TRUE & (trueGamma <= M.CI_64_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of the parameter estimate gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: gamma = ',  MLE2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: gamma = ',  MLE4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: gamma = ',  MLE64_gcov, '%'))




