## Wynand van Staden
## Model - Data simulation & Maximum Likelihood estimation
## Copyright 2019


## Please note this file makes use of the sirModel() functions in the 1. Data simulate for Least Squares.r file

library(stats4)
library(dplyr)
library(gmp)

#####################################################################
## Model Calibration
#####################################################################

betaGamma <- c(0.2, 0.02)
calibModelRuns <- 1000
samSize <- 10000

#Initializing variable to store results
MLEinfPop2 <- matrix(c(0, 0), calibModelRuns, 2)
MLEinfPop4 <- matrix(replicate(4, 0), calibModelRuns, 4)
MLEinfPop64 <- matrix(replicate(64, 0), calibModelRuns, 64)

M.result2 <- list()
M.result4 <- list()
M.result64 <- list()


#to just store beta and gamma parameters from the reuslt
MLE.result2 <- matrix(c(0, 0), calibModelRuns, 2)
MLE.result4 <- matrix(c(0, 0), calibModelRuns, 2)
MLE.result64 <- matrix(c(0, 0), calibModelRuns, 2)


print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {
  
  MLEinfPop2[i,] <- sirModel2(betaGamma[1], betaGamma[2])
  MLEinfPop4[i,] <- sirModel4(betaGamma[1], betaGamma[2])
  MLEinfPop64[i,] <- sirModel64(betaGamma[1], betaGamma[2])
  
  neg.logL2 <- function(p1, p2)  {
    p <- sirModel2(p1, p2)
    p <- ifelse(p==0, 0.0001, p)
    x2 <- MLEinfPop2[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x2[i])) + x2[i]*log(p[i]) + (samSize-x2[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    }
    
    as.numeric(-log.lik)
  }
  
  samSize = 100
  x2 = 10
  p = 0.09
  
  log(chooseZ(samSize, x2)) + x2*log(p) + (samSize-x2)*log(1-p)
  
  dbinom()
  
  neg.logL4 <- function(p1, p2)  {
    p <- sirModel4(p1,p2)
    p <- ifelse(p==0, 0.0001, p)
    x4 <- MLEinfPop4[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x4[i])) + x4[i]*log(p[i]) + (samSize-x4[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    
    }
    
    as.numeric(-log.lik)
  }
  
  neg.logL64 <- function(p1, p2)  {
    p <- sirModel64(p1,p2)
    p <- ifelse(p==0, 0.0001, p)
    x64 <- MLEinfPop64[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x64[i])) + x64[i]*log(p[i]) + (samSize-x64[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    }
    
    as.numeric(-log.lik)
  }
  
  M.result2[[i]] <- mle(neg.logL2, start=list(p1 = runif(1, 0.01, 0.5), p2 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.result2[i,] <- M.result2[[i]]@details$par
  
  M.result4[[i]] <- mle(neg.logL4, start=list(p1 = runif(1, 0.01, 0.5), p2 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.result4[i,] <- M.result4[[i]]@details$par
  
  M.result64[[i]] <- mle(neg.logL64, start=list(p1 = runif(1, 0.01, 0.5), p2 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.result64[i,] <- M.result64[[i]]@details$par
  
  print(i) #to check number of runs
  
}

#####################################################################
## Measuring performance
#####################################################################

# Calculating the average calibrated parameteres
MLE.avgPar2 <- c(round(mean(MLE.result2[,1]), 3) , round(mean(MLE.result2[, 2]), 3))
MLE.avgPar4 <- c(round(mean(MLE.result4[,1]), 3) , round(mean(MLE.result4[, 2]), 3))
MLE.avgPar64 <- c(round(mean(MLE.result64[,1]), 3) , round(mean(MLE.result64[, 2]), 3))

print("Average parameter estimates of beta and gamma:")
print(paste0('Model with 2 target features: beta = ', MLE.avgPar2[1] , ', gamma = ', MLE.avgPar2[2]))
print(paste0('Model with 4 target features: beta = ', MLE.avgPar4[1] , ', gamma = ', MLE.avgPar4[2]))
print(paste0('Model with 64 target features: beta = ', MLE.avgPar64[1] , ', gamma = ', MLE.avgPar64[2]))



#Calculating the Bias
MLE.bgBias2 <- MLE.avgPar2 - betaGamma
MLE.bgBias4 <- MLE.avgPar4 - betaGamma
MLE.bgBias64 <- MLE.avgPar64 - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 target features: beta Bias = ', MLE.bgBias2[1]/betaGamma[1] *100, '%, gamma Bias = ',
             MLE.bgBias2[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 4 target features: beta Bias = ', MLE.bgBias4[1]/betaGamma[1] *100, '%, gamma Bias = ',
             MLE.bgBias4[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 64 target features: beta Bias = ', MLE.bgBias64[1]/betaGamma[1] *100, '%, gamma Bias = ',
             MLE.bgBias64[2]/betaGamma[2] *100, '%' ))



#Calculating the accuracy using the Root Mean Square Error
MLE.bgAccu2 <- c(sqrt((sum((MLE.result2[, 1] - betaGamma[1])^2)/calibModelRuns)), 
                 sqrt((sum((MLE.result2[, 2] - betaGamma[2])^2)/calibModelRuns)))

MLE.bgAccu4 <- c(sqrt((sum((MLE.result4[, 1] - betaGamma[1])^2)/calibModelRuns)), 
                 sqrt((sum((MLE.result4[, 2] - betaGamma[2])^2)/calibModelRuns)))

MLE.bgAccu64 <- c(sqrt((sum((MLE.result64[, 1] - betaGamma[1])^2)/calibModelRuns)), 
                  sqrt((sum((MLE.result64[, 2] - betaGamma[2])^2)/calibModelRuns)))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: beta = ', round(MLE.bgAccu2[1], 3), ' gamma = ',  round(MLE.bgAccu2[2], 3)))
print(paste0('RMSE for Model with 4 target features: beta = ', round(MLE.bgAccu4[1], 3), ' gamma = ',  round(MLE.bgAccu4[2], 3)))
print(paste0('RMSE for Model with 64 target features: beta = ', round(MLE.bgAccu64[1], 3), ' gamma = ',  round(MLE.bgAccu64[2], 3)))


#Calculating the coverage using confidence intervals
ML.standError2 <- matrix(c(0, 0), calibModelRuns, 2)
ML.standError4 <- matrix(c(0, 0), calibModelRuns, 2)
ML.standError64 <- matrix(c(0, 0), calibModelRuns, 2)



for(i in 1:calibModelRuns){
  ML.standError2[i,] <- sqrt(abs(diag(solve(M.result2[[i]]@details$hessian))))
  ML.standError4[i,] <- sqrt(abs(diag(solve(M.result4[[i]]@details$hessian))))
  ML.standError64[i,] <- sqrt(abs(diag(solve(M.result64[[i]]@details$hessian))))
}


# Now confidence intervals of each of the parameter estimates
M.CI_2_b <- matrix(c(0, 0), calibModelRuns, 2)
M.CI_2_g <- matrix(c(0, 0), calibModelRuns, 2)

M.CI_4_b <- matrix(c(0, 0), calibModelRuns, 2)
M.CI_4_g <- matrix(c(0, 0), calibModelRuns, 2)

M.CI_64_b <- matrix(c(0, 0), calibModelRuns, 2)
M.CI_64_g <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){
  
  M.CI_2_b[i,] <- c(MLE.result2[i,1] - 1.96*ML.standError2[i,1], MLE.result2[i,1] + 1.96*ML.standError2[i,1])
  M.CI_2_g[i,] <- c(MLE.result2[i,2] - 1.96*ML.standError2[i,2], MLE.result2[i,2] + 1.96*ML.standError2[i,2])
  
  M.CI_4_b[i,] <- c(MLE.result4[i,1] - 1.96*ML.standError4[i,1], MLE.result4[i,1] + 1.96*ML.standError4[i,1])
  M.CI_4_g[i,] <- c(MLE.result4[i,2] - 1.96*ML.standError4[i,2], MLE.result4[i,2] + 1.96*ML.standError4[i,2])
  
  M.CI_64_b[i,] <- c(MLE.result64[i,1] - 1.96*ML.standError64[i,1], MLE.result64[i,1] + 1.96*ML.standError64[i,1])
  M.CI_64_g[i,] <- c(MLE.result64[i,2] - 1.96*ML.standError64[i,2], MLE.result64[i,2] + 1.96*ML.standError64[i,2])
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
MLE2_bcov <- sum((betaGamma[1] >= M.CI_2_b[,1]) == TRUE & (betaGamma[1] <= M.CI_2_b[,2]) == TRUE)/calibModelRuns * 100
MLE2_gcov <- sum((betaGamma[2] >= M.CI_2_g[,1]) == TRUE & (betaGamma[2] <= M.CI_2_g[,2]) == TRUE)/calibModelRuns * 100

MLE4_bcov <- sum((betaGamma[1] >= M.CI_4_b[,1]) == TRUE & (betaGamma[1] <= M.CI_4_b[,2]) == TRUE)/calibModelRuns * 100
MLE4_gcov <- sum((betaGamma[2] >= M.CI_4_g[,1]) == TRUE & (betaGamma[2] <= M.CI_4_g[,2]) == TRUE)/calibModelRuns * 100

MLE64_bcov <- sum((betaGamma[1] >= M.CI_64_b[,1]) == TRUE & (betaGamma[1] <= M.CI_64_b[,2]) == TRUE)/calibModelRuns * 100
MLE64_gcov <- sum((betaGamma[2] >= M.CI_64_g[,1]) == TRUE & (betaGamma[2] <= M.CI_64_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: beta = ', MLE2_bcov, '%, gamma = ',  MLE2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: beta = ', MLE4_bcov, '%, gamma = ',  MLE4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: beta = ', MLE64_bcov, '%, gamma = ',  MLE64_gcov, '%'))




