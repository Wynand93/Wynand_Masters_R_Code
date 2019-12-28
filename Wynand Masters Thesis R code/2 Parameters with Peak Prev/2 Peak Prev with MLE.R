## Wynand van Staden
## Maximum Likelihood estimation Calibration method with peak Prevalence
## Copyright 2019

## Please note this file makes use of the sirModelPeakPrev() functions in the 1.Peak Prev with LS.r file

library(stats4)
library(dplyr)
library(gmp)

#####################################################################
## Model Calibration
#####################################################################

betaGamma <- c(0.2, 0.02)
calibModelRuns <- 1000
samSize <- 1000

#Initialize variables to store results
MLEinfPopPeak <- matrix(c(0, 0, 0), calibModelRuns, 3)

M.resultPeak <- list()


#to just store beta and gamma parameters from the reuslt
MLE.resultPeak <- matrix(c(0, 0), calibModelRuns, 2)

print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {
  
  MLEinfPopPeak[i,] <- sirModelPeakPrev(betaGamma[1], betaGamma[2])
  
  neg.logLPeak <- function(p1, p2)  {
    p <- sirModelPeakPrev(p1, p2)
    p <- ifelse(p==0, 0.0001, p)
    x2 <- MLEinfPopPeak[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x2[i])) + x2[i]*log(p[i]) + (samSize-x2[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    }
    
    as.numeric(-log.lik)
  }
  
  M.resultPeak[[i]] <- mle(neg.logLPeak, start=list(p1 = runif(1, 0.01, 0.5), p2 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000))
  MLE.resultPeak[i,] <- M.resultPeak[[i]]@details$par
 
  print(i) #to check number of runs
  
}

#####################################################################
## Measuring performance
#####################################################################


# Calculating the average calibrated parameteres
MLE.avgParPeak <- c(round(mean(MLE.resultPeak[,1]), 3) , round(mean(MLE.resultPeak[, 2]), 3))

print("Average parameter estimates of beta and gamma:")
print(paste0('Model with 2 + Peak target features: beta = ', MLE.avgParPeak[1] , ', gamma = ', MLE.avgParPeak[2]))



#Calculating the Bias
MLE.bgBiasPeak <- MLE.avgParPeak - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 target features: beta Bias = ', MLE.bgBiasPeak[1]/betaGamma[1] *100, '%, gamma Bias = ',
             MLE.bgBiasPeak[2]/betaGamma[2] *100, '%' ))



#Calculating the accuracy using the Root Mean Square Error
MLE.bgAccuPeak <- c(sqrt((sum((MLE.resultPeak[, 1] - betaGamma[1])^2)/calibModelRuns)), 
                 sqrt((sum((MLE.resultPeak[, 2] - betaGamma[2])^2)/calibModelRuns)))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: beta = ', round(MLE.bgAccuPeak[1], 3), ' gamma = ',  round(MLE.bgAccuPeak[2], 3)))


#Calculating the coverage using confidence intervals
ML.standErrorPeak <- matrix(c(0, 0), calibModelRuns, 2)

for(i in 1:calibModelRuns){
  ML.standErrorPeak[i,] <- sqrt(abs(diag(solve(M.resultPeak[[i]]@details$hessian))))
}

# Now confidence intervals of each of the parameter estimates
M.CI_Peak_b <- matrix(c(0, 0), calibModelRuns, 2)
M.CI_Peak_g <- matrix(c(0, 0), calibModelRuns, 2)

for(i in 1:calibModelRuns){
  
  M.CI_Peak_b[i,1] <- MLE.resultPeak[i,1] - 1.96*ML.standErrorPeak[i,1]
  M.CI_Peak_b[i,2] <- MLE.resultPeak[i,1] + 1.96*ML.standErrorPeak[i,1]
  
  M.CI_Peak_g[i,1] <- MLE.resultPeak[i,2] - 1.96*ML.standErrorPeak[i,2]
  M.CI_Peak_g[i,2] <- MLE.resultPeak[i,2] + 1.96*ML.standErrorPeak[i,2]
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
MLEPeak_bcov <- sum((betaGamma[1] >= M.CI_Peak_b[,1]) == TRUE & (betaGamma[1] <= M.CI_Peak_b[,2]) == TRUE)/calibModelRuns * 100
MLEPeak_gcov <- sum((betaGamma[2] >= M.CI_Peak_g[,1]) == TRUE & (betaGamma[2] <= M.CI_Peak_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 + Peak target features: beta = ', MLEPeak_bcov, '%, gamma = ',  MLEPeak_gcov, '%'))



