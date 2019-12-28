## Wynand van Staden
## Maximum Likelihood estimation Calibration method with peak Prevalence
## Copyright 2019

## Please note: This file makes use of the sirModelPeakPrev() function in the One parameter 1 Peak Prev with LS.r file

library(stats4)
library(dplyr)
library(gmp)

#####################################################################
## Model Calibration
#####################################################################


trueGamma <- 0.02
calibModelRuns <- 1000

#Sample size is needed for the likelihood calculations
samSize <- 1000

MLEinfPopPeak <- matrix(c(0, 0, 0), calibModelRuns, 3)
M.resultPeak <- list()

#to gamma parameters from the reuslt
MLE.resultPeakPar <- c()


print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {
  
  MLEinfPopPeak[i,] <- sirModelPeakPrev(trueGamma)
  
  #calculating the likelihood of the explored parameters
  neg.logLPeak <- function(p1)  {
    p <- sirModelPeakPrev(p1)
    p <- ifelse(p==0, 0.0001, p)
    x2 <- MLEinfPopPeak[i, ] * samSize
    
    log.lik <- 0
    for(i in 1:length(p)){
      log.likCalculate <- log(chooseZ(samSize, x2[i])) + x2[i]*log(p[i]) + (samSize-x2[i])*log(1-p[i])
      log.lik <- log.lik + log.likCalculate
    }
    
    as.numeric(-log.lik)
  }
  
  M.resultPeak[[i]] <- mle(neg.logLPeak, start=list(p1 = runif(1, 0.01, 0.1)),   method = "Nelder-Mead", lower = 0.01, upper = 0.1, control=list(maxit=1000))
  MLE.resultPeakPar[i] <- M.resultPeak[[i]]@details$par
 
  print(i) #to check number of runs
  
}

#####################################################################
## Measuring performance
#####################################################################


# Calculating the average calibrated parameteres
MLE.avgParPeak <- round(mean(MLE.resultPeakPar), 3)

print("Average parameter estimate of gamma:")
print(paste0('Model with 2 + Peak target features: gamma = ', MLE.avgParPeak))

#Calculating the Bias
MLE.gBiasPeak <- MLE.avgParPeak - trueGamma

print("The Percentage Bias of the parameter estimate of gamma:")
print(paste0('Bias for Model with 2 target features + Peak prev: gamma Bias = ', MLE.gBiasPeak/trueGamma *100, '%'))



#Calculating the accuracy using the Root Mean Square Error
MLE.gAccuPeak <- sqrt((sum((MLE.resultPeakPar - trueGamma)^2)/calibModelRuns))

print("The accuracy of the parameter estimates of gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features + Peak prev: gamma = ',  round(MLE.gAccuPeak, 3)))


#Calculating the coverage using confidence intervals
ML.standErrorPeak <- c()

for(i in 1:calibModelRuns){
  ML.standErrorPeak[i] <- sqrt(abs(diag(solve(M.resultPeak[[i]]@details$hessian))))
}


# Now confidence intervals of each of the parameter estimates
M.CI_Peak_g <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){

  M.CI_Peak_g[i,] <- c(MLE.resultPeakPar[i] - 1.96*ML.standErrorPeak[i], MLE.resultPeakPar[i] + 1.96*ML.standErrorPeak[i])
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
MLEPeak_gcov <- sum((trueGamma >= M.CI_Peak_g[,1]) == TRUE & (trueGamma <= M.CI_Peak_g[,2]) == TRUE)/calibModelRuns * 100


print("The coverage of the parameter estimates gamma given the CI's:")
print(paste0('Coverage for Model with 2 + Peak target features: gamma = ',  MLEPeak_gcov, '%'))



