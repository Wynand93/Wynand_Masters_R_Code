## Wynand van Staden
## Least-Squares Claibrtion method using the Peak prevalence data
## Copyright 2019

# sir model returning infectious prevalences at 2 time points + the peak prevalence time point
sirModelPeakPrev <- function(beta, gamma, N = 10000, inf = 0.1, sampleSize = 1000){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= beta, gamma=gamma)
  result <- run(model)
  
  peakIPrev <- which.max(result@U[2,])
  
  individualsPeak <- c(rep("S", result@U[1, peakIPrev]), rep("I", result@U[2, peakIPrev]), rep("R", result@U[3, peakIPrev]))
  individuals50 <- c(rep("S", result@U[1, 50]), rep("I", result@U[2, 50]), rep("R", result@U[3, 50]))
  individuals65 <- c(rep("S", result@U[1, 65]), rep("I", result@U[2, 65]), rep("R", result@U[3, 65]))
  
  
  samplePop <- c(summary(as.factor(sample(individualsPeak, size = sampleSize))), summary(as.factor(sample(individuals50, size = sampleSize))), summary(as.factor(sample(individuals65, size = sampleSize)))) 
  
  pop <- samplePop[names(samplePop) == "I"] #let it return zero instead of numeric(0)
  
  for(i in 1:3){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
}


#####################################################################
## Model Calibration
#####################################################################

betaGamma <- c(0.2, 0.02)
calibModelRuns <- 1000

#Initializing variables to store results
infPopPeak <- matrix(c(0, 0, 0), calibModelRuns, 3)

resultPeak <- list()

#to just store beta and gamma parameters from the reuslt
LS.resultPeak <- matrix(c(0, 0), calibModelRuns, 2)

#This step takes quite a long time
print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {

  infPopPeak[i,] <- sirModelPeakPrev(betaGamma[1], betaGamma[2])
  
  f.optPeak <- function(params) sum((sirModelPeakPrev(params[1],params[2])- infPopPeak[i,])^2 )
  
  resultPeak[[i]] <- optim(c(runif(1, 0.01, 0.5), runif(1, 0.01, 0.1)), f.optPeak,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000), hessian = T)
  LS.resultPeak[i, ] <- resultPeak[[i]]$par 
  
  print(i) #to check number of runs
  
}


#####################################################################
## Measuring performance
#####################################################################


# Calculating the average calibrated parameteres
avgParPeak <- c(round(mean(LS.resultPeak[,1]), 3) , round(mean(LS.resultPeak[, 2]), 3))

print("Average parameter estimates of beta and gamma:")
print(paste0('Model with 2 + Peak target features: beta = ', avgParPeak[1] , ', gamma = ', avgParPeak[2]))

#Calculating the Bias
bgBiasPeak <- avgParPeak - betaGamma

print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 + Peak target features: beta Bias = ', bgBiasPeak[1]/betaGamma[1] *100, '%, gamma Bias = ',
             bgBiasPeak[2]/betaGamma[2] *100, '%' ))

#Calculating the accuracy using the Root Mean Square Error
bgAccuPeak <- c(sqrt((sum((LS.resultPeak[, 1] - betaGamma[1])^2)/calibModelRuns)), 
             sqrt((sum((LS.resultPeak[, 2] - betaGamma[2])^2)/calibModelRuns)))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 + Peak target features: beta = ', round(bgAccuPeak[1], 3), ' gamma = ',  round(bgAccuPeak[2], 3)))



#Calculating the coverage using confidence intervals
standErrorPeak <- matrix(c(0, 0), calibModelRuns, 2)

for(i in 1:calibModelRuns){
  standErrorPeak[i,] <- sqrt(abs(diag(solve(resultPeak[[i]]$hessian))))
}


# Now confidence intervals of each of the parameter estimates
CI_Peak_b <- matrix(c(0, 0), calibModelRuns, 2)
CI_Peak_g <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){
  CI_Peak_b[i,1] <- LS.resultPeak[i,1] - 1.96*standErrorPeak[i,1]
  CI_Peak_b[i,2] <- LS.resultPeak[i,1] + 1.96*standErrorPeak[i,1]
  
  CI_Peak_g[i,1] <- LS.resultPeak[i,2] - 1.96*standErrorPeak[i,2]
  CI_Peak_g[i,2] <- LS.resultPeak[i,2] + 1.96*standErrorPeak[i,2]

}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
LSPeak_bcov <- sum((betaGamma[1] >= CI_Peak_b[,1]) == TRUE & (betaGamma[1] <= CI_Peak_b[,2]) == TRUE)/calibModelRuns * 100
LSPeak_gcov <- sum((betaGamma[2] >= CI_Peak_g[,1]) == TRUE & (betaGamma[2] <= CI_Peak_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 + Peak target features: beta = ', LSPeak_bcov, '%, gamma = ',  LSPeak_gcov, '%'))

