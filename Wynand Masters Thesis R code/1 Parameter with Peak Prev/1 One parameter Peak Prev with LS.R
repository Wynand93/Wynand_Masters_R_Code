## Wynand van Staden
## Least-Squares Claibrtion method using the Peak prevalence data
## Copyright 2019


samSize <- 1000

# sir model returning infectious prevalences at 2 time points + the peak prevalence time point
sirModelPeakPrev <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
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

trueGamma <- 0.02
calibModelRuns <- 1000

infPopPeak <- matrix(c(0, 0, 0), calibModelRuns, 3)

resultPeak <- list()

#to just store beta and gamma parameters from the reuslt
LS.resultPeak <- c()

#This step takes quite a long time
print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {

  infPopPeak[i,] <- sirModelPeakPrev(trueGamma)
  
  f.optPeak <- function(params) sum((sirModelPeakPrev(params)- infPopPeak[i,])^2 )
  
  resultPeak[[i]] <- optim(c(runif(1, 0.01, 0.1)), f.optPeak,  method = "Nelder-Mead",lower = 0.01, upper = 0.1, control=list(maxit=1000), hessian = T)
  LS.resultPeak[i] <- resultPeak[[i]]$par 
  
  print(i) #to check number of runs
  
}

#####################################################################
## Measuring performance
#####################################################################

# Calculating the average calibrated parameteres
avgParPeak <- round(mean(LS.resultPeak), 3)

print("Average of the parameter estimates of gamma:")
print(paste0('Model with 2 + Peak target features: gamma = ', avgParPeak))

#Calculating the Bias
gBiasPeak <- avgParPeak - trueGamma

print("The Percentage Bias of the parameter estimates of gamma:")
print(paste0('Bias for Model with 2 + Peak target features: gamma Bias = ', gBiasPeak/trueGamma *100, '%' ))


#Calculating the accuracy using the Root Mean Square Error
gAccuPeak <- sqrt((sum((LS.resultPeak - trueGamma)^2)/calibModelRuns))

print("The accuracy of the parameter estimates for gamma using RMSE:")
print(paste0('RMSE for Model with 2 + Peak target features: gamma = ',  round(gAccuPeak, 3)))


#Calculating the coverage using confidence intervals

standErrorPeak <- c()

for(i in 1:calibModelRuns){
  standErrorPeak[i] <- sqrt(abs(diag(solve(resultPeak[[i]]$hessian))))
}

# Now confidence intervals of the parameter estimates
CI_Peak_g <- matrix(c(0, 0), calibModelRuns, 2)

for(i in 1:calibModelRuns){
  
  CI_Peak_g[i,] <- c(LS.resultPeak[i] - 1.96*standErrorPeak[i], LS.resultPeak[i] + 1.96*standErrorPeak[i])

}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
LSPeak_gcov <- sum((trueGamma >= CI_Peak_g[,1]) == TRUE & (trueGamma <= CI_Peak_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of each parameter estimates of gamma given the CI's:")
print(paste0('Coverage for Model with 2 + Peak target features: gamma = ',  LSPeak_gcov, '%'))

