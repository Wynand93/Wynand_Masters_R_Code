## Wynand van Staden
## Rejection ABC calibration method with Peak prev
## Copyright 2019

library(EasyABC)
library(SimInf)

#SIR model for ABC_rejection function with 2  + Peak target features
abcsirModelPeak <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = 1000){
  
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


#Creating a function to use the ABC-r method with the abcSirModelPeak model
abcModelRunPeak <- function(simRep, trueGamma, tol){
  
  #Creating the observed data
  trueSIRprevPeak <- abcsirModelPeak(trueGamma)
  
  tolp <- 1   
  
  ABC_rejPeak <- ABC_rejection(model = abcsirModelPeak, prior=list(c("unif",0.01,0.1)), 
                            summary_stat_target= trueSIRprevPeak, nb_simul=simRep,
                            tol=tolp, progress_bar = T)
 
 
  abcPeak_0.01 <- abc(target = c(trueSIRprevPeak),
                   param = ABC_rejPeak$param,
                   sumstat = ABC_rejPeak$stats,
                   tol = tol,
                   method="rejection") 
  
  return(abcPeak_0.01)
  
}

#####################################################################
## Model Calibration
#####################################################################

modelRuns <- 1000
trueGamma <- 0.02
simRep <- 1000
tol <- 0.1

ABC <- list()
AbcPeak <- list()

#Now to run the calibration
for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  AbcPeak[[i]] <- abcModelRunPeak(simRep, trueGamma, tol)
  
}

#Calculating the median of each of the posterior distributions
Abc.g.median.Peak <- c()

for(i in 1:modelRuns){
  
  Abc.g.median.Peak[i] <-median(AbcPeak[[i]]$unadj.values)
  
}

#####################################################################
## Measuring performance
#####################################################################

avgAbcPeak <- round(mean(Abc.g.median.Peak), 3)

print("The parameter estiamtes found by using the rejection ABC method")
print(paste0("Model with 2 + Peak target statistics: gamma = ", avgAbcPeak))

#Calculating the Bias
ABC.gBiasPeak <- avgAbcPeak - trueGamma


print("The Percentage Bias of the parameter estimate of gamma:")
print(paste0('Bias for Model with 2 + Peak target features: gamma Bias = ',ABC.gBiasPeak/trueGamma * 100, '%' ))

#Calculating the accuracy using the Root Mean Square Error
ABC.gAccuPeak <- sqrt((sum((Abc.g.median.Peak - trueGamma)^2)/modelRuns))

print("The accuracy of the parameter estimates of gamma using RMSE:")
print(paste0('RMSE for Model with 2 + Peak target features: gamma = ',  round(ABC.gAccuPeak, 3)))

## Calculating credible intervals
posteriorSizePeak <- simRep*tol

ABC_2.5 <-  posteriorSizePeak * 0.025
ABC_97.5 <- posteriorSizePeak * 0.975

ABC.CI_Peak_b <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_Peak_g <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  ABC.CI_Peak_g[i,] <- c(sort(AbcPeak[[i]]$unadj.values)[ABC_2.5], sort(AbcPeak[[i]]$unadj.values)[ABC_97.5])
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
ABC.Peak_gcov <- sum((trueGamma >= ABC.CI_Peak_g[,1]) == TRUE & (trueGamma <= ABC.CI_Peak_g[,2]) == TRUE)/modelRuns * 100


print("The coverage of the parameter estimates of gamma given the CI's:")
print(paste0('Coverage for Model with 2 + Peak target features: gamma = ',  ABC.Peak_gcov, '%'))




