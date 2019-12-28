## Wynand van Staden
## Rejection ABC calibration method with Peak prev
## Copyright 2019

library(EasyABC)
library(SimInf)

#SIR model for ABC_rejection function with 2  + Peak target features
abcsirModelPeak <- function(betaGamma, N = 10000, inf = 0.1, sampleSize = 1000){
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= betaGamma[1], gamma=betaGamma[2])
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


abcModelRunPeak <- function(simRep, betaGamma, tol){
  
  #Creating the observed data
  trueSIRprevPeak <- abcsirModelPeak(betaGamma)
  
  tolp <- 1   
  
  ABC_rejPeak <- ABC_rejection(model = abcsirModelPeak, prior=list(c("unif",0.01,0.5),c("unif",0.01,0.1)), 
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
betaGamma <- c(0.2, 0.02)
simRep <- 1000
tol <- 0.1

ABC <- list()
AbcPeak <- list()


for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  AbcPeak[[i]] <- abcModelRunPeak(simRep, betaGamma, tol)
  
}

Abc.median.Peak <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  Abc.median.Peak[i,] <- c(median(AbcPeak[[i]]$unadj.values[,1]), median(AbcPeak[[i]]$unadj.values[,2]))
  
}

#####################################################################
## Measuring performance
#####################################################################


avgAbcPeak <- c(round(mean(Abc.median.Peak[,1]), 3), round(mean(Abc.median.Peak[,2]), 3))

print("The parameter estiamtes found by using the rejection ABC method")
print(paste0("Model with 2 + Peak target statistics: beta= ", avgAbcPeak[1], ", gamma = ", avgAbcPeak[2]))

#Calculating the Bias
ABC.bgBiasPeak <- avgAbcPeak - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 + Peak target features: beta Bias = ', ABC.bgBiasPeak[1]/betaGamma[1] *100, '%, gamma Bias = ',
             ABC.bgBiasPeak[2]/betaGamma[2] *100, '%' ))


#Calculating the accuracy using the Root Mean Square Error
ABC.bgAccuPeak <- c(sqrt((sum((Abc.mode.Peak[, 1] - betaGamma[1])^2)/modelRuns)), 
                 sqrt((sum((Abc.mode.Peak[, 2] - betaGamma[2])^2)/modelRuns)))


print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 + Peak target features: beta = ', round(ABC.bgAccuPeak[1], 3), ' gamma = ',  round(ABC.bgAccuPeak[2], 3)))


## Calculating credible intervals
posteriorSize <- simRep*tol

ABC_2.5 <- posteriorSize * 0.025
ABC_97.5 <- posteriorSize * 0.975

ABC.CI_Peak_b <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_Peak_g <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  ABC.CI_Peak_b[i,] <- c(sort(AbcPeak[[i]]$unadj.values[,1])[ABC_2.5], sort(AbcPeak[[i]]$unadj.values[,1])[ABC_97.5])
  ABC.CI_Peak_g[i,] <- c(sort(AbcPeak[[i]]$unadj.values[,2])[ABC_2.5], sort(AbcPeak[[i]]$unadj.values[,2])[ABC_97.5])
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
ABC.Peak_bcov <- sum((betaGamma[1] >= ABC.CI_Peak_b[,1]) == TRUE & (betaGamma[1] <= ABC.CI_Peak_b[,2]) == TRUE)/modelRuns * 100
ABC.Peak_gcov <- sum((betaGamma[2] >= ABC.CI_Peak_g[,1]) == TRUE & (betaGamma[2] <= ABC.CI_Peak_g[,2]) == TRUE)/modelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 + Peak target features: beta = ', ABC.Peak_bcov, '%, gamma = ',  ABC.Peak_gcov, '%'))




