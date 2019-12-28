## Wynand van Staden
## Rejection ABC simulation and calibration
## Copyright 2019

library(EasyABC)
library(SimInf)

#SIR model for ABC_rejection function with 2 target features
abcsirModel2 <- function(betaGamma, N = 10000, inf = 0.1, sampleSize = 1000){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= betaGamma[1], gamma=betaGamma[2])
  result <- run(model)
  
  
  individuals50 <- c(rep("S", result@U[1, 50]), rep("I", result@U[2, 50]), rep("R", result@U[3, 50]))
  individuals65 <- c(rep("S", result@U[1, 65]), rep("I", result@U[2, 65]), rep("R", result@U[3, 65]))
  
  
  sampleIPop <- c(summary(as.factor(sample(individuals50, size = sampleSize))), summary(as.factor(sample(individuals65, size = sampleSize)))) 
  
  pop <- sampleIPop[names(sampleIPop) == "I"]
  
  for(i in 1:2){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
}


#SIR Model for ABC_rejection function with 4 target features
abcsirModel4 <- function(betaGamma, N = 10000, inf = 0.1, sampleSize = 1000){
  library(SimInf)
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= betaGamma[1], gamma=betaGamma[2])
  result <- run(model)
  
  individuals30 <- c(rep("S", result@U[1, 30]), rep("I", result@U[2, 30]), rep("R", result@U[3, 30]))
  individuals45 <- c(rep("S", result@U[1, 45]), rep("I", result@U[2, 45]), rep("R", result@U[3, 45]))
  individuals60 <- c(rep("S", result@U[1, 60]), rep("I", result@U[2, 60]), rep("R", result@U[3, 60]))
  individuals75 <- c(rep("S", result@U[1, 75]), rep("I", result@U[2, 75]), rep("R", result@U[3, 75]))
  
  
  sampleIPop <- c(summary(as.factor(sample(individuals30, size = sampleSize))), summary(as.factor(sample(individuals45, size = sampleSize))), summary(as.factor(sample(individuals60, size = sampleSize))), summary(as.factor(sample(individuals75, size = sampleSize)))) 
  
  pop <- sampleIPop[names(sampleIPop) == "I"]
  
  
  for(i in 1:4){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }

  return(pop/sampleSize)

}


#SIR Model for ABC_rejection function with 64 target features
abcsirModel64 <- function(betaGamma, N = 10000, inf = 0.1, sampleSize = 1000){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= betaGamma[1], gamma=betaGamma[2])
  result <- run(model)
  
  individuals1.64 <- list()
  
  for(i in 1:64){
    
    individuals1.64[[i]] <- c(rep("S", result@U[1, i]), rep("I", result@U[2, i]), rep("R", result@U[3, i]))
    
  }
  
  sampleIPop <- c()
  for(i in 1:64){
    
    sampleIPop <- c(sampleIPop, c(summary(as.factor(sample(individuals1.64[[i]], size = sampleSize)))))
    
  }
  
  pop <- sampleIPop[names(sampleIPop) == "I"]
  
  
  for(i in 1:64){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
}


abcModelRun <- function(simRep, betaGamma, tol){

#Creating the observed data
trueSIRprev2 <- abcsirModel2(betaGamma)
trueSIRprev4 <- abcsirModel4(betaGamma)
trueSIRprev64 <- abcsirModel64(betaGamma)
  
  
#setting tolerance value
tolp <- 1

#Running the abc-r method
ABC_rej2 <- ABC_rejection(model = abcsirModel2, prior=list(c("unif",0.01,0.5),c("unif",0.01,0.1)), 
                        summary_stat_target= trueSIRprev2, nb_simul=simRep,
                        tol=tolp, progress_bar = T)

ABC_rej4 <- ABC_rejection(model = abcsirModel4, prior=list(c("unif",0.01,0.5),c("unif",0.01,0.1)), 
                          summary_stat_target= trueSIRprev4, nb_simul=simRep,
                          tol=tolp, progress_bar = T)

ABC_rej64 <- ABC_rejection(model = abcsirModel64, prior=list(c("unif",0.01,0.5),c("unif",0.01,0.1)), 
                          summary_stat_target= trueSIRprev64, nb_simul=simRep,
                          tol=tolp, progress_bar = T)




abc2_0.01 <- abc(target = c(trueSIRprev2),
               param = ABC_rej2$param,
               sumstat = ABC_rej2$stats,
               tol = tol,
               method="rejection") 

abc4_0.01 <- abc(target = c(trueSIRprev4),
                 param = ABC_rej4$param,
                 sumstat = ABC_rej4$stats,
                 tol = tol,
                 method="rejection") 

abc64_0.01 <- abc(target = c(trueSIRprev64),
                 param = ABC_rej64$param,
                 sumstat = ABC_rej64$stats,
                 tol = tol,
                 method="rejection") 


return(list(abc2_0.01, abc4_0.01, abc64_0.01))

}

#####################################################################
## Model Calibration
#####################################################################

modelRuns <- 1000
betaGamma <- c(0.2, 0.02)
simRep <- 1000
tol <- 0.1

ABC <- list()
Abc2 <- list()
Abc4 <- list()
Abc64 <- list()

# This runs the calibration method modelRuns times for ModelRuns amount of posterior distributions
for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  ABC[[i]] <- abcModelRun(simRep, betaGamma, tol)
  
  Abc2[i] <- ABC[[i]][1]
  Abc4[i] <- ABC[[i]][2]
  Abc64[i] <- ABC[[i]][3]
}


#Now to calculate the median value of the posterior distributions
Abc.median.2 <- matrix(c(0, 0), modelRuns, 2)
Abc.median.4 <- matrix(c(0, 0), modelRuns, 2)
Abc.median.64 <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  Abc.median.2[i,] <- c(median(Abc2[[i]]$unadj.values[,1]), median(Abc2[[i]]$unadj.values[,2]))
  Abc.median.4[i,] <- c(median(Abc4[[i]]$unadj.values[,1]), median(Abc4[[i]]$unadj.values[,2]))
  Abc.median.64[i,] <- c(median(Abc64[[i]]$unadj.values[,1]), median(Abc64[[i]]$unadj.values[,2]))
  
}


#####################################################################
## Measuring performance
#####################################################################


#Calculating the mean of the medians found
avgAbc2 <- c(round(mean(Abc.median.2[,1]), 3), round(mean(Abc.median.2[,2]), 3))
avgAbc4 <- c(round(mean(Abc.median.4[,1]), 3), round(mean(Abc.median.4[,2]), 3))
avgAbc64 <- c(round(mean(Abc.median.64[,1]), 3), round(mean(Abc.median.64[,2]), 3))


print("The mean parameter estiamtes found of the median from the posterior distributions by using the rejection ABC method")
print(paste0("Model with 2 target statistics: beta= ", avgAbc2[1], ", gamma = ", avgAbc2[2]))
print(paste0("Model with 4 target statistics: beta= ", avgAbc4[1], ", gamma = ", avgAbc4[2]))
print(paste0("Model with 64 target statistics: beta= ", avgAbc64[1], ", gamma = ", avgAbc64[2]))

#Calculating the Bias
ABC.bgBias2 <- avgAbc2 - betaGamma
ABC.bgBias4 <- avgAbc4 - betaGamma
ABC.bgBias64 <- avgAbc64 - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 target features: beta Bias = ', ABC.bgBias2[1]/betaGamma[1] *100, '%, gamma Bias = ',
             ABC.bgBias2[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 4 target features: beta Bias = ', ABC.bgBias4[1]/betaGamma[1] *100, '%, gamma Bias = ',
             ABC.bgBias4[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 64 target features: beta Bias = ', ABC.bgBias64[1]/betaGamma[1] *100, '%, gamma Bias = ',
             ABC.bgBias64[2]/betaGamma[2] *100, '%' ))



#Calculating the accuracy using the Root Mean Square Error
ABC.bgAccu2 <- c(sqrt((sum((Abc.median.2[, 1] - betaGamma[1])^2)/modelRuns)), 
                 sqrt((sum((Abc.median.2[, 2] - betaGamma[2])^2)/modelRuns)))

ABC.bgAccu4 <- c(sqrt((sum((Abc.median.4[, 1] - betaGamma[1])^2)/modelRuns)), 
                 sqrt((sum((Abc.median.4[, 2] - betaGamma[2])^2)/modelRuns)))

ABC.bgAccu64 <- c(sqrt((sum((Abc.median.64[, 1] - betaGamma[1])^2)/modelRuns)), 
                  sqrt((sum((Abc.median.64[, 2] - betaGamma[2])^2)/modelRuns)))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: beta = ', round(ABC.bgAccu2[1], 3), ' gamma = ',  round(ABC.bgAccu2[2], 3)))
print(paste0('RMSE for Model with 4 target features: beta = ', round(ABC.bgAccu4[1], 3), ' gamma = ',  round(ABC.bgAccu4[2], 3)))
print(paste0('RMSE for Model with 64 target features: beta = ', round(ABC.bgAccu64[1], 3), ' gamma = ',  round(ABC.bgAccu64[2], 3)))


## Calculating credible intervals
posteriorSize <- simRep*tol

ABC_2.5 <- posteriorSize * 0.025
ABC_97.5 <- posteriorSize * 0.975

ABC.CI_2_b <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_2_g <- matrix(c(0, 0), modelRuns, 2)

ABC.CI_4_b <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_4_g <- matrix(c(0, 0), modelRuns, 2)

ABC.CI_64_b <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_64_g <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  ABC.CI_2_b[i,] <- c(sort(Abc2[[i]]$unadj.values[,1])[ABC_2.5], sort(Abc2[[i]]$unadj.values[,1])[ABC_97.5])
  ABC.CI_2_g[i,] <- c(sort(Abc2[[i]]$unadj.values[,2])[ABC_2.5], sort(Abc2[[i]]$unadj.values[,2])[ABC_97.5])
  
  ABC.CI_4_b[i,] <- c(sort(Abc4[[i]]$unadj.values[,1])[ABC_2.5], sort(Abc4[[i]]$unadj.values[,1])[ABC_97.5])
  ABC.CI_4_g[i,] <- c(sort(Abc4[[i]]$unadj.values[,2])[ABC_2.5], sort(Abc4[[i]]$unadj.values[,2])[ABC_97.5])
  
  ABC.CI_64_b[i,] <- c(sort(Abc64[[i]]$unadj.values[,1])[ABC_2.5], sort(Abc64[[i]]$unadj.values[,1])[ABC_97.5])
  ABC.CI_64_g[i,] <- c(sort(Abc64[[i]]$unadj.values[,2])[ABC_2.5], sort(Abc64[[i]]$unadj.values[,2])[ABC_97.5])

  }


# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
ABC.2_bcov <- sum((betaGamma[1] >= ABC.CI_2_b[,1]) == TRUE & (betaGamma[1] <= ABC.CI_2_b[,2]) == TRUE)/modelRuns * 100
ABC.2_gcov <- sum((betaGamma[2] >= ABC.CI_2_g[,1]) == TRUE & (betaGamma[2] <= ABC.CI_2_g[,2]) == TRUE)/modelRuns * 100

ABC.4_bcov <- sum((betaGamma[1] >= ABC.CI_4_b[,1]) == TRUE & (betaGamma[1] <= ABC.CI_4_b[,2]) == TRUE)/modelRuns * 100
ABC.4_gcov <- sum((betaGamma[2] >= ABC.CI_4_g[,1]) == TRUE & (betaGamma[2] <= ABC.CI_4_g[,2]) == TRUE)/modelRuns * 100

ABC.64_bcov <- sum((betaGamma[1] >= ABC.CI_64_b[,1]) == TRUE & (betaGamma[1] <= ABC.CI_64_b[,2]) == TRUE)/modelRuns * 100
ABC.64_gcov <- sum((betaGamma[2] >= ABC.CI_64_g[,1]) == TRUE & (betaGamma[2] <= ABC.CI_64_g[,2]) == TRUE)/modelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: beta = ', ABC.2_bcov, '%, gamma = ',  ABC.2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: beta = ', ABC.4_bcov, '%, gamma = ',  ABC.4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: beta = ', ABC.64_bcov, '%, gamma = ',  ABC.64_gcov, '%'))





      