## Wynand van Staden
## Rejection ABC simulation and calibration
## Copyright 2019

library(EasyABC)
library(SimInf)

samSize <- 1000

#SIR model for ABC_rejection function with 2 target features
abcsirModel2 <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta = beta, gamma = gamma)
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
abcsirModel4 <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta = beta, gamma = gamma)
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
abcsirModel64 <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta = beta, gamma = gamma)
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


abcModelRun <- function(simRep, trueGamma, tol){

# Parameters For printing the prior and posterior distributions
simRep = 1000
trueGamma = 0.02
tol = 0.1
  
#Creating the observed data
trueSIRprev2 <- abcsirModel2(trueGamma)
trueSIRprev4 <- abcsirModel4(trueGamma)
trueSIRprev64 <- abcsirModel64(trueGamma)
  
#Accepting all parameters in the intial calibration (tolerance will be set higher in the code to follow)
tolp <- 1

#Running the ABC-r model
ABC_rej2 <- ABC_rejection(model = abcsirModel2, prior=list(c("unif",0.01,0.1)), 
                        summary_stat_target= trueSIRprev2, nb_simul=simRep,
                        tol=tolp, progress_bar = T)

ABC_rej4 <- ABC_rejection(model = abcsirModel4, prior=list(c("unif",0.01,0.1)), 
                          summary_stat_target= trueSIRprev4, nb_simul=simRep,
                          tol=tolp, progress_bar = T)

ABC_rej64 <- ABC_rejection(model = abcsirModel64, prior=list(c("unif",0.01,0.1)), 
                          summary_stat_target= trueSIRprev64, nb_simul=simRep,
                          tol=tolp, progress_bar = T)


#Now finding the posterior distribrution with tol% best parameters
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

###############################################################################################
## Model Calibration
###############################################################################################

#This code runs the ABC-r calibration method
modelRuns <- 1000
trueGamma <- 0.02
simRep <- 1000
tol <- 0.1

ABC <- list()
Abc2 <- list()
Abc4 <- list()
Abc64 <- list()

#Run the method for modelRuns times to get modelRuns amount of posterior distributions
for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  ABC[[i]] <- abcModelRun(simRep, trueGamma, tol)
  
  Abc2[i] <- ABC[[i]][1]
  Abc4[i] <- ABC[[i]][2]
  Abc64[i] <- ABC[[i]][3]
}


#Now to find the median values of the posterior distributions
Abc.g.median.2 <- c()
Abc.g.median.4 <- c()
Abc.g.median.64 <- c()

for(i in 1:modelRuns){
  
  Abc.g.median.2[i] <- median(Abc2[[i]]$unadj.values[,1])
  Abc.g.median.4[i] <- median(Abc4[[i]]$unadj.values[,1])
  Abc.g.median.64[i] <- median(Abc64[[i]]$unadj.values[,1])
  
}

#####################################################################
## Measuring performance
#####################################################################


#Calculating the average of the mdeian values
avgAbc2 <- round(mean(Abc.g.median.2), 3)
avgAbc4 <- round(mean(Abc.g.median.4), 3)
avgAbc64 <- round(mean(Abc.g.median.64), 3)


print("The parameter estiamtes found by using the rejection ABC method")
print(paste0("Model with 2 target statistics: gamma = ", avgAbc2))
print(paste0("Model with 4 target statistics: gamma = ", avgAbc4))
print(paste0("Model with 64 target statistics: gamma = ", avgAbc64))

#Calculating the Bias
ABC.gBias2 <- avgAbc2 - trueGamma
ABC.gBias4 <- avgAbc4 - trueGamma
ABC.gBias64 <- avgAbc64 - trueGamma


print("The Percentage Bias of the parameter estimate gamma:")
print(paste0('Bias for Model with 2 target features: gamma Bias = ', ABC.gBias2/trueGamma *100, '%' ))
print(paste0('Bias for Model with 4 target features: gamma Bias = ', ABC.gBias4/trueGamma *100, '%' ))
print(paste0('Bias for Model with 64 target features: gamma Bias = ', ABC.gBias64/trueGamma *100, '%' ))



#Calculating the accuracy using the Root Mean Square Error
ABC.gAccu2 <- sqrt((sum((Abc.g.median.2 - trueGamma)^2)/modelRuns)) 
ABC.gAccu4 <- sqrt((sum((Abc.g.median.4 - trueGamma)^2)/modelRuns)) 
ABC.gAccu64 <- sqrt((sum((Abc.g.median.64 - trueGamma)^2)/modelRuns)) 

print("The accuracy of the parameter estimate gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: gamma = ',  round(ABC.gAccu2, 3)))
print(paste0('RMSE for Model with 4 target features: gamma = ',  round(ABC.gAccu4, 3)))
print(paste0('RMSE for Model with 64 target features: gamma = ',  round(ABC.gAccu64, 3)))


## Calculating credible intervals
posteriorSize <- simRep*tol

ABC_2.5 <- posteriorSize * 0.025
ABC_97.5 <- posteriorSize * 0.975

ABC.CI_2_g <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_4_g <- matrix(c(0, 0), modelRuns, 2)
ABC.CI_64_g <- matrix(c(0, 0), modelRuns, 2)


for(i in 1:modelRuns){
  
  ABC.CI_2_g[i,] <- c(sort(Abc2[[i]]$unadj.values)[ABC_2.5], sort(Abc2[[i]]$unadj.values)[ABC_97.5])
  ABC.CI_4_g[i,] <- c(sort(Abc4[[i]]$unadj.values)[ABC_2.5], sort(Abc4[[i]]$unadj.values)[ABC_97.5])
  ABC.CI_64_g[i,] <- c(sort(Abc64[[i]]$unadj.values)[ABC_2.5], sort(Abc64[[i]]$unadj.values)[ABC_97.5])

  }

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
ABC.2_gcov <- sum((trueGamma >= ABC.CI_2_g[,1]) == TRUE & (trueGamma <= ABC.CI_2_g[,2]) == TRUE)/modelRuns * 100
ABC.4_gcov <- sum((trueGamma >= ABC.CI_4_g[,1]) == TRUE & (trueGamma <= ABC.CI_4_g[,2]) == TRUE)/modelRuns * 100
ABC.64_gcov <- sum((trueGamma >= ABC.CI_64_g[,1]) == TRUE & (trueGamma <= ABC.CI_64_g[,2]) == TRUE)/modelRuns * 100

print("The coverage of the parameter estimate gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: gamma = ',  ABC.2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: gamma = ',  ABC.4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: gamma = ',  ABC.64_gcov, '%'))





      