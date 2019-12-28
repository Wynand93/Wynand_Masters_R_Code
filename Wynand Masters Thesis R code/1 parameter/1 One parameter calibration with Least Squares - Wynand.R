## Wynand van Staden
## Model - Data simulation & Least Squares estimation to find one parameter
## Copyright 2019

#sampleSize of the data in the study
samSize <- 1000

#SIR Model with 2 target features
sirModel2 <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)

  model <- SIR(u0, tspan = seq(0,100,by=1), beta= beta, gamma=gamma)
  result <- run(model)
 
  individuals50 <- c(rep("S", result@U[1, 50]), rep("I", result@U[2, 50]), rep("R", result@U[3, 50]))
  individuals65 <- c(rep("S", result@U[1, 65]), rep("I", result@U[2, 65]), rep("R", result@U[3, 65]))
  
  
  samplePop <- c(summary(as.factor(sample(individuals50, size = sampleSize))), summary(as.factor(sample(individuals65, size = sampleSize)))) 
  
  pop <- samplePop[names(samplePop) == "I"] #let it return zero instead of numeric(0)
  
  for(i in 1:2){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
}


#SIR Model with 4 target features
sirModel4 <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= beta, gamma=gamma)
  result <- run(model)
  
  individuals30 <- c(rep("S", result@U[1, 30]), rep("I", result@U[2, 30]), rep("R", result@U[3, 30]))
  individuals45 <- c(rep("S", result@U[1, 45]), rep("I", result@U[2, 45]), rep("R", result@U[3, 45]))
  individuals60 <- c(rep("S", result@U[1, 60]), rep("I", result@U[2, 60]), rep("R", result@U[3, 60]))
  individuals75 <- c(rep("S", result@U[1, 75]), rep("I", result@U[2, 75]), rep("R", result@U[3, 75]))
  
  
  samplePop <- c(summary(as.factor(sample(individuals30, size = sampleSize))), summary(as.factor(sample(individuals45, size = sampleSize))), summary(as.factor(sample(individuals60, size = sampleSize))), summary(as.factor(sample(individuals75, size = sampleSize)))) 
  
  pop <- samplePop[names(samplePop) == "I"]
  
  for(i in 1:4){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
  
}

#SIR Model with 64 target features
sirModel64 <- function(gamma, beta = 0.2, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)
  
  model <- SIR(u0, tspan = seq(0,100,by=1), beta= beta, gamma=gamma)
  result <- run(model)
  
  individuals1.64 <- list()
  
  for(i in 1:64){
    
    individuals1.64[[i]] <- c(rep("S", result@U[1, i]), rep("I", result@U[2, i]), rep("R", result@U[3, i]))
  
  }

  samplePop <- c()
  for(i in 1:64){
  
    samplePop <- c(samplePop, c(summary(as.factor(sample(individuals1.64[[i]], size = sampleSize)))))
    
  }

  pop <- samplePop[names(samplePop) == "I"]
  
  for(i in 1:64){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
}


###############################################################################################
## Model Calibration
###############################################################################################
trueGamma <- 0.02
calibModelRuns <- 1000

#Initializing variables to store results
infPop2 <- matrix(c(0, 0), calibModelRuns, 2)
infPop4 <- matrix(replicate(4, 0), calibModelRuns, 4)
infPop64 <- matrix(replicate(64, 0), calibModelRuns, 64)

result2 <- list()
result4 <- list()
result64 <- list()


#to just store beta and gamma parameters from the reuslt
LS.result2 <- matrix(c(0), calibModelRuns, 1)
LS.result4 <- matrix(c(0), calibModelRuns, 1)
LS.result64 <- matrix(c(0), calibModelRuns, 1)

print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {
  
  infPop2[i,] <- sirModel2(trueGamma)
  infPop4[i,] <- sirModel4(trueGamma)
  infPop64[i,] <- sirModel64(trueGamma)
  
  f.opt2 <- function(params) sum((sirModel2(params)- infPop2[i,])^2 )
  f.opt4 <- function(params) sum((sirModel4(params)- infPop4[i,])^2 )
  f.opt64 <- function(params) sum((sirModel64(params)- infPop64[i,])^2 )
  
  result2[[i]] <- optim(runif(1, 0.01, 0.1), f.opt2,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000), hessian = T)
  LS.result2[i, ] <- result2[[i]]$par 

  result4[[i]] <- optim(runif(1, 0.01, 0.1), f.opt4,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000),hessian = T)
  LS.result4[i, ] <- result4[[i]]$par

  result64[[i]] <- optim(runif(1, 0.01, 0.1), f.opt64,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000), hessian = T)
  LS.result64[i, ] <- result64[[i]]$par


  print(i) #to check number of runs
  
}

##############################################################################
## Measuring performance
##############################################################################

# Calculating the average calibrated parameteres
avgPar2 <- round(mean(LS.result2), 3) 
avgPar4 <- round(mean(LS.result4), 3)
avgPar64 <- round(mean(LS.result64), 3)

print("Average parameter estimates of gamma:")
print(paste0('Model with 2 target features: gamma = ', avgPar2))
print(paste0('Model with 4 target features: gamma = ', avgPar4))
print(paste0('Model with 64 target features: gamma = ', avgPar64))

#Calculating the Bias
gBias2 <- avgPar2 - trueGamma
gBias4 <- avgPar4 - trueGamma
gBias64 <- avgPar64 - trueGamma

print("The Percentage Bias of the parameter estimate gamma:")
print(paste0('Bias for Model with 2 target features: ','gamma Bias = ', gBias2/trueGamma *100, '%' ))
print(paste0('Bias for Model with 4 target features: ','gamma Bias = ', gBias4/trueGamma *100, '%' ))
print(paste0('Bias for Model with 64 target features: ','gamma Bias = ', gBias64/trueGamma *100, '%' ))



#Calculating the accuracy using the Root Mean Square Error
gAccu2 <- sqrt((sum((LS.result2 - trueGamma)^2)/calibModelRuns)) 
gAccu4 <- sqrt((sum((LS.result4 - trueGamma)^2)/calibModelRuns)) 
gAccu64 <- sqrt((sum((LS.result64 - trueGamma)^2)/calibModelRuns)) 


print("The accuracy of the parameter estimate gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: ', 'gamma = ',  round(gAccu2, 3)))
print(paste0('RMSE for Model with 4 target features: ', 'gamma = ',  round(gAccu4, 3)))
print(paste0('RMSE for Model with 64 target features: ', 'gamma = ',  round(gAccu64, 3)))



#Calculating the coverage using confidence intervals

standError2 <- matrix(c(0), calibModelRuns, 1)
standError4 <- matrix(c(0), calibModelRuns, 1)
standError64 <- matrix(c(0), calibModelRuns, 1)


for(i in 1:calibModelRuns){
  standError2[i,] <- sqrt(abs(diag(solve(result2[[3]]$hessian))))
  standError4[i,] <- sqrt(abs(diag(solve(result4[[i]]$hessian))))
  standError64[i,] <- sqrt(abs(diag(solve(result64[[i]]$hessian))))
}


# Now confidence intervals of each of the parameter estimates
CI_2_g <- matrix(c(0, 0), calibModelRuns, 2)
CI_4_g <- matrix(c(0, 0), calibModelRuns, 2)
CI_64_g <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){
  
  CI_2_g[i,1] <- LS.result2[i] - 1.96*standError2[i]
  CI_2_g[i,2] <- LS.result2[i] + 1.96*standError2[i]

  CI_4_g[i,1] <- LS.result4[i] - 1.96*standError4[i]
  CI_4_g[i,2] <- LS.result4[i] + 1.96*standError4[i]
  
  CI_64_g[i,1] <- LS.result64[i] - 1.96*standError64[i]
  CI_64_g[i,2] <- LS.result64[i] + 1.96*standError64[i]
  
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
LS2_gcov <- sum((trueGamma >= CI_2_g[,1]) == TRUE & (trueGamma <= CI_2_g[,2]) == TRUE)/calibModelRuns * 100

LS4_gcov <- sum((trueGamma >= CI_4_g[,1]) == TRUE & (trueGamma <= CI_4_g[,2]) == TRUE)/calibModelRuns * 100

LS64_gcov <- sum((trueGamma >= CI_64_g[,1]) == TRUE & (trueGamma <= CI_64_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of the parameter estimate gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: gamma = ',  LS2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: gamma = ',  LS4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: gamma = ',  LS64_gcov, '%'))






