## Wynand van Staden
## Model - Data simulation & Least Squares estimation
## Copyright 2019

samSize <- 1000

#SIR Model with 2 target features
sirModel2 <- function(beta, gamma, N = 10000, inf = 0.1, sampleSize = samSize){
  library(SimInf)
  
  u0 <- data.frame(S=N*(1-inf), I=N*inf, R=0)

  model <- SIR(u0, tspan = seq(0,100,by=1), beta= beta, gamma=gamma)
  result <- run(model)
 
  individuals50 <- c(rep("S", result@U[1, 50]), rep("I", result@U[2, 50]), rep("R", result@U[3, 50]))
  individuals65 <- c(rep("S", result@U[1, 65]), rep("I", result@U[2, 65]), rep("R", result@U[3, 65]))
  
  
  samplePop <- c(summary(as.factor(sample(individuals50, size = sampleSize))), summary(as.factor(sample(individuals65, size = sampleSize)))) 
  
  pop <- samplePop[names(samplePop) == "I"] #lret it return zero instead of numeric(0)
  
  for(i in 1:2){
    if(is.na(pop[i])){
      pop[i] <- 0
    }
  }
  
  
  return(pop/sampleSize)
  
}

#SIR Model with 4 target features
sirModel4 <- function(beta, gamma, N = 10000, inf = 0.1, sampleSize = samSize){
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
sirModel64 <- function(beta, gamma, N = 10000, inf = 0.1, sampleSize = samSize){
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


#####################################################################
## Model Calibration
#####################################################################

betaGamma <- c(0.2, 0.02)
calibModelRuns <- 1000

#Initializing variables to store results
infPop2 <- matrix(c(0, 0), calibModelRuns, 2)
infPop4 <- matrix(replicate(4, 0), calibModelRuns, 4)
infPop64 <- matrix(replicate(64, 0), calibModelRuns, 64)

result2 <- list()
result4 <- list()
result64 <- list()


#to just store beta and gamma parameters from the reuslt
LS.result2 <- matrix(c(0, 0), calibModelRuns, 2)
LS.result4 <- matrix(c(0, 0), calibModelRuns, 2)
LS.result64 <- matrix(c(0, 0), calibModelRuns, 2)

print("Data Simulation Run Counter")
for (i in 1:calibModelRuns) {
  
  infPop2[i,] <- sirModel2(betaGamma[1], betaGamma[2])
  infPop4[i,] <- sirModel4(betaGamma[1], betaGamma[2])
  infPop64[i,] <- sirModel64(betaGamma[1], betaGamma[2])
  
  f.opt2 <- function(params) sum((sirModel2(params[1],params[2])- infPop2[i,])^2 )
  f.opt4 <- function(params) sum((sirModel4(params[1],params[2])- infPop4[i,])^2 )
  f.opt64 <- function(params) sum((sirModel64(params[1],params[2])- infPop64[i,])^2 )
  
  result2[[i]] <- optim(c(runif(1, 0.01, 0.5), runif(1, 0.01, 0.1)), f.opt2,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000), hessian = T)
  LS.result2[i, ] <- result2[[i]]$par 
  
  result4[[i]] <- optim(c(runif(1, 0.01, 0.5), runif(1, 0.01, 0.1)), f.opt4,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000),hessian = T)
  LS.result4[i, ] <- result4[[i]]$par
  
  result64[[i]] <- optim(c(runif(1, 0.01, 0.5), runif(1, 0.01, 0.1)), f.opt64,  method = "Nelder-Mead",lower = c(0.01,0.01), upper = c(0.5, 0.1), control=list(maxit=1000), hessian = T)
  LS.result64[i, ] <- result64[[i]]$par 
    
  
  print(i) #to check number of runs
  
}

#####################################################################
## Measuring performance
#####################################################################

# Calculating the average calibrated parameteres
avgPar2 <- c(round(mean(LS.result2[,1]), 3) , round(mean(LS.result2[, 2]), 3))
avgPar4 <- c(round(mean(LS.result4[,1]), 3) , round(mean(LS.result4[, 2]), 3))
avgPar64 <- c(round(mean(LS.result64[,1]), 3) , round(mean(LS.result64[, 2]), 3))

print("Average parameter estimates of beta and gamma:")
print(paste0('Model with 2 target features: beta = ', avgPar2[1] , ', gamma = ', avgPar2[2]))
print(paste0('Model with 4 target features: beta = ', avgPar4[1] , ', gamma = ', avgPar4[2]))
print(paste0('Model with 64 target features: beta = ', avgPar64[1] , ', gamma = ', avgPar64[2]))



#Calculating the Bias
bgBias2 <- avgPar2 - betaGamma
bgBias4 <- avgPar4 - betaGamma
bgBias64 <- avgPar64 - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 target features: beta Bias = ', bgBias2[1]/betaGamma[1] *100, '%, gamma Bias = ',
             bgBias2[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 4 target features: beta Bias = ', bgBias4[1]/betaGamma[1] *100, '%, gamma Bias = ',
             bgBias4[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 64 target features: beta Bias = ', bgBias64[1]/betaGamma[1] *100, '%, gamma Bias = ',
             bgBias64[2]/betaGamma[2] *100, '%' ))


#Calculating the accuracy using the Root Mean Square Error
bgAccu2 <- c(sqrt((sum((LS.result2[, 1] - betaGamma[1])^2)/calibModelRuns)), 
             sqrt((sum((LS.result2[, 2] - betaGamma[2])^2)/calibModelRuns)))

bgAccu4 <- c(sqrt((sum((LS.result4[, 1] - betaGamma[1])^2)/calibModelRuns)), 
             sqrt((sum((LS.result4[, 2] - betaGamma[2])^2)/calibModelRuns)))

bgAccu64 <- c(sqrt((sum((LS.result64[, 1] - betaGamma[1])^2)/calibModelRuns)), 
             sqrt((sum((LS.result64[, 2] - betaGamma[2])^2)/calibModelRuns)))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: beta = ', round(bgAccu2[1], 3), ' gamma = ',  round(bgAccu2[2], 3)))
print(paste0('RMSE for Model with 4 target features: beta = ', round(bgAccu4[1], 3), ' gamma = ',  round(bgAccu4[2], 3)))
print(paste0('RMSE for Model with 64 target features: beta = ', round(bgAccu64[1], 3), ' gamma = ',  round(bgAccu64[2], 3)))



#Calculating the coverage using confidence intervals
standError2 <- matrix(c(0, 0), calibModelRuns, 2)
standError4 <- matrix(c(0, 0), calibModelRuns, 2)
standError64 <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){
  standError2[i,] <- sqrt(abs(diag(solve(result2[[i]]$hessian))))
  standError4[i,] <- sqrt(abs(diag(solve(result4[[i]]$hessian))))
  standError64[i,] <- sqrt(abs(diag(solve(result64[[i]]$hessian))))
}


# Now confidence intervals of each of the parameter estimates
CI_2_b <- matrix(c(0, 0), calibModelRuns, 2)
CI_2_g <- matrix(c(0, 0), calibModelRuns, 2)

CI_4_b <- matrix(c(0, 0), calibModelRuns, 2)
CI_4_g <- matrix(c(0, 0), calibModelRuns, 2)

CI_64_b <- matrix(c(0, 0), calibModelRuns, 2)
CI_64_g <- matrix(c(0, 0), calibModelRuns, 2)


for(i in 1:calibModelRuns){
  
  CI_2_b[i,] <- c(LS.result2[i,1] - 1.96*standError2[i,1], LS.result2[i,1] + 1.96*standError2[i,1])
  CI_2_g[i,] <- c(LS.result2[i,2] - 1.96*standError2[i,2], LS.result2[i,2] + 1.96*standError2[i,2])
  
  CI_4_b[i,] <- c(LS.result4[i,1] - 1.96*standError4[i,1], LS.result4[i,1] + 1.96*standError4[i,1])
  CI_4_g[i,] <- c(LS.result4[i,2] - 1.96*standError4[i,2], LS.result4[i,2] + 1.96*standError4[i,2])
  
  CI_64_b[i,] <- c(LS.result64[i,1] - 1.96*standError64[i,1], LS.result64[i,1] + 1.96*standError64[i,1])
  CI_64_g[i,] <- c(LS.result64[i,2] - 1.96*standError64[i,2], LS.result64[i,2] + 1.96*standError64[i,2])
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
LS2_bcov <- sum((betaGamma[1] >= CI_2_b[,1]) == TRUE & (betaGamma[1] <= CI_2_b[,2]) == TRUE)/calibModelRuns * 100
LS2_gcov <- sum((betaGamma[2] >= CI_2_g[,1]) == TRUE & (betaGamma[2] <= CI_2_g[,2]) == TRUE)/calibModelRuns * 100

LS4_bcov <- sum((betaGamma[1] >= CI_4_b[,1]) == TRUE & (betaGamma[1] <= CI_4_b[,2]) == TRUE)/calibModelRuns * 100
LS4_gcov <- sum((betaGamma[2] >= CI_4_g[,1]) == TRUE & (betaGamma[2] <= CI_4_g[,2]) == TRUE)/calibModelRuns * 100

LS64_bcov <- sum((betaGamma[1] >= CI_64_b[,1]) == TRUE & (betaGamma[1] <= CI_64_b[,2]) == TRUE)/calibModelRuns * 100
LS64_gcov <- sum((betaGamma[2] >= CI_64_g[,1]) == TRUE & (betaGamma[2] <= CI_64_g[,2]) == TRUE)/calibModelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: beta = ', LS2_bcov, '%, gamma = ',  LS2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: beta = ', LS4_bcov, '%, gamma = ',  LS4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: beta = ', LS64_bcov, '%, gamma = ',  LS64_gcov, '%'))


