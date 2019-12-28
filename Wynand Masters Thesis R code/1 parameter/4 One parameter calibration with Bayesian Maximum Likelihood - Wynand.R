## Wynand can Staden
## Data simulation with Bayesian Maximum Likelihood Estimation
## Copyright 2019

# Please note this file uses the sirModel() functions in the 
# 1. One parameter calibration with Least Squares.r file

library(gmp)        #for chooseZ() function
library(dplyr)      # for sample_n function
library(ggplot2)

baysianML <- function(randDraw, trueGamma, samSize = 1000){
  
  #initializing variable to store results
  BML2.loglik2 <- c()
  BML4.loglik4 <- c()
  BML64.loglik64 <- c()
  
  gammaPrior <- runif(randDraw, min = 0.01, max = 0.1)
  
  p2 <- matrix(c(0, 0), randDraw, 2)
  p4 <- matrix(replicate(4, 0), randDraw, 4)
  p64 <- matrix(replicate(64, 0), randDraw, 64)

  #model calibration
  for(i in 1:randDraw){
    
    p2[i,] <- sirModel2(gammaPrior[i])
    p4[i,] <- sirModel4(gammaPrior[i])
    p64[i,] <- sirModel64(gammaPrior[i])
    
    x2 <- sirModel2(trueGamma) * samSize
    x4 <- sirModel4(trueGamma) * samSize
    x64 <- sirModel64(trueGamma)* samSize
    
    #3.
    # L(p) = p^x*(1-p)^(n - x)
    # log(L) = xlog(p) + (n-x)log(1-p)
    
    #using the negative log-likelihood as likelihoods for the parameter combinations
    loglik2 <- c()
    loglik4 <- c()
    loglik64 <- c()

    for(j in 1:length(x2)){
      p2[i, j] <- ifelse(p2[i, j]==0, 0.0001, p2[i, j])


      loglik2[j] <- log(chooseZ(samSize, x2[j])) + (x2[j])*log(p2[i,j]) + (samSize-(x2[j]))*log(1-p2[i,j])

    }

    for(j in 1:length(x4)){
      p4[i, j] <- ifelse(p4[i, j]==0, 0.0001, p4[i, j])

      loglik4[j] <- log(chooseZ(samSize, x4[j])) + (x4[j])*log(p4[i,j]) + (samSize-(x4[j]))*log(1-p4[i,j])
    }

    for(j in 1:length(x64)){
      p64[i, j] <- ifelse(p64[i, j]==0, 0.0001, p64[i, j])

      loglik64[j] <- log(chooseZ(samSize, x64[j])) + (x64[j])*log(p64[i,j]) + (samSize-(x64[j]))*log(1-p64[i,j])
    }
    
    BML2.loglik2[i] <- sum(loglik2)  
    BML4.loglik4[i] <- sum(loglik4)
    BML64.loglik64[i] <- sum(loglik64)
    
    cat(paste0(i, ", "))
   
  }

  BMLE.result2 <- data.frame(gammaPrior, BML2.loglik2)
  BMLE.result4 <- data.frame(gammaPrior, BML4.loglik4)
  BMLE.result64 <- data.frame(gammaPrior, BML64.loglik64)
  
  
  return(list(BMLE.result2, BMLE.result4, BMLE.result64))

}


#####################################################################
## Model calibration
#####################################################################

# Now to run the calibration method testing
modelRuns <- 1000
randDraw <- 1000
trueGamma <- 0.02

BML.post <- list()
BMLE2 <- list()
BMLE4 <- list()
BMLE64 <- list()


for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  BML.post[[i]] <- baysianML(randDraw, trueGamma) 
  
  BMLE2[i] <- BML.post[[i]][1]
  BMLE4[i] <- BML.post[[i]][2]
  BMLE64[i] <- BML.post[[i]][3]
  
}


#4. Resampling from each individual posterior distribution using weights
BMLE2.weight2 <- list()
BMLE4.weight4 <- list()
BMLE64.weight64 <- list()

nameCols2 <- c("gammaPrior", "weight2")
nameCols4 <- c("gammaPrior", "weight4")
nameCols64 <- c("gammaPrior", "weight64")


## Weight calculation Method 
for(i in 1:modelRuns){
  
  weight2 <- exp(BMLE2[[i]]$BML2.loglik2)
  weight2 <- weight2/sum(weight2)
  BMLE2.weight2[[i]] <- data.frame(BMLE2[[i]]$gammaPrior, weight2)
  colnames(BMLE2.weight2[[i]]) <- nameCols2
  
  weight4 <- exp(BMLE4[[i]]$BML4.loglik4)
  weight4 <- weight4/sum(weight4)
  BMLE4.weight4[[i]] <- data.frame(BMLE4[[i]]$gammaPrior, weight4)
  colnames(BMLE4.weight4[[i]]) <- nameCols4
  
  weight64 <- exp(BMLE64[[i]]$BML64.loglik64)
  weight64 <- weight64/sum(weight64)
  BMLE64.weight64[[i]] <- data.frame(BMLE64[[i]]$gammaPrior, weight64)
  colnames(BMLE64.weight64[[i]]) <- nameCols64
  
}

#ReSample step and finding the median for each random Draw
BMLE.g.post.2 <- list()
BMLE.g.post.4 <- list()
BMLE.g.post.64 <- list()

resampleSize <- 1000

for(i in 1:modelRuns){
  
  BMLE.g.post.2[[i]] <- sample(BMLE2.weight2[[i]]$gammaPrior, size = resampleSize, replace = T, prob = BMLE2.weight2[[i]]$weight2) 
  BMLE.g.post.4[[i]] <- sample(BMLE4.weight4[[i]]$gammaPrior, size = resampleSize, replace = T, prob = BMLE4.weight4[[i]]$weight4) 
  BMLE.g.post.64[[i]] <- sample(BMLE64.weight64[[i]]$gammaPrior, size = resampleSize, replace = T, prob = BMLE64.weight64[[i]]$weight64) 
  
}

#####################################################################
## Measuring performance
#####################################################################



#calculating the median for the posterior distributions of gamma
BMLE.g.median.2 <- c()
BMLE.g.median.4 <- c()
BMLE.g.median.64 <- c()

for(i in 1:modelRuns){
  BMLE.g.median.2[i] <- median(BMLE.g.post.2[[i]])
  BMLE.g.median.4[i] <- median(BMLE.g.post.4[[i]])
  BMLE.g.median.64[i] <- median(BMLE.g.post.64[[i]])
}

# Calculating the average value of the median values 
BMLE.avgPar2 <- round(mean(BMLE.g.median.2), 3)
BMLE.avgPar4 <- round(mean(BMLE.g.median.4), 3)
BMLE.avgPar64 <- round(mean(BMLE.g.median.64), 3)

print("Average parameter estimate gamma:")
print(paste0('Model with 2 target features: gamma = ', BMLE.avgPar2))
print(paste0('Model with 4 target features: gamma = ', BMLE.avgPar4))
print(paste0('Model with 64 target features: gamma = ', BMLE.avgPar64))



#Calculating the Bias
BMLE.gBias2 <- BMLE.avgPar2 - trueGamma
BMLE.gBias4 <- BMLE.avgPar4 - trueGamma
BMLE.gBias64 <- BMLE.avgPar64 - trueGamma


print("The Percentage Bias of the parameter estimate gamma:")
print(paste0('Bias for Model with 2 target features:gamma Bias = ', BMLE.gBias2/trueGamma *100, '%' ))
print(paste0('Bias for Model with 4 target features:gamma Bias = ', BMLE.gBias4/trueGamma *100, '%' ))
print(paste0('Bias for Model with 64 target features:gamma Bias = ', BMLE.gBias64/trueGamma *100, '%' ))



#Calculating the accuracy using the Root Mean Square Error
BMLE.gAccu2 <- sqrt((sum((BMLE.g.median.2 - trueGamma)^2)/modelRuns))
BMLE.gAccu4 <- sqrt((sum((BMLE.g.median.4 - trueGamma)^2)/modelRuns))
BMLE.gAccu64 <- sqrt((sum((BMLE.g.median.64 - trueGamma)^2)/modelRuns))

print("The accuracy of the parameter estimate gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: gamma = ',  round(BMLE.gAccu2, 3)))
print(paste0('RMSE for Model with 4 target features: gamma = ',  round(BMLE.gAccu4, 3)))
print(paste0('RMSE for Model with 64 target features: gamma = ',  round(BMLE.gAccu64, 3)))


# Calculating credible intervals of the re-sampled posterior distributions
BML_2.5 <- resampleSize * 0.025
BML_97.5 <- resampleSize * 0.975

BML.CI_2_g <- matrix(c(0, 0), modelRuns, 2)
BML.CI_4_g <- matrix(c(0, 0), modelRuns, 2)
BML.CI_64_g <- matrix(c(0, 0), modelRuns, 2)


for(i in 1:resampleSize){
  BML.CI_2_g[i,] <- c(sort(BMLE.g.post.2[[i]])[BML_2.5], sort(BMLE.g.post.2[[i]])[BML_97.5])
  BML.CI_4_g[i,] <- c(sort(BMLE.g.post.4[[i]])[BML_2.5], sort(BMLE.g.post.4[[i]])[BML_97.5])
  BML.CI_64_g[i,] <- c(sort(BMLE.g.post.64[[i]])[BML_2.5], sort(BMLE.g.post.64[[i]])[BML_97.5])
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
BML.2_gcov <- sum((trueGamma >= BML.CI_2_g[,1]) == TRUE & (trueGamma <= BML.CI_2_g[,2]) == TRUE)/modelRuns * 100
BML.4_gcov <- sum((trueGamma >= BML.CI_4_g[,1]) == TRUE & (trueGamma <= BML.CI_4_g[,2]) == TRUE)/modelRuns * 100
BML.64_gcov <- sum((trueGamma >= BML.CI_64_g[,1]) == TRUE & (trueGamma <= BML.CI_64_g[,2]) == TRUE)/modelRuns * 100


print("The coverage of the parameter estimate gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: gamma = ',  BML.2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: gamma = ',  BML.4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: gamma = ',  BML.64_gcov, '%'))



