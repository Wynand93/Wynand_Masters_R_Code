## Wynand can Staden
## Data simulation with Bayesian Maximum Likelihood Estimation
## Copyright 2019

## Please note this file makes use of the sirModel functions in the 1. Data simulate for Least Squares.r file


library(gmp)        #for chooseZ() function
library(dplyr)      # for sample_n function

#This function generates a prior distribution for the parameters and assigns likelihood values to the parameter combinations
baysianML <- function(randDraw, betaGamma, samSize = 1000){
  
  #Initialzing variables to store results
  BML2.loglik2 <- c()
  BML4.loglik4 <- c()
  BML64.loglik64 <- c()
  
  betaPrior <- runif(randDraw, min = 0.01, max = 0.5)
  gammaPrior <- runif(randDraw, min = 0.01, max = 0.1)
  
  p2 <- matrix(c(0, 0), randDraw, 2)
  p4 <- matrix(replicate(4, 0), randDraw, 4)
  p64 <- matrix(replicate(64, 0), randDraw, 64)

  
  for(i in 1:randDraw){
  
      
    p2[i,] <- sirModel2(betaPrior[i], gammaPrior[i])
    p4[i,] <- sirModel4(betaPrior[i], gammaPrior[i])
    p64[i,] <- sirModel64(betaPrior[i], gammaPrior[i])
    
    x2 <- sirModel2(betaGamma[1], betaGamma[2]) * samSize
    x4 <- sirModel4(betaGamma[1], betaGamma[2]) * samSize
    x64 <- sirModel64(betaGamma[1], betaGamma[2])* samSize
    
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
  
  #Now to assosciate each parameter combination its log-likelihood
  BMLE.result2 <- data.frame(betaPrior,gammaPrior, BML2.loglik2)
  BMLE.result4 <- data.frame(betaPrior,gammaPrior, BML4.loglik4)
  BMLE.result64 <- data.frame(betaPrior, gammaPrior, BML64.loglik64)
  
  
  
  return(list(BMLE.result2, BMLE.result4, BMLE.result64))

}


#####################################################################
## Model Calibration
#####################################################################

modelRuns <- 1000
randDraw <- 1000
betaGamma <- c(0.2, 0.02)


BML.post <- list()
BMLE2 <- list()
BMLE4 <- list()
BMLE64 <- list()


for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  BML.post[[i]] <- baysianML(randDraw, betaGamma) 
  
  BMLE2[i] <- BML.post[[i]][1]
  BMLE4[i] <- BML.post[[i]][2]
  BMLE64[i] <- BML.post[[i]][3]
  
}


#4. Resampling from each individual posterior distribution using weights
BMLE2.weight2 <- list()
BMLE4.weight4 <- list()
BMLE64.weight64 <- list()

nameCols2 <- c("betaPrior", "gammaPrior", "weight2")
nameCols4 <- c("betaPrior", "gammaPrior", "weight4")
nameCols64 <- c("betaPrior", "gammaPrior", "weight64")


## Weight calculation Method 
for(i in 1:modelRuns){

  weight2 <- exp(BMLE2[[i]]$BML2.loglik2)
  weight2 <- weight2/sum(weight2)
  BMLE2.weight2[[i]] <- data.frame(BMLE2[[i]]$betaPrior, BMLE2[[i]]$gammaPrior, weight2)
  colnames(BMLE2.weight2[[i]]) <- nameCols2

  weight4 <- exp(BMLE4[[i]]$BML4.loglik4)
  weight4 <- weight4/sum(weight4)
  BMLE4.weight4[[i]] <- data.frame(BMLE4[[i]]$betaPrior, BMLE4[[i]]$gammaPrior, weight4)
  colnames(BMLE4.weight4[[i]]) <- nameCols4

   weight64 <- exp(BMLE64[[i]]$BML64.loglik64)
   weight64 <- weight64/sum(weight64)
   BMLE64.weight64[[i]] <- data.frame(BMLE64[[i]]$betaPrior, BMLE64[[i]]$gammaPrior, weight64)
   colnames(BMLE64.weight64[[i]]) <- nameCols64
   
}
## Please note: with samsize = 10000, the weight64 calculations produces Inf values

#ReSample step and finding the median for each random Draw
BMLE.post.2 <- list()
BMLE.post.4 <- list()
BMLE.post.64 <- list()

resampleSize <- 1000

for(i in 1:modelRuns){

  BMLE.post.2[[i]] <- sample_n(BMLE2.weight2[[i]], size = resampleSize, replace = T, weight = BMLE2.weight2[[i]]$weight2) 
  BMLE.post.4[[i]] <- sample_n(BMLE4.weight4[[i]], size = resampleSize, replace = T, weight = BMLE4.weight4[[i]]$weight4)
  BMLE.post.64[[i]] <- sample_n(BMLE64.weight64[[i]], size = resampleSize, replace = T, weight = BMLE64.weight64[[i]]$weight64) 
   
}


# The following code will calculate the median of the posterior distributions 
BMLE.median.2 <- matrix(c(0, 0), modelRuns, 2)
BMLE.median.4 <- matrix(c(0, 0), modelRuns, 2)
# BMLE.median.64 <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  BMLE.median.2[i,] <- c(median(BMLE.post.2[[i]]$betaPrior), median(BMLE.post.2[[i]]$gammaPrior))
  BMLE.median.4[i,] <- c(median(BMLE.post.4[[i]]$betaPrior), median(BMLE.post.4[[i]]$gammaPrior))
  BMLE.median.64[i,] <- c(median(BMLE.post.64[[i]]$betaPrior), median(BMLE.post.64[[i]]$gammaPrior))
   
}


#####################################################################
## Measuring performance
#####################################################################


# Calculating the average of the estimated median parameters
BMLE.avgPar2 <- c(round(mean(BMLE.median.2[,1]), 3) , round(mean(BMLE.median.2[,2]), 3))
BMLE.avgPar4 <- c(round(mean(BMLE.median.4[,1]), 3) , round(mean(BMLE.median.4[,2]), 3))
BMLE.avgPar64 <- c(round(mean(BMLE.median.64[,1]), 3) , round(mean(BMLE.median.64[,2]), 3))

print("Mean parameter estimates of beta and gamma:")
print(paste0('Model with 2 target features: beta = ', BMLE.avgPar2[1] , ', gamma = ', BMLE.avgPar2[2]))
print(paste0('Model with 4 target features: beta = ', BMLE.avgPar4[1] , ', gamma = ', BMLE.avgPar4[2]))
print(paste0('Model with 64 target features: beta = ', BMLE.avgPar64[1] , ', gamma = ', BMLE.avgPar64[2]))



#Calculating the Bias
BMLE.bgBias2 <- BMLE.avgPar2 - betaGamma
BMLE.bgBias4 <- BMLE.avgPar4 - betaGamma
BMLE.bgBias64 <- BMLE.avgPar64 - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 target features: beta Bias = ', BMLE.bgBias2[1]/betaGamma[1] *100, '%, gamma Bias = ',
             BMLE.bgBias2[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 4 target features: beta Bias = ', BMLE.bgBias4[1]/betaGamma[1] *100, '%, gamma Bias = ',
             BMLE.bgBias4[2]/betaGamma[2] *100, '%' ))
print(paste0('Bias for Model with 64 target features: beta Bias = ', BMLE.bgBias64[1]/betaGamma[1] *100, '%, gamma Bias = ',
             BMLE.bgBias64[2]/betaGamma[2] *100, '%' ))


#Calculating the accuracy using the Root Mean Square Error
BMLE.bgAccu2 <- c(sqrt((sum((BMLE.median.2[,1] - betaGamma[1])^2))/modelRuns), 
                 sqrt((sum((BMLE.median.2[,2] - betaGamma[2])^2))/modelRuns))

BMLE.bgAccu4 <- c(sqrt((sum((BMLE.median.4[,1] - betaGamma[1])^2))/modelRuns), 
                 sqrt((sum((BMLE.median.4[,2] - betaGamma[2])^2))/modelRuns))

BMLE.bgAccu64 <- c(sqrt((sum((BMLE.median.64[,1] - betaGamma[1])^2))/modelRuns), 
                  sqrt((sum((BMLE.median.64[,2] - betaGamma[2])^2))/modelRuns))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 target features: beta = ', round(BMLE.bgAccu2[1], 3), ' gamma = ',  round(BMLE.bgAccu2[2], 3)))
print(paste0('RMSE for Model with 4 target features: beta = ', round(BMLE.bgAccu4[1], 3), ' gamma = ',  round(BMLE.bgAccu4[2], 3)))
print(paste0('RMSE for Model with 64 target features: beta = ', round(BMLE.bgAccu64[1], 3), ' gamma = ',  round(BMLE.bgAccu64[2], 3)))


# Calculating credible intervals (within the re-sampled posterior distributions)
BML_2.5 <- resampleSize * 0.025
BML_97.5 <- resampleSize * 0.975

BML.CI_2_b <- matrix(c(0, 0), modelRuns, 2)
BML.CI_2_g <- matrix(c(0, 0), modelRuns, 2)

BML.CI_4_b <- matrix(c(0, 0), modelRuns, 2)
BML.CI_4_g <- matrix(c(0, 0), modelRuns, 2)

BML.CI_64_b <- matrix(c(0, 0), modelRuns, 2)
BML.CI_64_g <- matrix(c(0, 0), modelRuns, 2)

length(BMLE.post.2[[1]]$betaPrior)

for(i in 1:modelRuns){
  BML.CI_2_b[i,] <- c(sort(BMLE.post.2[[i]]$betaPrior)[BML_2.5], sort(BMLE.post.2[[i]]$betaPrior)[BML_97.5])
  BML.CI_2_g[i,] <- c(sort(BMLE.post.2[[i]]$gammaPrior)[BML_2.5], sort(BMLE.post.2[[i]]$gammaPrior)[BML_97.5])

  BML.CI_4_b[i,] <- c(sort(BMLE.post.4[[i]]$betaPrior)[BML_2.5], sort(BMLE.post.4[[i]]$betaPrior)[BML_97.5])
  BML.CI_4_g[i,] <- c(sort(BMLE.post.4[[i]]$gammaPrior)[BML_2.5], sort(BMLE.post.4[[i]]$gammaPrior)[BML_97.5])

  BML.CI_64_b[i,] <- c(sort(BMLE.post.64[[i]]$betaPrior)[BML_2.5], sort(BMLE.post.64[[i]]$betaPrior)[BML_97.5])
  BML.CI_64_g[i,] <- c(sort(BMLE.post.64[[i]]$gammaPrior)[BML_2.5], sort(BMLE.post.64[[i]]$gammaPrior)[BML_97.5])

}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
BML.2_bcov <- sum((betaGamma[1] >= BML.CI_2_b[,1]) == TRUE & (betaGamma[1] <= BML.CI_2_b[,2]) == TRUE)/modelRuns * 100
BML.2_gcov <- sum((betaGamma[2] >= BML.CI_2_g[,1]) == TRUE & (betaGamma[2] <= BML.CI_2_g[,2]) == TRUE)/modelRuns * 100

BML.4_bcov <- sum((betaGamma[1] >= BML.CI_4_b[,1]) == TRUE & (betaGamma[1] <= BML.CI_4_b[,2]) == TRUE)/modelRuns * 100
BML.4_gcov <- sum((betaGamma[2] >= BML.CI_4_g[,1]) == TRUE & (betaGamma[2] <= BML.CI_4_g[,2]) == TRUE)/modelRuns * 100

BML.64_bcov <- sum((betaGamma[1] >= BML.CI_64_b[,1]) == TRUE & (betaGamma[1] <= BML.CI_64_b[,2]) == TRUE)/modelRuns * 100
BML.64_gcov <- sum((betaGamma[2] >= BML.CI_64_g[,1]) == TRUE & (betaGamma[2] <= BML.CI_64_g[,2]) == TRUE)/modelRuns * 100

print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: beta = ', BML.2_bcov, '%, gamma = ',  BML.2_gcov, '%'))
print(paste0('Coverage for Model with 4 target features: beta = ', BML.4_bcov, '%, gamma = ',  BML.4_gcov, '%'))
print(paste0('Coverage for Model with 64 target features: beta = ', BML.64_bcov, '%, gamma = ',  BML.64_gcov, '%'))



