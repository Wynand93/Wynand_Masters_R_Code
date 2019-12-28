## Wynand can Staden
## Bayesian Maximum Likelihood Estimation calibration method with Peak prevalence
## Copyright 2019

#Please note this file uses the sirModelPeakPrev function in the 1 Peak Prev with LS.r file

library(gmp)        #for chooseZ() function
library(dplyr)      # for sample_n function

#Creating a function to call and use the BMLE method
baysianMLPeak <- function(randDraw, betaGamma, samSize = 1000){
  
  BMLPeak <- c()
  
  betaPrior <- runif(randDraw, min = 0.01, max = 0.5)
  gammaPrior <- runif(randDraw, min = 0.01, max = 0.1)
  
  pPeak <- matrix(c(0, 0, 0), randDraw, 3)
  
  for(i in 1:randDraw){
    
    pPeak[i,] <- sirModelPeakPrev(betaPrior[i], gammaPrior[i])
    
    xPeak <- sirModelPeakPrev(betaGamma[1], betaGamma[2]) * samSize
    
    loglikPeak <- c()
    
    for(j in 1:length(xPeak)){
      
      pPeak[i, j] <- ifelse(pPeak[i, j]==0, 0.0001, pPeak[i, j])
      loglikPeak[j] <- log(chooseZ(samSize, xPeak[j])) + (xPeak[j])*log(pPeak[i,j]) + (samSize-(xPeak[j]))*log(1-pPeak[i,j])
      
    }
    
    BMLPeak[i] <- sum(loglikPeak)  
    
    cat(paste0(i, ", "))
    
  }
  
  #now to assign negative log-likelihood values to each of the parameter combinations
  BMLE.resultPeak <- data.frame(betaPrior, gammaPrior, BMLPeak)
  
  return(BMLE.resultPeak)
  
}


#####################################################################
## Model Calibration
#####################################################################

modelRuns <- 1000
randDraw <- 1000
betaGamma <- c(0.2, 0.02)

BMLEPeak <- list()

for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  BMLEPeak[[i]] <- baysianMLPeak(randDraw, betaGamma) 
  
}


BMLEPeak.weight <- list()

nameColsPeak <- c("betaPrior", "gammaPrior", "weightPeak")

## Weight calculation Method 
for(i in 1:modelRuns){
  
  weightPeak <- exp(BMLEPeak[[i]]$BMLPeak)
  BMLEPeak.weight[[i]] <- data.frame(BMLEPeak[[i]]$betaPrior, BMLEPeak[[i]]$gammaPrior, weightPeak)
  colnames(BMLEPeak.weight[[i]]) <- nameColsPeak
 
}

#ReSampling from prior distribution using weights for parameter combinations
BMLE.postPeak <- list()

resampleSize <- 1000

for(i in 1:modelRuns){
  
  BMLE.postPeak[[i]] <- sample_n(BMLEPeak.weight[[i]], size = resampleSize, replace = T, weight = BMLEPeak.weight[[i]]$weightPeak) 
  
}

#Now calculating the median of the weighted posterior distributions
BMLE.median.Peak <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  BMLE.median.Peak[i,] <- c(median(BMLE.postPeak[[i]]$betaPrior), median(BMLE.postPeak[[i]]$gammaPrior))  
  
}


#####################################################################
## Measuring performance
#####################################################################

# Calculating the average of the median of the estimated posterior distributions
BMLE.avgParPeak <- c(round(mean(BMLE.median.Peak[,1]), 3) , round(mean(BMLE.median.Peak[,2]), 3))

print("Average parameter estimates of beta and gamma:")
print(paste0('Model with 2 + Peak target features: beta = ', BMLE.avgParPeak[1] , ', gamma = ', BMLE.avgParPeak[2]))


#Calculating the Bias
BMLE.bgBiasPeak <- BMLE.avgParPeak - betaGamma


print("The Percentage Bias of each parameter estimates of beta and gamma:")
print(paste0('Bias for Model with 2 + Peak target features: beta Bias = ', BMLE.bgBiasPeak[1]/betaGamma[1] *100, '%, gamma Bias = ',
             BMLE.bgBiasPeak[2]/betaGamma[2] *100, '%' ))


#Calculating the accuracy using the Root Mean Square Error
BMLE.bgAccuPeak <- c(sqrt((sum((BMLE.median.Peak[,1] - betaGamma[1])^2)/modelRuns)), 
                  sqrt((sum((BMLE.median.Peak[,2] - betaGamma[2])^2)/modelRuns)))

print("The accuracy of each parameter estimates of beta and gamma using RMSE:")
print(paste0('RMSE for Model with 2 + Peak target features: beta = ', round(BMLE.bgAccuPeak[1], 3), ' gamma = ',  round(BMLE.bgAccuPeak[2], 3)))


# Calculating credible intervals
BML_2.5 <- resampleSize * 0.025
BML_97.5 <- resampleSize * 0.975

BML.CI_Peak_b <- matrix(c(0, 0), modelRuns, 2)
BML.CI_Peak_g <- matrix(c(0, 0), modelRuns, 2)


for(i in 1:modelRuns){
  BML.CI_Peak_b[i,] <- c(sort(BMLE.postPeak[[i]]$betaPrior)[BML_2.5], sort(BMLE.postPeak[[i]]$betaPrior)[BML_97.5])
  BML.CI_Peak_g[i,] <- c(sort(BMLE.postPeak[[i]]$gammaPrior)[BML_2.5], sort(BMLE.postPeak[[i]]$gammaPrior)[BML_97.5])
  
  }

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
BML.Peak_bcov <- sum((betaGamma[1] >= BML.CI_Peak_b[,1]) == TRUE & (betaGamma[1] <= BML.CI_Peak_b[,2]) == TRUE)/modelRuns * 100
BML.Peak_gcov <- sum((betaGamma[2] >= BML.CI_Peak_g[,1]) == TRUE & (betaGamma[2] <= BML.CI_Peak_g[,2]) == TRUE)/modelRuns * 100


print("The coverage of each parameter estimates of beta and gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: beta = ', BML.Peak_bcov, '%, gamma = ',  BML.Peak_gcov, '%'))



