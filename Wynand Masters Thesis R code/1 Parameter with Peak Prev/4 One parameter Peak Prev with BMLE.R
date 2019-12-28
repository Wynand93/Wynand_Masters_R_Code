## Wynand can Staden
## Bayesian Maximum Likelihood Estimation calibration method with Peak prevalence
## Copyright 2019

## Please note: This file makes use of the sirModelPeakPrev function in the One parameter 1 Peak Prev with LS.r file

library(gmp)

#creating a function to call and use the BMLE method
baysianMLPeak <- function(randDraw, trueGamma, samSize = 1000){
  
  BMLPeak <- c()
  
  gammaPrior <- runif(randDraw, min = 0.01, max = 0.1)
  
  pPeak <- matrix(c(0, 0, 0), randDraw, 3)

  # L(p) = p^x*(1-p)^(n - x)
  # log(L) = xlog(p) + (n-x)log(1-p)
  
  for(i in 1:randDraw){
    
    pPeak[i,] <- sirModelPeakPrev(gammaPrior[i])
    xPeak <- sirModelPeakPrev(trueGamma) * samSize
    loglikPeak <- c()
    
    for(j in 1:length(xPeak)){
      
      pPeak[i, j] <- ifelse(pPeak[i, j]==0, 0.0001, pPeak[i, j])
      loglikPeak[j] <- log(chooseZ(samSize, xPeak[j])) + (xPeak[j])*log(pPeak[i,j]) + (samSize-(xPeak[j]))*log(1-pPeak[i,j])
      
    }
    
    BMLPeak[i] <- sum(loglikPeak)  
    
    cat(paste0(i, ", ")) #to count the numer of runs
    
  }
  
  BMLE.resultPeak <- data.frame(gammaPrior, BMLPeak)
  
  return(BMLE.resultPeak)
  
}

#####################################################################
## Model Calibration
#####################################################################

modelRuns <- 1000
randDraw <- 1000
trueGamma <- 0.02

BMLE.g.Peak <- list()

for(i in 1:modelRuns){
  print(paste0("Model run: ", i  ))
  
  BMLE.g.Peak[[i]] <- baysianMLPeak(randDraw, trueGamma) 
  
}


#4.
BMLE.g.Peak.weight <- list()

nameCols.g.Peak <- c("gammaPrior", "weightPeak")

## Weight calculation Method 
for(i in 1:modelRuns){
  
  weight.g.Peak <- exp(BMLE.g.Peak[[i]]$BMLPeak)
  BMLE.g.Peak.weight[[i]] <- data.frame(BMLE.g.Peak[[i]]$gammaPrior, weight.g.Peak)
  colnames(BMLE.g.Peak.weight[[i]]) <- nameCols.g.Peak
  
}

#4.
#ReSampling from prior distribution using weights for parameter combinations
BMLE.post.g.Peak <- list()

resampleSize <- 1000

for(i in 1:modelRuns){
  
  BMLE.post.g.Peak[[i]] <- sample(BMLE.g.Peak.weight[[i]]$gammaPrior, size = resampleSize, replace = T, prob = BMLE.g.Peak.weight[[i]]$weightPeak) 
  
}


#Now to calculate the mdeian of the resampled posterior distributions
BMLE.g.medianPeak <- c()

for(i in 1:modelRuns){
  
  BMLE.g.medianPeak[i] <- median(BMLE.post.g.Peak[[i]])
}

#####################################################################
## Measuring performance
#####################################################################


# Calculating the average calibrated parameteres
BMLE.avgParPeak <- round(mean(BMLE.g.medianPeak), 3)

print("Average parameter estimategamma:")
print(paste0('Model with 2 + Peak target features: gamma = ', BMLE.avgParPeak))


#Calculating the Bias
BMLE.gBiasPeak <- BMLE.avgParPeak - trueGamma


print("The Percentage Bias of the parameter estimate of gamma:")
print(paste0('Bias for Model with 2 + Peak target features: gamma Bias = ', BMLE.gBiasPeak/trueGamma *100, '%' ))


#Calculating the accuracy using the Root Mean Square Error
BMLE.gAccuPeak <- sqrt((sum((BMLE.g.medianPeak - trueGamma)^2)/modelRuns))

print("The accuracy of the parameter estimates of gamma using RMSE:")
print(paste0('RMSE for Model with 2 + Peak target features: gamma = ',  round(BMLE.gAccuPeak, 3)))


# Calculating credible intervals
BML_2.5 <- resampleSize * 0.025
BML_97.5 <- resampleSize * 0.975

BML.CI_Peak_g <- matrix(c(0, 0), modelRuns, 2)

for(i in 1:modelRuns){
  
  BML.CI_Peak_g[i,] <- c(sort(BMLE.post.g.Peak[[i]])[BML_2.5], sort(BMLE.post.g.Peak[[i]])[BML_97.5])
  
}

# Now to calculate coverage of the true estimate given the confidence intervals of the parameter estimates
BML.Peak_gcov <- sum((trueGamma >= BML.CI_Peak_g[,1]) == TRUE & (trueGamma <= BML.CI_Peak_g[,2]) == TRUE)/modelRuns * 100


print("The coverage of the parameter estimates of gamma given the CI's:")
print(paste0('Coverage for Model with 2 target features: gamma = ',  BML.Peak_gcov, '%'))



