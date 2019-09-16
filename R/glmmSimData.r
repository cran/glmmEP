########## R function: glmmSimData ##########

# Simulates data from a generalized linear mixed model.

# Last changed: 22 JAN 2018

glmmSimData <- function(seed=12345)
{
   # Set seed:

   set.seed(seed)

   # Set sample sizes:

   m <- 250
   nVec <- sample(20:30,m,replace=TRUE)

   # Set true values of parameters:

   beta0True <- 0.37
   beta1True <- 0.93
   beta2True <- -0.46
   beta3True <- 0.08
   beta4True <- -1.34
   beta5True <- 1.09

   SigmaTrue <- rbind(c(0.53,-0.36),
                      c(-0.36,0.92))
  
   sigma1True <- sqrt(SigmaTrue[1,1])
   sigma2True <- sqrt(SigmaTrue[2,2])
   rhoTrue <- SigmaTrue[1,2]/(sigma1True*sigma2True)  

   SVDobj <- svd(SigmaTrue)
   sqrtSigmaTrue <- t(SVDobj$v%*%(t(SVDobj$u)*sqrt(SVDobj$d)))
   uMat <- matrix(rnorm(2*m),m,2)%*%sqrtSigmaTrue

   x1 <- rep(NA,sum(nVec))
   x2 <- rep(NA,sum(nVec))
   x3 <- rep(NA,sum(nVec))
   x4 <- rep(NA,sum(nVec))
   x5 <- rep(NA,sum(nVec))

   y <- rep(NA,sum(nVec))
   idNum <- rep(NA,sum(nVec))

   sttPos <- 1
   for (i in 1:m)
   {
      endPos <- sttPos + nVec[i] -1

      x1Curr <- runif(nVec[i]) ; x1[sttPos:endPos] <- x1Curr
      x2Curr <- runif(nVec[i]) ; x2[sttPos:endPos] <- x2Curr
      x3Curr <- runif(nVec[i]) ; x3[sttPos:endPos] <- x3Curr
      x4Curr <- runif(nVec[i]) ; x4[sttPos:endPos] <- x4Curr
      x5Curr <- runif(nVec[i]) ; x5[sttPos:endPos] <- x5Curr

      linPred <- (beta0True + uMat[i,1] + (beta1True + uMat[i,2])*x1Curr
                  + beta2True*x2Curr + beta3True*x3Curr + beta4True*x4Curr
                  + beta5True*x5Curr)

      y[sttPos:endPos] <- rbinom(nVec[i],1,pnorm(linPred))
      idNum[sttPos:endPos] <- i
      sttPos <- endPos + 1
   }

   Xfixed <- cbind(1,x1,x2,x3,x4,x5)
   Xrandom <- cbind(1,x1)

   # Return data vectors:

   return(list(y=y,Xfixed=Xfixed,Xrandom=Xrandom,idNum=idNum))
}

############ End of glmmSimData ############

