########## R function: glmmEP ##########

# For performing frequentist generalised linear mixed model
# analysis via expectation propagation.

# Last changed: 01 FEB 2018

glmmEP <- function(y,Xfixed=NULL,Xrandom,idNum,control=glmmEP.control())
{
   # Convert 'idNum' to sequential integers:

   idNumInp <- idNum
   idNum <- match(idNumInp,unique(idNumInp))
   
   # Set dimension variables:

   m <- length(unique(idNum))
   nVec <- rep(NA,m)
   for (i in 1:m)
      nVec[i] <- length(idNum[idNum==i])
   numObs <- length(y)

   # Compute vector of starting indices for each group:

   diffVec <- c(0,diff(idNum))
   indStt <- c(1,(1:length(diffVec))[diffVec==1])

   # Obtain the dimensions of the fixed and random effects components:

   if (is.null(Xfixed)) dF <- 0
   if (!is.null(Xfixed)) dF <- ncol(Xfixed)
   dR <- ncol(Xrandom)

   # Unpack the control values:

   confLev <- control$confLev
   BFGSmaxit <- control$BFGSmaxit
   BFGSreltol <- control$BFGSreltol
   EPmaxit <- control$EPmaxit
   EPreltol <- control$EPreltol
   NMmaxit <- control$NMmaxit
   NMreltol <- control$NMreltol
   quiet <- control$quiet
   preTransfData <- control$preTransfData

   # Do intercept checking:

   if (dF>0)
   {
      XfCol1 <- as.vector(Xfixed[,1])
      if (any(XfCol1!=1)) stop("The first column of Xfixed must have all entries equal to 1.")
   }
   XrCol1 <- as.vector(Xrandom[,1])
   if (any(XrCol1!=1)) stop("The first column of Xrandom must have all entries equal to 1.")

   # Make sure that first entry of colnames(Xfixed) is "intercept":

   if (colnames(Xfixed)[1]!="intercept") colnames(Xfixed)[1] <- "intercept"

   if (preTransfData)
   {
      # Scale all predictors to the unit interval for fitting:

      if (dF>1)
      {
         XfixedOrig <- as.matrix(Xfixed)

         minVecF <- apply(as.matrix(XfixedOrig[,-1]),2,min)
         maxVecF <- apply(as.matrix(XfixedOrig[,-1]),2,max)  
         denVecF <- maxVecF - minVecF

         if (any(denVecF==0)) stop("A non-intercept column of Xfixed has all entries identical.")

         minMat <- matrix(rep(minVecF,each=nrow(Xfixed)),nrow(Xfixed),dF-1)
         denMat <- matrix(rep(denVecF,each=nrow(Xfixed)),nrow(Xfixed),dF-1)
   
         Xfixed[,-1] <- (XfixedOrig[,-1] - minMat)/denMat
      }

      if (dR>1)
      {
         XrandomOrig <- as.matrix(Xrandom)

         minVecR <- apply(as.matrix(XrandomOrig[,-1]),2,min)
         maxVecR <- apply(as.matrix(XrandomOrig[,-1]),2,max)  
         denVecR <- maxVecR - minVecR

         if (any(denVecR==0)) stop("A non-intercept column of Xrandom has all entries identical.")
   
         minMat <- matrix(rep(minVecR,each=nrow(Xrandom)),nrow(Xrandom),dR-1)
         denMat <- matrix(rep(denVecR,each=nrow(Xrandom)),nrow(Xrandom),dR-1)

         Xrandom[,-1] <- (XrandomOrig[,-1] - minMat)/denMat
      }
   }

   # Obtain duplication matrices:

   if (dR==1)
   {
      Dd <- 1   ;   DdPlus <- 1
   }
   if (dR>1)
   {
      Dd <- duplication.matrix(dR)
      DdPlus <- solve(crossprod(Dd),t(Dd))
   }
 
   # Obtain initial values via penalised quasi-likelihood 
   # and lme4:::glmer():

   if (dF==0)
   {
      uHatLapl <- matrix(0,m,dR)
      parmVecInit <- 0
   }

   if (dF>0)
   {
      allData <- data.frame(y,Xfixed,Xrandom,idNum)
      fitLapl <- glmer(y ~ Xfixed[,-1] + (-1+Xrandom|idNum),data=allData,family=binomial(link = "probit"))
      betaLapl <- as.numeric(coefficients(summary(fitLapl))[,1])
      uHatLapl <- as.matrix(coefficients(fitLapl)$idNum)[,1:dR]
      sdVec <- attr(summary(fitLapl)$varcor$idNum,"stddev")
      corrMat <- attr(summary(fitLapl)$varcor$idNum,"correlation")
      if (dR==1) SigmaLapl <- as.numeric((sdVec^2)*corrMat)
      if (dR>1) SigmaLapl <- diag(sdVec)%*%corrMat%*%diag(sdVec)
      eigObj <- eigen(SigmaLapl)
      logSigmaLapl <- crossprod(t(eigObj$vector),(0.5*log(eigObj$values)*t(eigObj$vector)))
      thetaLapl <- vech((logSigmaLapl + t(logSigmaLapl))/2)
      parmVecInit <- c(betaLapl,thetaLapl)
   }

   # Perform search for maximum of expectation propagation approximate
   # log-likelihood:

   if (!quiet) cat("\n\n\n Starting Nelder-Mead phase:\n\n\n")

   options(warn=-1)

   optimObjNM <- optim(parmVecInit,EPlogLik,y=y,Xfixed=Xfixed,Xrandom=Xrandom,
                       dF=dF,dR=dR,m=m,nVec=nVec,numObs=numObs,indStt=indStt, 
                       uHat=uHatLapl,EPmaxit=EPmaxit,EPreltol=EPreltol,
                       method="Nelder-Mead",control=list(fnscale=(-1),
                       maxit=NMmaxit,trace=(1-as.numeric(quiet)),REPORT=1,
                       reltol=NMreltol))
   options(warn=1)
   if (!quiet)
   {
      cat("\n\n\n Finished Nelder-Mead phase.\n")
      cat("\n\n\n Starting Broyden-Fletcher-Goldfarb-Shanno phase:\n\n\n")
   }

   parmVecInit <- optimObjNM$par
   optimObjBFGS <- optim(parmVecInit,EPlogLik,y=y,Xfixed=Xfixed,Xrandom=Xrandom,
                         dF=dF,dR=dR,m=m,nVec=nVec,numObs=numObs,indStt=indStt, 
                         uHat=uHatLapl,EPmaxit=EPmaxit,EPreltol=EPreltol,method="BFGS",
                         control=list(fnscale=(-1),maxit=BFGSmaxit,
                         trace=(1-as.numeric(quiet)),REPORT=1,reltol=BFGSreltol))

   if (!quiet) cat("\n\n\n Finished Broyden-Fletcher-Goldfarb-Shanno phase.\n\n")

   # Unpack maximising parameter vector:

   if (dF==0)
   {
      betaEP <- NULL
      thetaEP <- optimObjBFGS$par
   } 
   if (dF>0)
   {
      betaEP <- optimObjBFGS$par[1:dF]
      thetaEP <- optimObjBFGS$par[-(1:dF)]
   }

   # Obtain `omega' vector corresponding to thetaEP:

   omegaEP <- as.vector(thetaTOomega(thetaEP))

   # Form parameter vector corresponding the 'omega' 
   # parameterisation:

   parmVecFORomegaEP <- c(betaEP,omegaEP)

   # Obtain Hessian with respect to (beta,omega)
   # at the maximum:

   if (!quiet) cat("\n\n\n Starting Hessian computation phase:\n\n\n")

   optimObjHess <- optim(parmVecFORomegaEP,EPlogLikFORomega,y=y,Xfixed=Xfixed,
                         Xrandom=Xrandom,dF=dF,dR=dR,m=m,nVec=nVec,numObs=numObs,
                         indStt=indStt,uHat=uHatLapl,EPmaxit=EPmaxit,EPreltol=EPreltol,
                         method="BFGS",control=list(fnscale=(-1),
                         maxit=BFGSmaxit,trace=(1-as.numeric(quiet)),REPORT=1),hessian=TRUE)

   if (!quiet) cat("\n\n\n Finished Hessian computation phase.\n\n")

   # Extract estimates and standard deviations for inference concerning
   # (beta,omega):

   betaomegaHat <- optimObjHess$par
   betaomegaSD <- sqrt(-diag(solve(optimObjHess$hessian)))

   # Obtain confidence interval limits for (beta,omega):

   cVal <- qnorm((1+confLev)/2)
   if (dF==0)
   {
      CIomegaLow <- betaomegaHat - cVal*betaomegaSD
      CIomegaUpp <- betaomegaHat + cVal*betaomegaSD
   }
   if (dF>0)
   {
      CIbetaLow <- betaomegaHat[1:dF] - cVal*betaomegaSD[1:dF]
      CIbetaUpp <- betaomegaHat[1:dF] + cVal*betaomegaSD[1:dF]
      CIomegaLow <- betaomegaHat[-(1:dF)] - cVal*betaomegaSD[-(1:dF)]
      CIomegaUpp <- betaomegaHat[-(1:dF)] + cVal*betaomegaSD[-(1:dF)]
   } 

   # Obtain `xi' vector estimates and confidence intervals:

   CIstr <- paste(100*confLev,"%"," C.I",sep="")
   if (dR==1)
   {
      xiEP <- exp(omegaEP)
      CIxiLow <- exp(CIomegaLow)
      CIxiUpp <- exp(CIomegaUpp)
 
      # Form the output matrix of inferential summaries of
      # all interpretable paramters:

      parSummMat <- matrix(NA,dF+1,3)
      if (dF>0) parSummMat[1:dF,] <- c(CIbetaLow[1:dF],betaEP[1:dF],CIbetaUpp[1:dF])
      parSummMat[dF+1,] <- c(CIxiLow,xiEP,CIxiUpp)

      parSummMat <- as.data.frame(parSummMat)
 
      dimnames1vec <- c("sigma")
      if (dF>0) dimnames1vec <- c(colnames(Xfixed),dimnames1vec)

      rownames(parSummMat) <- dimnames1vec
      colnames(parSummMat) <- c(paste(CIstr,"low"),"estimate",paste(CIstr,"upp"))
   }

   if (dR>1)
   {
      # Obtain summary matrix for model parameters:

      xiEP <- c(exp(omegaEP[1:dR]),tanh(omegaEP[-(1:dR)]))
   
      CIxiLow <- c(exp(CIomegaLow[1:dR]),tanh(CIomegaLow[-(1:dR)]))
      CIxiUpp <- c(exp(CIomegaUpp[1:dR]),tanh(CIomegaUpp[-(1:dR)]))
      
      # Form the output matrix of inferential summaries of
      # all interpretable parameters:

      if (dR==0) parSummMat <- NULL
      if (dF>0) parSummMat <- cbind(CIbetaLow[1:dF],betaEP[1:dF],CIbetaUpp[1:dF])
      
      parSummMat <- rbind(parSummMat,cbind(CIxiLow,xiEP,CIxiUpp))

      parSummMat <- as.data.frame(parSummMat)
      
      dimnames1vec <- paste("sigma",1:dR,sep="")
      for (i in 1:(dR-1))
         for (j in (i+1):dR)
            dimnames1vec <- c(dimnames1vec,paste("rho",i,j,sep=""))

      if (dF>0) dimnames1vec <- c(colnames(Xfixed),dimnames1vec)

      rownames(parSummMat) <- dimnames1vec
      colnames(parSummMat) <- c(paste(CIstr,"low"),"estimate",paste(CIstr,"upp"))
   }

   # Obtain expectation propagation-approximate best predictions for random effects:

   EPobjOpt <- EPlogLik(optimObjBFGS$par,y,Xfixed,Xrandom,dF,dR,m,nVec,numObs,
                        indStt,uHat=uHatLapl,EPmaxit,EPreltol)

   etaOut <- attr(EPobjOpt,"etaOut")

   BPmat <- matrix(NA,m,dR)
   for (i in 1:m)
   {
      Eta2Mat <- vecInverse(crossprod(DdPlus,etaOut[i,-(1:dR)])) 
      BPmat[i,] <- -0.5*solve(Eta2Mat,etaOut[i,1:dR])
   } 
   colnames(BPmat) <- colnames(Xrandom)

   if (preTransfData)  # Convert all results to correspond to the original units: 
   {
      if (dF>1)
      {
         # Save the coefficient vector corresponding to the
         # transformed data fit:

         betaStarVec <- parSummMat[1:dF,2]

         # First do the adjustment for the non-intercept fixed effect 
         # regression coefficients:

         parSummMat[2:dF,] <- parSummMat[2:dF,]/denVecF

         # Next, do the adjustment for the fixed effect intercept:

         for (j in 2:dF)
            parSummMat[1,] <- parSummMat[1,] - betaStarVec[j]*minVecF[j-1]/denVecF[j-1]
      }

      if (dR>1)
      {
         # Do the adjustment for the standard deviation parameters:

         parSummMat[(dF+2):(dF+dR),] <-  parSummMat[(dF+2):(dF+dR),]*denVecR
     
         # Do the adjustment for the random effects predictions:

         BPmat[2:dR,] <- BPmat[2:dR,]*denVecR
      } 
   }

   outObj <- list(parameters=parSummMat,randomEffects=BPmat)

   class(outObj) <- "glmmEP"

   return(outObj)
}

############ End of glmmEP ############
