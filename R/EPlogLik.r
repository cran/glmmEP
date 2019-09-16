########## R function: EPlogLik ##########

# For obtaining the expectation propagation approximate 
# log-likelihood.

# Last changed: 16 MAY 2019

EPlogLik <- function(parmVec,y,Xfixed=NULL,Xrandom,dF,dR,m,nVec,numObs,indStt,uHat=NULL,
                     EPmaxit,EPreltol)
{

   # Set the 'y dagger' response variable:

   yDagg <- 2*y - 1

   # Unpack "parmVec":

   if (dF>0) beta <- parmVec[1:dF]
   indsVechSigma <- (dF+1):length(parmVec)
   dmnVechSigma <- 0.5*dR*(dR + 1)
   parmVechLogSqrtSigma <- parmVec[indsVechSigma]
   eigObj <- eigen(vechInverse(parmVechLogSqrtSigma))
   Sigma <- crossprod(t(eigObj$vector),(exp(2*eigObj$values)*t(eigObj$vector)))

   # Obtain required duplication matrix and Moore-Penrose inverse:

   if (dR==1)
   {
      Dd <- 1 ; DdPlus <- 1
   }
   if (dR>1)
   {
      Dd <- duplication.matrix(dR)
      DdPlus <- solve(crossprod(Dd),t(Dd))
   }

   # Set up inputs for .Fortran() call to "epllk":

   idF <- dF  
   idR <- dR                  ;    idRsq <- idR^2
   lena2 <- idR + 0.5*idR*(idR-1)   ;    lena <- idR + lena2
   nMax <- max(nVec)           ;    nlena <- nMax*lena                
   
   etaSg0 <- -0.5*(idR*log(2*pi) + determinant(Sigma)$modulus)
   etaSg <- c(rep(0,idR),(-0.5*crossprod(Dd,vec(solve(Sigma)))))

   # Allocate arrays used in iterations:

   yDcur <- rep(0,nMax)
   XrCur <- matrix(0,nMax,idR)
   Xbeta <- rep(0,nMax)
   XuHat <- rep(0,nMax)

   etaFtS <- matrix(0,nMax,lena)
   etaStF <- matrix(0,nMax,lena)
   
   SUMlt <- rep(0,lena)
   c1Cur <- rep(0,idR)
   etaIN1 <- rep(0,idR)
   etaIN2 <- rep(0,lena2)
   
   eINa1 <- rep(0,idR)
   eINa2 <- rep(0,lena2)
   
   eINb1 <- rep(0,idR)
   eINb2 <- rep(0,lena2)
   
   etaPvF <- rep(0,nlena)
   etaPvS <- rep(0,nlena)
   etaCrF <- rep(0,nlena)
   etaCrS <- rep(0,nlena)
   
   # Allocate additional arrays used in call to "kpbt":
   
   wk1 <- rep(0,idRsq)
   ipvt <- rep(0,idR)
   A2ina1 <- rep(0,idR)  
   A2inc1 <- rep(0,idR)      
   A2mat <- matrix(0,idR,idR)
   A2str <- matrix(0,idR,idR)
   R2comp <- matrix(0,idR,idR)
   wk2 <- rep(0,idR)
   R5 <- matrix(0,idR,idR)
   R5TA2 <- matrix(0,idR,idR)
   vR5TA2 <- rep(0,idRsq)
   xkpan1 <- rep(0,idR)
   xkpan2 <- rep(0,lena2)
   
   # Allocate additional arrays used in call to "cpbt":
   
   wka <- rep(0,idRsq)
   wkb <- rep(0,idRsq)
   B2inb1 <- rep(0,idR)  
   det <- rep(0,2)
   work <- matrix(0,idR,idR)
   A2neg <- matrix(0,idR,idR)
   B2mat <- matrix(0,idR,idR)
   B2neg <- matrix(0,idR,idR)
   
   # Allocate additional arrays used in call to "asn":

   xm2A2 <- matrix(0,idR,idR)
   wkv <- rep(0,idRsq)
   
   # Allocate array for outputting key natural parameter
   # vectors:
   
   etaOut <- matrix(0,m,lena)
   
   # Allocate final vector of length 3 of miscellaneous quantities.
   # This `ugly' coding is due to the fact that 65 or fewer arguments
   # can be passed in a .Fortran() call.
   
   xmiscl <- c(EPmaxit,EPreltol,0)
   
   # Set up matrices for the fixed effects part, but
   # use meaningless dummy matrices containing zeroes
   # for the case where there are no fixed effects:
   
   if (dF>0)
      XfCur <- matrix(0,nMax,idF)
   
   if (dF==0)
   {
      beta <- 0
      Xfixed <- matrix(0,numObs,1)
      XfCur <- matrix(0,nMax,1)
   } 
   
   # Call "epllk" via .Fortran(): 
   
   F77obj <- .Fortran("epllk",as.double(beta),as.double(etaSg0),
                      as.double(etaSg),as.integer(m),as.integer(nVec),as.integer(nMax),
                      as.integer(numObs),as.integer(indStt),as.integer(idF),as.integer(idR),
                      as.integer(idRsq),as.integer(lena),as.integer(lena2),as.integer(nlena),
                      as.double(yDagg),as.double(Xfixed),as.double(Xrandom),as.double(yDcur),
                      as.double(XfCur),as.double(XrCur),as.double(uHat),
                      as.double(Dd),as.double(DdPlus),as.double(Xbeta),
                      as.double(XuHat),as.double(etaFtS),as.double(etaStF),
                      as.double(SUMlt),as.double(c1Cur),as.double(etaIN1),as.double(etaIN2),
                      as.double(eINa1),as.double(eINa2),as.double(eINb1),as.double(eINb2),
                      as.double(etaPvF),as.double(etaPvS),as.double(etaCrF),as.double(etaCrS),
                      as.double(wk1),as.integer(ipvt),as.double(A2ina1),as.double(A2inc1),
                      as.double(A2mat),as.double(A2str),as.double(R2comp),as.double(wk2),
                      as.double(R5),as.double(R5TA2),as.double(vR5TA2),as.double(xkpan1),
                      as.double(xkpan2),as.double(B2inb1),as.double(work),as.double(A2neg),
                      as.double(B2mat),as.double(B2neg),as.double(xm2A2),as.double(det),
                      as.double(wka),as.double(wkb),as.double(wkv),xmiscl=as.double(xmiscl),
                      etaOut=as.double(etaOut))

   # Return EP-approximate log-likelihood value for input parameters:
 
   ans <- F77obj$xmiscl[3]
   attr(ans,"etaOut") <- matrix(F77obj$etaOut,m,lena)
  
   return(ans)
}

############ End of EPlogLik ############
