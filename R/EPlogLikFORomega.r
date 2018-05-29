########## R function: EPlogLikFORomega ##########

# The function EPlogLik() with `omega' rather than
# 'theta' inputted for the covariance matrix component.

# Last changed: 11 JAN 2018

EPlogLikFORomega <- function(parmVecFORomega,y,Xfixed,Xrandom,dF,dR,
                             m,nVec,numObs,indStt,uHat,EPmaxit,EPreltol)
{
   # Unpack "parmVecFORomega":

   if (dF==0) 
   {
      beta <- NULL
      omega <- parmVecFORomega
   }
   if (dF>0) 
   {
      beta <- parmVecFORomega[1:dF]
      omega <- parmVecFORomega[-(1:dF)]
   }

   # Convert `omega' vector to its corresponding `theta'
   # vector:

   theta <- omegaTOtheta(omega)

   # Construct "parmVecFORtheta":

   parmVecFORtheta <- c(beta,theta)
  
   # Return expectation approximate log-likelihood:

   return(EPlogLik(parmVecFORtheta,y,Xfixed,Xrandom,dF,dR,m,nVec,numObs,
                   indStt,uHat,EPmaxit,EPreltol))
}

############ End of EPlogLikFORomega ############
