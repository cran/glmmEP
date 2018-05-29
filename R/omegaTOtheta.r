########## R script: omegaTOtheta ##########

# For conversion from an `omega' vector to 
# its corresponding `theta' vector.

# Last changed: 05 JAN 2018

omegaTOtheta <- function(omega)
{
   # Determine dimension value:

   dmn <- (sqrt(8*length(omega)+1) - 1)/2
   is.wholenumber <- function(x,tol=sqrt(.Machine$double.eps))
      return(abs(x-round(x))<tol)
   if (!is.wholenumber(dmn))
      stop("input vector not of legal length")

   # Obtain the `Sigma' matrix:

   if (dmn==1) Sigma <- exp(2*omega)
   
   if (dmn>1)
   {
      Sigma <- matrix(NA,dmn,dmn)

      omega1 <- omega[1:dmn]   ;    omega2 <- omega[-(1:dmn)] 

      diag(Sigma) <- exp(2*omega1)

      vecbdSigma <- tanh(omega2)*vecbd(tcrossprod(exp(omega1)))

      iPos <- 1 
      for (j in 1:(dmn-1))
         for (i in (j+1):dmn)
         { 
            Sigma[i,j] <- vecbdSigma[iPos]
            Sigma[j,i] <- vecbdSigma[iPos]
            iPos <- iPos + 1
         }
   }

   # Obtain spectral decomposition of Sigma:
   
   eigObj <- eigen(Sigma)

   # Obtain corresponding `theta' vector:

   theta <- vech(crossprod(t(eigObj$vector),(0.5*log(eigObj$values)*t(eigObj$vector))))

   return(as.numeric(theta))
}

############ End of omegaTOtheta ############
