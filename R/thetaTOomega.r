########## R script: thetaTOomega ##########

# For conversion from a `theta' vector to 
# its corresponding `omega' vector.

# Last changed: 08 JAN 2018

thetaTOomega <- function(theta)
{
   # Determine dimension value:

   dmn <- (sqrt(8*length(theta)+1) - 1)/2
   is.wholenumber <- function(x,tol=sqrt(.Machine$double.eps))
      return(abs(x-round(x))<tol)
   if (!is.wholenumber(dmn))
      stop("The input vector is not of legal length.")

   # Obtain spectral decomposition of vech^{-1}(theta):
   
   eigObj <- eigen(vechInverse(theta))

   # Obtain corresponding `Sigma' matrix:

   Sigma <- crossprod(t(eigObj$vector),(exp(2*eigObj$values)*t(eigObj$vector)))

   # Obtain `omega' vector:

   if (dmn==1) omega <- log(sqrt(Sigma[1,1]))

   if (dmn>1) 
   {
      dgSigma <- diag(Sigma)
      omega <- c(log(sqrt(dgSigma)),atanh(vecbd(Sigma)/sqrt(vecbd(tcrossprod(dgSigma)))))
   }

   return(as.numeric(omega))
}

############ End of thetaTOomega ############
