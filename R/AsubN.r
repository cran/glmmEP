########## R function: AsubN ##########

# For computation of the log-partition function
# of the Multivariate Normal distribution.

# Last changed: 18 JAN 2018

AsubN <- function(a)
{
   dmn <- (sqrt(8*length(a) + 9) - 3)/2
   if (dmn==1)
   {
      Dd <- 1 ; DdPlus <- 1
   }
   if (dmn>1)
   {
      Dd <- duplication.matrix(dmn)
      DdPlus <- solve(crossprod(Dd),t(Dd))
   }

   a1 <- a[1:dmn] ; a2 <- a[-(1:dmn)]
   A2 <- as.matrix(vecInverse(crossprod(DdPlus,a2)))

   return(-0.25*crossprod(a1,solve(A2,a1)) - 0.5*determinant((-2*A2))$modulus)
}

############ End of AsubN ############
