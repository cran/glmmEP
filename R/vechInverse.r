########### R function: vechInverse ##########

# Obtains the inverse vecorisation-half of a vector.

# Last changed: 08 JAN 2018

vechInverse <- function(a)
{
   is.wholenumber <- function(x, tol = sqrt(.Machine$double.eps))
        return(abs(x - round(x)) < tol)

   dimA <- (sqrt(8*length(a) + 1) - 1)/2
   if (!is.wholenumber(dimA)) 
      stop("The input vector does not have a legal length.")

   out <- matrix(0,nrow=dimA,ncol=dimA)
   out[lower.tri(out, diag = TRUE)] <- a
   out <- out + t(out)
   diag(out) <- diag(out)/2
   return(out)
}

########## End of vechInverse ############


