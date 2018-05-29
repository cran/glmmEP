########## R function: vech ##########

# Obtains the vectorisation-half of a matrix.

# Last changed: 08 JAN 2018

vech <- function(A)
{
   if (nrow(A) != ncol(A))
      stop("A must be a square matrix")
   return(A[row(A) >= col(A)])
}

############ End of vech ############
