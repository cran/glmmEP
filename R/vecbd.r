########## R script: vecbd ##########

# For evaluation of vecbd(A) for a square
# matrix A - the vector of entries of A
# below the diagonal.

# Last changed: 05 JAN 2018

vecbd <- function(A)
{
   # Check the legality of the input matrix:

   if (nrow(A)!=ncol(A)) stop("A must be a square matrix.\n")
   if (nrow(A)==1) stop("A must have at least 2 rows.\n")
   
   diag(A) <- NA
   
   clean <- function(x)
       return(x[is.na(x)==0])

   return(clean(vech(A)))
}

############ End of vecbd ############
