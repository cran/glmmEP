########## R script: vecInverse ##########

# For obtaining the inverse vectorisation of a matrix.

# Last changed: 08 JAN 2018

vecInverse <- function(a,dmnVec=NULL)
{
   is.wholenumber <- function(x, tol = sqrt(.Machine$double.eps))
        return(abs(x - round(x)) < tol)
   
   a <- as.vector(a)

   # Check legality of inputs:
   
   if (is.null(dmnVec)) 
   {
      rtlena <- sqrt(length(a))
      if (!is.wholenumber(sqrt(length(a))))
         stop("The input vector is not of legal length.")
      dmnVec <- rep(round(rtlena),2)
   }
   if (!is.null(dmnVec))
   {    
      if (length(a)!=prod(dmnVec))
         stop("The input vector is not of legal length.")
   }

   dim(a) <- dmnVec

   return(a)
}

############ End of vecInverse ############
