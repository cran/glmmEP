########## R function: glmmEP.control ##########

# Control function for glmmEP().

# Last changed: 30 JAN 2018

glmmEP.control <- function(confLev=0.95,BFGSmaxit=500,BFGSreltol=1e-10,
                           EPmaxit=100,EPreltol=1e-5,NMmaxit=100,NMreltol=1e-10,
                           quiet=FALSE,preTransfData=TRUE)
{
   if (confLev<=0) 
   {
      warning("The inputted confidence level is zero or negative.\n
               The default value of 0.95 was used instead.")
      confLev <- 0.95
   }

   if (confLev>=1) 
   {
      warning("The inputted confidence level is 1 or higher.\n
               The default value of 0.95 was used instead.")
      confLev <- 0.95
   }

   if (BFGSmaxit<=0) 
   {
      warning("The inputted Broyden-Fletcher-Goldfarb-Shanno maximum number of iterations value is 0 or negative.\n
               The default value of 100 was used instead.")
      BFGSmaxit <- 500
   }

   if (BFGSreltol<=0) 
   {
      warning("The inputted Broyden-Fletcher-Goldfarb-Shanno relative tolerance value is 0 or negative.\n
               The default value of 1e-10 was used instead.")
      BFGSreltol <- 1e-10
   }

   if (EPmaxit<=0) 
   {
      warning("The inputted expectation propagation maximum number of iterations value is 0 or negative.\n
               The default value of 100 was used instead.")
      EPSmaxit <- 100
   }

   if (EPreltol<=0) 
   {
      warning("The inputted expectation propagation relative tolerance value is 0 or negative.\n
               The default value of 1e-5 was used instead.")
      EPSreltol <- 1e-5
   }

   if (NMmaxit<=0) 
   {
      warning("The inputted Nelder-Mead maximum number of iterations value is 0 or negative.\n
               The default value of 100 was used instead.")
      NMmaxit <- 100
   }

   if (NMreltol<=0) 
   {
      warning("The inputted Nelder-Mead relative tolerance value is 0 or negative.\n
               The default value of 1e-10 was used instead.")
      NMreltol <- 1e-10
   }

   if ((quiet!=TRUE)&(quiet!=FALSE))
   {
      warning("The inputted quiet value is non-Boolean. The default value of FALSE was used instead.\n")
   }

   if ((preTransfData!=TRUE)&(preTransfData!=FALSE))
   {
      warning("The inputted preTransfData value is non-Boolean. The default value of FALSE was used instead.\n")
   }

   return(list(confLev=confLev,BFGSmaxit=BFGSmaxit,BFGSreltol=BFGSreltol,EPmaxit=EPmaxit,
               EPreltol=EPreltol,NMmaxit=NMmaxit,NMreltol=NMreltol,quiet=quiet,
               preTransfData=preTransfData))
}

############ End of glmmEP.control ############
