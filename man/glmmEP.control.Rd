\name{glmmEP.control}
\alias{glmmEP.control}
\title{Controlling generalized linear mixed model fitting via expectation propagation}
\description{Function for optional use in calls to \code{glmmEP()} to control convergence values and other specifications for expectation propagation-based fitting of generalized linear mixed models.}
\usage{
glmmEP.control(confLev=0.95,BFGSmaxit=500,BFGSreltol=1e-10,
               EPmaxit=100,EPreltol=1e-5,NMmaxit=100,NMreltol=1e-10,
               quiet=FALSE,preTransfData=TRUE)
}
\arguments{
\item{confLev}{Confidence level of confidence intervals expressed as a proportion (i.e. a number between 0 and 1). The default is 0.95 corresponding to 95\% confidence intervals.}
\item{BFGSmaxit}{Positive integer specifying the maximum number of iterations in the Broyden-Fletcher-Goldfarb-Shanno optimization phase. The default is 500.}
\item{BFGSreltol}{Positive number specifying the relative tolerance value as defined in the R function optim() in the Broyden-Fletcher-Goldfarb-Shanno optimization phase. The default is 1e-10.}
\item{EPmaxit}{Positive integer specifying the maximum number of iterations in the expectation propagation message passing iterations. The default is 100.}
\item{EPreltol}{Positive number specifying the relative tolerance value for the expectation propagation message passing iterations. The default is 1e-5.}
\item{NMmaxit}{Positive integer specifying the maximum number of iterations in the Nelder-Mead optimization phase. The default is 100.}
\item{NMreltol}{Positive number specifying the relative tolerance value as defined in the R function optim() in the Nelder-Mead optimization phase. The default is 1e-10.}
\item{quiet}{Flag for specifying whether or not glmmEP() runs `quietly' or with progress reports printed to the screen. The default is FALSE.}
\item{preTransfData}{Flag for specifying whether or not the predictor data are pre-transformed to the unit interval for fitting, with estimates, predictions and confidence intervals transformed to match the units of the original data before. The default is TRUE.}
}
\author{Matt Wand\email{matt.wand@uts.edu.au} and James Yu\email{jamescfyu@gmail.com}}
\references{Hall, P., Johnstone, I.M., Ormerod, J.T., Wand, M.P. and Yu, J. (2018). 
Fast and accurate binary response mixed model analysis via expectation propagation. Submitted.}
\examples{
library(glmmEP)

# Obtain simulated data corresponding to the simulation study in Section 4.1.2. of 
# Hall et al. (2018):

dataObj <- glmmSimData(seed=54321)
y <- dataObj$y  
Xfixed <- dataObj$Xfixed
Xrandom <- dataObj$Xrandom  
idNum <- dataObj$idNum

# Obtain and summarise probit mixed model fit with user control
# of some of the parameters in glmmEP.control():

myNMmaxit <- 500 ; myBFGSreltol <- 0.001
\donttest{
fitSimData <- glmmEP(y,Xfixed,Xrandom,idNum,
              control=glmmEP.control(NMmaxit=myNMmaxit,BFGSreltol=myBFGSreltol,quiet=TRUE))
summary(fitSimData)}
}



