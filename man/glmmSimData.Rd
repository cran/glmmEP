\name{glmmSimData}
\alias{glmmSimData}
\title{Simulation of data from a generalized linear mixed model.}
\description{Section 4.1.2 of the refence below descries a simulation study with data generated from a probit mixed model with six fixed effects parameters and a bivariate random effects vector having a 2 by 2 symmetric positive definite covariance matrix. The function simulates a data set from this model with 2500 groups and the number of observation in each group being a random draw from {20,21,...,30}.}
\usage{
glmmSimData(seed=12345)
}
\arguments{
\item{seed}{A positive integer which acts the seed for random data generation.}
}
\author{Matt Wand\email{matt.wand@uts.edu.au} and James Yu\email{james.yu@student.uts.edu.au}}
\references{Hall, P.,Johnstone, I.M., Ormerod, J.T., Wand, M.P. and Yu, J. (2018). 
Fast and accurate binary response mixed model analysis via expectation propagation. <arXiv:1805.08423v1>.}
\examples{
# Obtain simulated data corresponding to the simulation study in Section 4.1.2. of 
# Hall et al. (2018):

library(glmmEP)
dataObj <- glmmSimData(seed=54321)
print(names(dataObj))
}

