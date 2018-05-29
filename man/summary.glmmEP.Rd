\name{summary.glmmEP}
\alias{summary.glmmEP}
\title{Inferential summary for a generalized linear mixed model with expectation propagation fitting.}
\description{Produces a table containing approximated maximum likelihood estimates and corresponding confidence interval limits for the fixed effect parameters and random effects standard deviation and correlation parameters.}
\usage{
\method{summary}{glmmEP}(object,...)
}
\arguments{
\item{object}{A \code{glmmEP()} fit object.}
\item{...}{Place-holder for additional arguments.}
}
\author{Matt Wand\email{matt.wand@uts.edu.au} and James Yu\email{jamescfyu@gmail.com}}
\examples{
dataObj <- glmmSimData(seed=54321)
y <- dataObj$y  
Xfixed <- dataObj$Xfixed
Xrandom <- dataObj$Xrandom  
idNum <- dataObj$idNum
\donttest{
fitSimData <- glmmEP(y,Xfixed,Xrandom,idNum)
summary(fitSimData)
}
}
