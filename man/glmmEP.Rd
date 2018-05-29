\name{glmmEP}
\alias{glmmEP}
\title{Approximate frequentist inference for generalized linear mixed models via expectation propagation.}
\description{The machine learning paradigm known as expectation propagation, typically used for approximation inference in Bayesian graphical models contexts, is adapted to the problem of frequentist inference in generalized linear mixed model analysis. The essence is rapid approximation evaluation of the log-likelihood surface by using expectation propagation projections onto the Multivariate Normal family that circumvent the need for numerical integration. The approach applies to any reasonable dimension of the random effects vector, but currently is limited to probit mixed models with one level of nesting. The reference below provides full details.}
\usage{
glmmEP(y,Xfixed=NULL,Xrandom,idNum,control)
}
\arguments{
\item{y}{Vector containing the (0/1) response data.}
\item{Xfixed}{Matrix with columns consisting of predictors that enter the model as fixed effects. Each entry of the first column equal the number 1 corresponding to the fixed effects intercept.}
\item{Xrandom}{Matrix with columns consisting of predictors that enter the model as random effects. Each entry of the  first column must equal the number 1 corresponding to the random effects intercept.}
\item{idNum}{Vector containing group membership identification numbers.}
\item{control}{Function for controlling convergence and other specifications.}
}
\author{Matt Wand\email{matt.wand@uts.edu.au} and James Yu\email{james.yu@student.uts.edu.au}}
\references{Hall, P., Johnstone, I.M., Ormerod, J.T., Wand, M.P. and Yu, J. (2018). 
Fast and accurate binary response mixed model analysis via expectation propagation. 
<arXiv:1805.08423v1>.}
\examples{
library(glmmEP)

# Simulate data corresponding to a
# a simple random intercept model:

set.seed(1)  ; m <- 100 ; n <- 10
beta0True <- 0.37 ; beta1True <- 0.93 ; sigmaTrue <- 0.73
uTrue <- rnorm(m)*sigmaTrue
x <- runif(m*n)
y <- rbinom(m*n,1,pnorm(beta0True + beta1True*x 
            + crossprod(t(kronecker(diag(m),rep(1,n))),uTrue)))
idNum <- rep(1:m,each=n)
Xfixed <- cbind(1,x)
Xrandom <- matrix(1,length(y),1)
colnames(Xfixed) <- c("intercept","x")

# Obtain and summarise glmmEP() fit:

fit <- glmmEP(y,Xfixed,Xrandom,idNum)
summary(fit)

# Obtain simulated data corresponding to the
# simulation study in Section 4.1.2. of 
# Hall et al. (2018):

dataObj <- glmmSimData(seed=54321)
y <- dataObj$y  
Xfixed <- dataObj$Xfixed
Xrandom <- dataObj$Xrandom  
idNum <- dataObj$idNum

# Obtain the probit mixed model fit and summaries:

\donttest{
fitSimData <- glmmEP(y,Xfixed,Xrandom,idNum)
summary(fitSimData)
plot(fitSimData$randomEffects) ; abline(h=0) ; abline(v=0)
}

if (require("mlmRev"))
{
   # Example involving the data frame `Contraception'
   # in the package `mlmRev', starting with
   # setting up the input vectors and matrices
   # for a random intercepts model:

   data(Contraception)
   y <- as.numeric(as.character(Contraception$use)=="Y")
   Xfixed <- cbind(1,as.numeric(as.character(Contraception$urban)=="Y"),
                   Contraception$age,
                   as.numeric(as.character(Contraception$livch)=="1"),
                   as.numeric(as.character(Contraception$livch)=="2"),
                   as.numeric(as.character(Contraception$livch)=="3+"))
   colnames(Xfixed) <- c("intercept","isUrban","age","livChEq1",
                         "livChEq2","livChGe3")
   Xrandom <- as.matrix(rep(1,length(y)))
   colnames(Xrandom) <- "intercept"
   idNum <- as.numeric(as.character(Contraception$district))

   # Obtain random intercepts probit mixed model fit and summaries:

   \donttest{
   fitContracRanInt <- glmmEP(y,Xfixed,Xrandom,idNum)
   summary(fitContracRanInt)
   hist(as.numeric(fitContracRanInt$randomEffects),
        xlab="random intercepts predicted values",probability=TRUE,
        col="dodgerblue",breaks=15,main="")
   abline(v=0,col="slateblue",lwd=2)
   }

   # Change `Xrandom' to correspond to a random intercepts and 
   # slopes model:

   Xrandom <- cbind(1,as.numeric(as.character(Contraception$urban)=="Y"))
   colnames(Xrandom) <- c("intercept","isUrban")

   # Obtain random intercepts and slopes probit mixed model fit and summaries:

   \donttest{
   fitContracRanIntAndSlp <- glmmEP(y,Xfixed,Xrandom,idNum)
   summary(fitContracRanIntAndSlp)
   plot(fitContracRanIntAndSlp$randomEffects,
        col="dodgerblue",xlab="random intercepts predicted values",
        ylab="random slopes predicted values",bty="l")
   abline(v=0,col="slateblue",lwd=2); abline(h=0,col="slateblue",lwd=2)
   }
}
}
