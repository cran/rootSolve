\name{gradient}
\alias{gradient}
\title{Estimates the gradient matrix for a simple function}
\description{Given a vector of variables (x), and a function (f) that estimates one function value or a
 set of function values (\eqn{f(x)}), estimates the gadient matrix, containing, on rows i and columns j
 \deqn{d(f(x)_i)/d(x_j)}
 \cr The gradient matrix is not necessarily square}
\usage{gradient(f, x, ...)}
\arguments{
  \item{f }{function returning one function value, or a vector of function values}
  \item{x }{either one value or a vector containing the x-value(s) at which the gradient matrix should be estimated}
  \item{... }{other arguments passed to function \code{f}}
}
\value{
   The gradient matrix where the number of rows equals the length of \code{f} and the number of columns equals the length of \code{x}.
    \cr the elements on i-th row and j-th column contain: \eqn{d((f(x))_i)/d(x_j)}
    
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{############################################################
# 1. Sensitivity analysis of the logistic differential equation 
# dN/dt = r*(1-N/K)*N  , N(t0)=N0.

# analytical solution of the logistic equation:
logistic <- function (x,times)
{
 with (as.list(x),
 {
  N=K/(1+(K-N0)/N0*exp(-r*times))
  return(N)
  })
}

# parameters for the US population from 1900
x=c(N0=76.1,r=0.02,K=500)

# Sensitivity function: SF: dfi/dxj at 
# output intervals from 1900 to 1950
SF<-gradient(f=logistic,x,times=0:50)

# mean sensitivity, scaled for the value of the parameter: 
# mean [dfi/(dxj/xj)]= mean (SF*x)
rowMeans(t(SF)*x)

############################################################
# 2. Stability of the budworm model, as a function of its
# rate of increase.
#
# Example from the book of Soetaert and Herman(2008)
# A practical guide to ecological modelling
# using R as a simulation platform. Springer
# code and theory are explained in this book

r   <- 0.05
K   <- 10
bet <- 0.1
alf <- 1

# density-dependent growth and sigmoid-type mortality rate
rate <- function(x,r=0.05) r*x*(1-x/K)-bet*x^2/(x^2+alf^2)

# Stability of a root ~ sign of eigenvalue of Jacobian 
stability <- function (r)
{
  Eq <- uniroot.all(rate,c(0,10),r=r)
  eig  <- vector()
  for (i in 1:length(Eq)) 
      eig[i] <- sign (gradient(rate,Eq[i],r=r))
  return(list(Eq=Eq,Eigen=eig))
}

# bifurcation diagram
rseq <- seq(0.01,0.07,by=0.0001)

plot(0,xlim=range(rseq),ylim=c(0,10),type="n",
     xlab="r",ylab="B*",main="Budworm model, bifurcation",
     sub="Example from book of Soetaert and Herman")

for (r in rseq) {
  st <- stability(r) 
  points(rep(r,length(st$Eq)),st$Eq,pch=22,
         col=c("darkblue","black","lightblue")[st$Eigen+2],
         bg =c("darkblue","black","lightblue")[st$Eigen+2]) 
}

legend("topleft",pch=22,pt.cex=2,c("stable","unstable"),
col=c("darkblue","lightblue"),pt.bg=c("darkblue","lightblue"))}
\details{the function \code{f} that estimates the function values will be called as  
  f(x, ...). If \code{x} is a vector, then the first argument passed to \code{f} should also be a vector. 
  
  The gradient is estimated numerically, by perturbing the x-values.}
\seealso{\code{\link{jacobian.full}}, for generating a full and \bold{square} gradient (jacobian) matrix and where the function call is more complex}
\note{\code{gradient} can be used to calculate so-called \code{sensitivity functions}, 
  that estimate the effect of parameters on output variables.}
\keyword{math}

