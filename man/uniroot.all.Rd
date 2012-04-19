\name{uniroot.all}
\alias{uniroot.all}
\title{
  Finds many (all) roots of one equation within an interval
}
\description{
  The function \code{uniroot.all} searches the interval from lower to upper
  for several roots (i.e., zero's) of a function \code{f} with respect to
  its first argument.
}
\usage{
uniroot.all(f, interval, lower = min(interval), upper = max(interval), 
            tol = .Machine$double.eps^0.2, maxiter = 1000, n = 100, ...)}
\arguments{
  \item{f }{the function for which the root is sought.
  }
  \item{interval }{a vector containing the end-points of the interval to
    be searched for the root.
  }
  \item{lower }{the lower end point of the interval to be searched.
  }
  \item{upper }{the upper end point of the interval to be searched.
  }
  \item{tol }{the desired accuracy (convergence tolerance).
  }
  \item{maxiter }{the maximum number of iterations.
  }
  \item{n }{number of subintervals in which the root is sought.
  }
  \item{... }{additional named or unnamed arguments to be passed to \code{f}
    (but beware of partial matching to other arguments).
  }
}
\value{
  a vector with the roots found in the interval
}
\author{
Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
## =======================================================================
##   Mathematical examples  
## =======================================================================

# a well-behaved case...
fun <- function (x) cos(2*x)^3

curve(fun(x), 0, 10,main = "uniroot.all")

All <- uniroot.all(fun, c(0, 10))
points(All, y = rep(0, length(All)), pch = 16, cex = 2)

# a difficult case...
f <- function (x) 1/cos(1+x^2)
AA <- uniroot.all(f, c(-5, 5))
curve(f(x), -5, 5, n = 500, main = "uniroot.all")
points(AA, rep(0, length(AA)), col = "red", pch = 16)

f(AA)  # !!!

## =======================================================================
## Ecological modelling example  
## =======================================================================

# Example from the book of Soetaert and Herman(2009)
# A practical guide to ecological modelling -
# using R as a simulation platform. Springer

r   <- 0.05
K   <- 10
bet <- 0.1
alf <- 1

# the model : density-dependent growth and sigmoid-type mortality rate
rate <- function(x, r = 0.05) r*x*(1-x/K) - bet*x^2/(x^2+alf^2)

# find all roots within the interval [0,10]
Eq   <- uniroot.all(rate, c(0, 10))

# jacobian evaluated at all roots: 
# This is just one value - and therefore jacobian = eigenvalue
# the sign of eigenvalue: stability of the root: neg=stable, 0=saddle, pos=unstable

eig <- vector()
for (i in 1:length(Eq)) 
  eig[i] <- sign (gradient(rate, Eq[i]))

curve(rate(x), ylab = "dx/dt", from = 0, to = 10,
      main = "Budworm model, roots", 
      sub = "Example from Soetaert and Herman, 2009")
abline(h = 0)
points(x = Eq, y = rep(0, length(Eq)), pch = 21, cex = 2,
       bg = c("grey", "black", "white")[eig+2] )
legend("topleft", pch = 22, pt.cex = 2,
       c("stable", "saddle", "unstable"),
       col = c("grey", "black", "white"), 
       pt.bg = c("grey", "black", "white"))}
\details{
  \code{f} will be called as \code{f(x, ...)} for a numeric value of \code{x}.

  Run \code{demo(Jacobandroots)} for an example of the use of \code{uniroot.all}
  for steady-state analysis.

  See also second example of \code{\link{gradient}}
  This example is discussed in the book by Soetaert and Herman (2009).
}
\seealso{
  \code{\link{uniroot}} for more information about input.
}
\note{
  The function calls \code{uniroot}, the basic R-function.

  It is not guaranteed that all roots will be recovered.

  This will depend on \code{n}, the number of subintervals in which the
  interval is divided.

  If the function "touches" the X-axis (i.e. the root is a saddle point),
  then this root will generally not be retrieved.
  (but chances of this are pretty small).
}
\keyword{optimize}

