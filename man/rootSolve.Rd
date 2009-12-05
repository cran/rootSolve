\name{rootSolve-package}
\alias{rootSolve-package}
\alias{rootSolve}
\docType{package}
\title{
  Roots and steady-states
}
\description{
Functions that:

  (1) generate gradient and Jacobian matrices (full and banded),

  (2) find roots of non-linear equations by the Newton-Raphson method,

  (3) estimate steady-state conditions of a system of (differential) equations
    in full, banded or sparse form, using the Newton-Raphson method or by
    a dynamic run,

  (4) solve the steady-state conditions for uni-and multicomponent 1-D, 2-D
    and 3-D partial differential equations, that have been converted to ODEs
    by numerical differencing (using the
    method-of-lines approach).

}
\details{
  \tabular{ll}{
     Package: \tab rootSolve\cr
     Type: \tab Package\cr
     Version: \tab 1.6\cr
     Date: \tab 2009-10-10\cr
     License: \tab GNU Public License 2 or above\cr
  }

rootSolve was created to solve the examples from chapter 7 (stability and
steady-state) from the book of Soetaert and Herman, 2009.

Please cite this work when using rootSolve.

}

\author{
  Karline Soetaert
}
\references{
  Soetaert, K and Herman, PMJ, 2009. A Practical Guide to Ecological Modelling.
  Using \R as a Simulation Platform. Springer, 372pp, ISBN 978-1-4020-8623-6.
  
  Soetaert K., 2009. rootSolve: Nonlinear root finding, equilibrium and 
  steady-state analysis of ordinary differential equations.
  R-package version 1.6
}
\seealso{

 \code{\link{uniroot.all}}, to solve for all roots of one (nonlinear) equation

 \code{\link{multiroot}}, to solve n roots of n (nonlinear) equations

 \code{\link{steady}}, for a general interface to most of the steady-state
   solvers

 \code{\link{steady.band}}, to find the steady-state of ODE models with a
   banded Jacobian

 \code{\link{steady.1D}}, \code{\link{steady.2D}},
   \code{\link{steady.3D}}, steady-state solvers for 1-D, 2-D and 3-D
   partial differential equations.

 \code{\link{stode}}, iterative steady-state solver for ODEs with full
   or banded Jacobian.

 \code{\link{stodes}}, iterative steady-state solver for ODEs with arbitrary
   sparse Jacobian.

 \code{\link{runsteady}}, steady-state solver by dynamically running to
   steady-state
 
 \code{\link{jacobian.full}}, \code{\link{jacobian.band}}, estimates the
   Jacobian matrix assuming a full or banded structure.
   
 \code{\link{gradient}}, \code{\link{hessian}}, estimates the gradient
    matrix or the Hessian.

 package vignette \code{rootSolve}
}

\examples{
\dontrun{

## run demos
demo("Jacobandroots")
demo("Steadystate")

## open the directory with documents
browseURL(paste(system.file(package="rootSolve"), "/doc", sep=""))

## main package vignette
vignette("rootSolve")
}
}

\keyword{ package }