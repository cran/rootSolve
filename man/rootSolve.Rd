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

  (3) estimate steady-state conditions of a system of (differential) equations in full, banded or sparse form, using the Newton-Raphson method or by a dynamic run,

  (4) solve the steady-state conditions for uni-and multicomponent 1-D and 2-D reactive transport models (boundary value problems of ODE) using the method-of-lines approach.

}
\details{
  \tabular{ll}{
    Package: \tab rootSolve\cr
    Type: \tab Package\cr
    Version: \tab 1.3\cr
    Date: \tab 2009-01-31\cr
    License: \tab GNU Public License 2 or above\cr
  }

rootSolve is designed for solving n roots of n nonlinear equations.


}
\author{
  Karline Soetaert
}
%\references{
%% to be filled
%}
\seealso{

 \code{\link{uniroot.all}}, to solve all roots of one equation

 \code{\link{multiroot}}, to solve n roots of n equations

 \code{\link{steady}}, \code{\link{steady.1D}}, \code{\link{steady.2D}} general steady-state solvers

 \code{\link{stode}}, \code{\link{stodes}}, steady-state solvers for full, banded or arbitrary sparse models (Newton-Raphson method)
 
 \code{\link{runsteady}} steady-state solver by dynamically running to steady-state
 
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