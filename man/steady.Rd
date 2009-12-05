\name{steady}
\alias{steady}
\title{
  General steady-state solver for a set of ordinary differential equations.
}
\description{
  Estimates the steady-state condition for a system of ordinary differential
  equations.

  This is a wrapper around steady-state solvers \code{stode}, \code{stodes}
  and \code{runsteady}.
}
\usage{
steady(y, time=0, func, parms=NULL, method="stode",...)}

\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector.
    If \code{y} has a name attribute, the names will be used to label the output matrix.
  }
  \item{time }{time for which steady-state is wanted; the default is time=0.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ode system (the model defininition) at time \code{time},
    or a character string giving the name of a compiled function in a
    dynamically loaded shared library.
    If \code{func}  is an \R-function, it must be defined as:
    \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
    in the integration, \code{y} is the current estimate of the variables
    in the ODE system.  If the initial values \code{y} has a names
    attribute, the names will be available inside \code{func}. \code{parms} is
    a vector or list of parameters; ... (optional) are any other arguments
    passed to the function.
  
    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{time}, and whose next elements are global values whose steady-state
    value is also required.
    
    The derivatives
    should be specified in the same order as the state variables \code{y}.
  }
  \item{parms }{parameters passed to \code{func}.
  }
  \item{method }{the solution method to use, one of \code{stode}, \code{stodes}
    or \code{runsteady}.
  }
  \item{... }{additional arguments passed to function \code{stode},
    \code{stodes} or \code{runsteady}.
  }
}
\value{
  A list containing 
  \item{y }{ a vector with the state variable values from the last iteration
    during estimation of steady-state condition of the system of equations.
    If \code{y} has a names attribute, it will be used to label the output values. }
  \item{... }{ the number of "global" values returned.
  }
  The output will have the attribute \code{steady}, which returns \code{TRUE},
  if steady-state has been reached and the attribute
  \code{precis} with the precision attained during each iteration.
}
\examples{
## =======================================================================
##  Bacteria (Bac) growing on a substrate (Sub)
## =======================================================================
  
  model <- function(t,state,pars)
  {
  with (as.list(c(state,pars)), {
  #       substrate uptake             death  respiration
  dBact = gmax*eff*Sub/(Sub+ks)*Bact - dB*Bact - rB*Bact
  dSub  =-gmax    *Sub/(Sub+ks)*Bact + dB*Bact          +input
  
  return(list(c(dBact,dSub)))
                                })
  }
  
  pars <- list(gmax =0.5,eff = 0.5,
               ks =0.5, rB =0.01, dB =0.01, input=0.1)
  # Newton-Raphson
  steady(y=c(Bact=0.1,Sub=0),time=0,
         func=model,parms=pars,pos=TRUE)  

  # Dynamic run to steady-state
  as.data.frame(steady(y=c(Bact=0.1,Sub=0),time=c(0,1e5),
         func=model,parms=pars,method="runsteady"))

}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\details{
  This is simply a wrapper around the various steady-state solvers.

  See package vignette for information about specifying the
  model in compiled code.

  See the selected solver for the additional options.
}

\seealso{
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

}

\keyword{math}

