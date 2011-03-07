\name{jacobian.full}
\alias{jacobian.full}
\title{
  Full square jacobian matrix for a system of ODEs (ordinary differential
  equations)
}
\description{
  Given a vector of (state) variables, and a function that estimates one
  function value for each (state) variable (e.g. the rate of change),
  estimates the Jacobian matrix (\eqn{d(f(x))/d(x)})

  Assumes a full and square Jacobian matrix
}

\usage{
jacobian.full(y, func, dy = NULL, time = 0, parms = NULL, 
              pert = 1e-8, ...)
}

\arguments{
  \item{y }{(state) variables, a vector; if \code{y} has a name attribute,
    the names will be used to label the Jacobian matrix columns.
  }
  \item{func }{function that calculates one function value for each element
    of \code{y}; if an ODE system, \code{func} calculates the rate of change
    (see details).
  }
  \item{dy }{reference function value; if not specified, it will be estimated
    by calling \code{func}.
  }
  \item{time }{time, passed to function \code{func}.
  }
  \item{parms }{parameter values, passed to function \code{func}.
  }
  \item{pert }{numerical perturbation factor; increase depending on precision
    of model solution.
  }
  \item{... }{other arguments passed to function \code{func}.
  }
}
\value{
  The square jacobian matrix;  the elements on i-th row and j-th column are
  given by: \eqn{d(f(x)_i)/d(x_j)}
}

\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{
## =======================================================================
## 1. Structure of the Jacobian
## =======================================================================
mod <- function (t = 0, y, parms = NULL,...)
{
 dy1<-  y[1] + 2*y[2]
 dy2<-3*y[1] + 4*y[2] + 5*y[3]
 dy3<-         6*y[2] + 7*y[3] + 8*y[4]
 dy4<-                  9*y[3] +10*y[4]
 return(as.list(c(dy1, dy2, dy3, dy4)))
}

jacobian.full(y = c(1, 2, 3, 4), func = mod)

## =======================================================================
## 2. Stability properties of a physical model
## =======================================================================
coriolis <- function (t, velocity, pars, f)
{
  dvelx <- f*velocity[2]
  dvely <- -f*velocity[1]
  list(c(dvelx, dvely))
}

# neutral stability; f is coriolis parameter
Jac <- jacobian.full(y = c(velx = 0, vely = 0), func = coriolis,
                     parms = NULL, f = 1e-4)
print(Jac)
eigen(Jac)$values

## =======================================================================
## 3. Type of equilibrium
## =======================================================================
## From Soetaert and Herman (2009). A practical guide to ecological 
## modelling. Using R as a simulation platform. Springer

eqn <- function (t, state, pars)
 {
  with (as.list(c(state, pars)),  {
  dx <- a*x + cc*y
  dy <- b*y + dd*x
  list(c(dx, dy))
  })
 }

# stable equilibrium
A <- eigen(jacobian.full(y = c(x = 0, y = 0), func = eqn,
              parms = c(a = -0.1, b = -0.3, cc = 0, dd = 0)))$values
# unstable equilibrium
B <- eigen(jacobian.full(y = c(x = 0, y = 0), func = eqn,
              parms = c(a = 0.2, b = 0.2, cc = 0.0, dd = 0.2)))$values
# saddle point
C <- eigen(jacobian.full(y = c(x = 0, y = 0), func = eqn,
              parms = c(a = -0.1, b = 0.1, cc = 0, dd = 0)))$values
# neutral stability
D <- eigen(jacobian.full(y = c(x = 0, y = 0), func = eqn,
              parms = c(a = 0, b = 0, cc = -0.1, dd = 0.1)))$values
# stable focal point
E <- eigen(jacobian.full(y = c(x = 0, y = 0), func = eqn,
              parms = c(a = 0, b = -0.1, cc = -0.1, dd = 0.1)))$values
# unstable focal point
F <- eigen(jacobian.full(y = c(x = 0, y = 0), func=eqn,
              parms = c(a = 0., b = 0.1, cc = 0.1, dd = -0.1)))$values

data.frame(type = c("stable", "unstable", "saddle", "neutral",
           "stable focus", "unstable focus"),
           eigenvalue_1 = c(A[1], B[1], C[1], D[1], E[1], F[1]),
           eigenvalue_2 = c(A[2], B[2], C[2], D[2], E[2], F[2]))

## =======================================================================
## 4. Limit cycles
## =======================================================================
## From Soetaert and Herman (2009). A practical guide to ecological 
## modelling. Using R as a simulation platform. Springer

eqn2 <- function (t, state, pars)
 {
  with (as.list(c(state, pars)),
  {
  dx<-  a*y   + e*x*(x^2+y^2-1)
  dy<-  b*x   + f*y*(x^2+y^2-1)
  list(c(dx, dy))
  })
 }

# stable limit cycle with unstable focus
eigen(jacobian.full(c(x = 0, y = 0), eqn2,
                    parms = c(a = -1, b = 1, e = -1, f = -1)))$values
# unstable limit cycle with stable focus
eigen(jacobian.full(c(x = 0, y = 0), eqn2,
                    parms = c(a = -1, b = 1, e = 1, f = 1)))$values
}

\details{
  The function \code{func} that estimates the rate of change of the state
  variables has to be consistent with functions called from
	R-package \code{deSolve}, which contains integration routines.
	
  This function call is as: \bold{function(time,y,parms,...)} where
  \itemize{
	  \item \code{y} : (state) variable values at which the Jacobian is estimated.
  	\item \code{parms}: parameter vector  - need not be used.
  	\item \code{time}: time at which the Jacobian is estimated - in general,
      \code{time} will not be used.
  	\item \code{...}: (optional) any other arguments.
  }
  The Jacobian is estimated numerically, by perturbing the x-values.
}

\seealso{
 \code{\link{jacobian.band}}, estimates the Jacobian matrix
   assuming a banded structure.

 \code{\link{hessian}}, estimates the Hessian matrix.

 \code{\link{gradient}}, for a full (not necessarily square) gradient matrix
   and where the function call is simpler.

 \code{\link{uniroot.all}}, to solve for all roots of one (nonlinear) equation

 \code{\link{multiroot}}, to solve n roots of n (nonlinear) equations

}

\note{
  This function is useful for stability analysis of ODEs, which start by
  estimating the Jacobian at equilibrium points.
  The type of equilibrium then depends on the eigenvalue of the Jacobian.
}
\keyword{math}

