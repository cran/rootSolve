\name{jacobian.band}
\alias{jacobian.band}
\title{
  Banded jacobian matrix for a system of ODEs (ordinary differential equations)
}
\description{
  Given a vector of (state) variables \code{y}, and a function that estimates a
  function value for each (state) variable (e.g. the rate of change),
  estimates the Jacobian matrix (\eqn{d(f(y))/d(y)}).

  Assumes a banded structure of the Jacobian matrix, i.e. where the non-zero
  elements are restricted to a number of bands above and below the diagonal.
}

\usage{
jacobian.band(y, func, bandup=1, banddown=1, 
  dy=NULL, time=0, parms=NULL, pert=1e-8, ...)
}

\arguments{
  \item{y }{(state) variables, a vector; if \code{y} has a name attribute,
    the names will be used to label the jacobian matrix columns.
  }
  \item{func }{function that calculates one function value for each element
    of \code{y}; if an ODE system, \code{func} calculates the rate of change
    (see details).
  }
  \item{bandup }{number of nonzero bands above the diagonal of the Jacobian
    matrix.
  }
  \item{banddown }{number of nonzero bands below the diagonal of the Jacobian
    matrix.
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
  Jacobian matrix, in banded format, i.e. only the nonzero bands near the
  diagonal form the rows of the Jacobian. 

  this matrix has \code{bandup}+\code{banddown}+1 rows, while the number of
  columns equal the length of \code{y}.

  Thus, if the full Jacobian is given by:
  \tabular{lllll}{
         \tab [,1],\tab [,2],\tab [,3],\tab [,4] \cr
    [,1] \tab 1    \tab 2    \tab 0    \tab 0 \cr
    [,2] \tab 3    \tab 4    \tab 5    \tab 0 \cr
    [,3] \tab 0    \tab 6    \tab 7    \tab 8 \cr
    [,4] \tab 0    \tab 0    \tab 9    \tab 10
  }
		
  the banded jacobian will be:
  \tabular{rllll}{
         \tab [,1],\tab [,2],\tab [,3],\tab [,4] \cr
    [,1] \tab 0    \tab 2    \tab 5    \tab 8 \cr
    [,2] \tab 1    \tab 4    \tab 7    \tab 10 \cr
    [,3] \tab 3    \tab 6    \tab 9    \tab 0
  }
}

\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{mod <- function (t=0,y, parms=NULL,...)
{
 dy1<-  y[1] + 2*y[2]
 dy2<-3*y[1] + 4*y[2] + 5*y[3]
 dy3<-         6*y[2] + 7*y[3] + 8*y[4]
 dy4<-                  9*y[3] +10*y[4]
 return(as.list(c(dy1,dy2,dy3,dy4)))
}

jacobian.band(y=c(1,2,3,4),func=mod)
}

\details{
  The function \code{func} that estimates the rate of change of the state
  variables has to be	consistent with functions called from R-package
  \code{deSolve}, which contains integration routines.
	
  This function call is as: \bold{function(time,y,parms,...)} where
  \itemize{
  	\item \code{y} : (state) variable values at which the Jacobian is estimated.
	  \item \code{parms}: parameter vector  - need not be used.
	  \item \code{time}: time at which the Jacobian is estimated - in general,
      \code{time} will not be used.
	  \item \code{...}: (optional) any other arguments
  }
  The Jacobian is estimated numerically, by perturbing the x-values.
}

\seealso{
 \code{\link{jacobian.full}}, estimates the Jacobian matrix
   assuming a full matrix.

 \code{\link{hessian}}, estimates the Hessian matrix.

 \code{\link{gradient}}, for a full (not necessarily square) gradient matrix
   and where the function call is simpler.

 \code{\link{uniroot.all}}, to solve for all roots of one (nonlinear) equation

 \code{\link{multiroot}}, to solve n roots of n (nonlinear) equations

}
\keyword{math}

