\name{hessian}
\alias{hessian}
\title{
  Estimates the hessian matrix
}
\description{
  Given a vector of variables (x), and a function (f) that estimates one
  function value, estimates the hessian matrix by numerical differencing.
  The hessian matrix is a square matrix of second-order partial derivatives
  of the function f with respect to x. It contains, on rows i and columns j
  \deqn{d^2(f(x))/d(x_i)/d(x_j)}
}
\usage{
hessian(f, x, centered = FALSE, pert=1e-8, ...)
}
\arguments{
  \item{f }{function returning one function value, or a vector of
    function values.
  }
  \item{x }{either one value or a vector containing the x-value(s) at which
    the hessian matrix should be estimated.
  }
  \item{centered }{if TRUE, uses a centered difference approximation, else
    a forward difference approximation.
  }
  \item{pert }{numerical perturbation factor; increase depending on
    precision of model solution.
  }
  \item{... }{other arguments passed to function \code{f}.
  }
}
\value{
   The gradient matrix where the number of rows equals the length of \code{f}
   and the number of columns equals the length of \code{x}.

   the elements on i-th row and j-th column contain: \eqn{d((f(x))_i)/d(x_j)}
    
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{
## =======================================================================
## the banana function
## =======================================================================
   fun <- function(x)  100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
   mm  <-nlm(fun, p=c(0,0))$estimate
   (Hes <- hessian(fun,mm))
   # can also be estimated by nlm(fun, p=c(0,0), hessian=TRUE)
   solve(Hes)   # estimate of parameter uncertainty
}
\details{
  Function \code{hessian(f,x)} returns a forward or centered difference
  approximation of the gradient, which itself is also estimated by differencing.
  Because of that, it is not very precise.
}
\seealso{
 \code{\link{gradient}}, for a full (not necessarily square) gradient matrix

}
\keyword{math}

