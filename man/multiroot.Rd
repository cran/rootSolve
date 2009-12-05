\name{multiroot}
\alias{multiroot}
\title{
  Solves for n roots of n (nonlinear) equations.
}
\description{
  Given a vector of n variables, and a set of n (nonlinear) equations in
  these variables,

  estimates the root of the equations, i.e. the variable values where all
  function values = 0.

  Assumes a full Jacobian matrix, uses the Newton-Raphson method.
}

\usage{
multiroot(f, start, maxiter=100, 
  rtol=1e-6, atol=1e-8, ctol=1e-8, 
  useFortran=TRUE, positive=FALSE,
  jacfunc=NULL, jactype="fullint",
  verbose=FALSE, bandup=1, banddown=1, ...)
}
\arguments{
  \item{f }{function for which the root is sought; it must return a vector
    with as many values as the length of \code{start}.
  }
  \item{start }{vector containing initial guesses for the unknown x;
    if \code{start} has a name attribute, the names will be used to label
    the output vector.
  }
  \item{maxiter }{maximal number of iterations allowed.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector, one
    value for each element in the unknown x.
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one
    value for each element in x.
  }
  \item{ctol }{a scalar. If between two iterations, the maximal change in
    the variable values is less than this amount, then it is assumed that
    the root is found.
  }
  \item{useFortran }{logical, if \code{FALSE}, then an \R -implementation
    of the Newton-Raphson method is used - see details.
  }
  \item{positive}{if \code{TRUE}, the unknowns y are forced to be
    non-negative numbers.
  }
  \item{jacfunc }{if not \code{NULL}, a user-supplied \R function that
    estimates the Jacobian of the system of differential equations
    dydot(i)/dy(j). In some circumstances, supplying \code{jacfunc}
    can speed up the computations.  The \R calling sequence for
    \code{jacfunc} is identical to that of \code{func}.

    If the Jacobian is a full matrix, \code{jacfunc} should return a matrix
    dydot/dy, where the ith row contains the derivative of \eqn{dy_i/dt}
    with respect to \eqn{y_j}, or a vector containing the matrix elements
    by columns.

    If the Jacobian is banded, \code{jacfunc} should return a matrix containing
    only the nonzero bands of the jacobian, (dydot/dy), rotated row-wise.
    
  }
  \item{jactype }{the structure of the Jacobian, one of "fullint", "fullusr",
    "bandusr", "bandint", or "sparse" - either full or banded and
      estimated internally or by the user, or arbitrary sparse.
     If the latter, then the solver will call, \code{\link{stodes}},
     else \code{\link{stode}}

     If the Jacobian is arbitrarily "sparse", then it will be calculated by
     the solver (i.e. it is not possible to also specify \code{jacfunc}).
  }
  \item{verbose }{if \code{TRUE}: full output to the screen, e.g. will output
    the steady-state settings.
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case the
    Jacobian is banded.
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case the
    jacobian is banded.
  }
  \item{... }{additional arguments passed to function \code{f}.
  }
}
\value{
  a list containing:
  \item{root }{the location (x-values) of the root.
  }
  \item{f.root }{the value of the function evaluated at the \code{root}.
  }
  \item{iter }{the number of iterations used.
  }
  \item{estim.precis }{the estimated precision for \code{root}. 
    It is defined as the mean of the absolute function values
    (mean(abs(\code{f.root}))).
  }
}

\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{
## =======================================================================
## example 1  
## 2 simultaneous equations
## =======================================================================

model <- function(x) c(F1=x[1]^2+ x[2]^2 -1,F2=x[1]^2- x[2]^2 +0.5)

(ss<-multiroot(f=model,start=c(1,1)))

## =======================================================================
## example 2
## 3 equations, two solutions
## =======================================================================

model <- function(x) c(F1= x[1] + x[2] + x[3]^2 - 12,
                       F2= x[1]^2 - x[2] + x[3] - 2,
                       F3= 2 * x[1] - x[2]^2 + x[3] - 1 )

# first solution
(ss<-multiroot(model,c(1,1,1),useFortran=FALSE))
(ss<-multiroot(f=model,start=c(1,1,1)))

# second solution; use different start values
(ss<-multiroot(model,c(0,0,0)))
model(ss$root)

## =======================================================================
## example 3: find a matrix
## =======================================================================

f2<-function(x)
 {
 X<-matrix(nr=5,x)
 X \%*\% X \%*\% X -matrix(nr=5,data=1:25,byrow=TRUE)
 }
x<-multiroot(f2, start= 1:25 )$root
X<-matrix(nr=5,x)

X\%*\%X\%*\%X
}
\details{
  \code{start} gives the initial guess for each variable; different initial
  guesses may return different roots.

  The input parameters \code{rtol}, and \code{atol} determine the \bold{error
  control} performed by the solver.  
  
  The solver will control the vector
  \bold{e} of estimated local errors in \bold{f}, according to an
  inequality of the form max-norm of ( \bold{e}/\bold{ewt} )
  \eqn{\leq}{<=} 1, where \bold{ewt} is a vector of positive error
  weights.  The values of \code{rtol} and \code{atol} should all be
  non-negative.

  The form of \bold{ewt} is:
  \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{f}) + \mathbf{atol}}{\bold{rtol} * abs(\bold{x}) + \bold{atol}}
  
  where multiplication of two vectors is element-by-element.

  In addition, the solver will stop if between two iterations, the maximal
  change in the values of \bold{x} is less than \code{ctol}.

  There is no checking whether the requested precision exceeds the capabilities
  of the machine.
}

\seealso{
 \code{\link{stode}}, which uses a different function call.

 \code{\link{uniroot.all}}, to solve for all roots of one (nonlinear) equation

 \code{\link{steady}}, \code{\link{steady.band}},  \code{\link{steady.1D}},
 \code{\link{steady.2D}}, \code{\link{steady.3D}}, steady-state solvers,
   which find the roots of ODEs or PDEs. The function call differs from
   \code{multiroot}.

 \code{\link{jacobian.full}}, \code{\link{jacobian.band}}, estimates the
   Jacobian matrix assuming a full or banded structure.

 \code{\link{gradient}}, \code{\link{hessian}}, estimates the gradient
   matrix or the Hessian.
}

\note{
  The Fortran implementation of the Newton-Raphson method function (the
  default) is generally faster than the \R implementation.
  The \R implementation has been included for didactic purposes.


  \code{multiroot} makes use of function \code{stode}.
  Technically, it is just a wrapper around function \code{stode}.
  If the sparsity structure of the Jacobian is known, it may be more efficiently
  to call \code{stode, stodes, steady, steady.1D, steady.2D, steady.3D}.
  

  It is NOT guaranteed that the method will converge to the root.
}
\keyword{optimize}

