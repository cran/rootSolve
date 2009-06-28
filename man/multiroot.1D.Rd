\name{multiroot.1D}
\alias{multiroot.1D}
\title{
  Solves for n roots of n (nonlinear) equations, created by discretizing
  ordinary differential equations.
}
\description{
  multiroot.1D finds the solution to boundary value problems of ordinary
  differential equations, which are approximated using the method-of-lines
  approach.

  Assumes a banded Jacobian matrix, uses the Newton-Raphson method.
}

\usage{
multiroot.1D(f, start, maxiter=100,
  rtol=1e-6, atol=1e-8, ctol=1e-8, 
  nspec = NULL, dimens = NULL,
  verbose=FALSE,positive=FALSE, names=NULL, ...)
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
  \item{nspec }{the number of *species* (components) in the model.
    If \code{NULL}, then \code{dimens} should be specified.
  }
  \item{dimens}{the number of *boxes* in the model. If NULL, then
    \code{nspec} should be specified.
  }
  \item{verbose }{if \code{TRUE}: full output to the screen, e.g. will output
    the steady-state settings.
  }
  \item{positive}{if \code{TRUE}, the unknowns y are forced to be
    non-negative numbers.
  }
  \item{names }{the names of the components; used to label the output, which
    will be written as a matrix.
  }
  \item{... }{additional arguments passed to function \code{f}.
  }
}
\value{
  a list containing:
  \item{root }{the values of the root.
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
#===============================================================================
# Example 1: simple standard problem
# solve the BVP ODE:
# d2y/dt^2=-3py/(p+t^2)^2
# y(t= -0.1)=-0.1/sqrt(p+0.01)
# y(t=  0.1)= 0.1/sqrt(p+0.01)
# where p = 1e-5
#
# analytical solution y(t) = t/sqrt(p + t^2).
#
#===============================================================================

bvp <- function(y) {
  dy2 <- diff(diff(c(ya,y,yb))/dx)/dx
  return(dy2+3*p*y/(p+x^2)^2)
}

dx<-0.0001
x<-seq(-0.1,0.1,by = dx)

p <- 1e-5
ya<- -0.1/sqrt(p+0.01)
yb<- 0.1/sqrt(p+0.01)

print(system.time(y<-multiroot.1D(start=runif(length(x)),f=bvp,nspec=1)))
plot(x,y$root, ylab="y", main="BVP test problem")

# add analytical solution
curve(x/sqrt(p+x*x),add=TRUE,type="l",col="red")

#===============================================================================
# Example 2: bvp test problem 28
# solve:
# xi*y'' + y*y' - y=0
# with boundary conditions:
# y0=1
# y1=3/2
#===============================================================================

prob28 <-function(y,xi) {
 dy2 = diff(diff(c(ya,y,yb))/dx)/dx             # y''
 dy  = 0.5*(diff(c(ya,y)) +diff(c(y,yb)))/dx    # y'  - centered differences

 xi*dy2 +dy*y-y
}

ya <- 1
yb <- 3/2
dx <- 0.001
x <- seq(0,1,by=dx)
print(system.time(
  Y1 <- multiroot.1D(f=prob28,start=runif(length(x)),nspec=1, xi=0.1)
))
Y2<- multiroot.1D(f=prob28,start=runif(length(x)),nspec=1, xi=0.01)
Y3<- multiroot.1D(f=prob28,start=runif(length(x)),nspec=1, xi=0.001)

plot(x,Y3$root,type="l",col="green",main="bvp test problem 28")
lines(x,Y2$root,col="red")
lines(x,Y1$root,col="blue")

}
\details{
  \code{multiroot.1D} is similar to \code{\link{steady.1D}}, except for the
  function specification which is simpler in \code{multiroot.1D}.
  
  It is to be used to solve (simple) boundary value problems of
  differential equations.

  The following differential equation:
  \deqn{0=f(x,y,\frac{dy}{dx},\frac{d^2y}{dx^2})}{0=f(x,y,y',y'')}

  with boundary conditions
  
  \eqn{y_{x=a}}{y(x=a)} = ya, at the start and \eqn{y_{x=b}}{y(x=b)}=yb
  at the end of the integration interval [a,b] is approximated

  as follows:
  
  1. First the integration interval x is discretized, for instance:

  \code{dx <- 0.01}
  
  \code{x <- seq(a,b,by=dx)}

  where \code{dx} should be small enough to keep numerical errors small.

  2. Then the first- and second-order
  derivatives are differenced on this numerical
  grid. R's \code{diff} function is very efficient in taking numerical
  differences, so it is used to approximate the first-, and second-order
  derivates as follows.
  

  A \emph{first-order derivative y'} can be approximated either as:
  \describe{
    \item{y'=\code{diff(c(ya,y))/dx}}{if only the initial condition ya is prescribed,}
    \item{y'=\code{diff(c(y,yb))/dx}}{if only the final condition, yb is prescribed,}
    \item{y'=\code{0.5*(diff(c(ya,y))/dx+diff(c(y,yb))/dx)}}{if initial, ya,
      and final condition, yb are prescribed.}
  }
  The latter (centered differences) is to be preferred.
  
  A \emph{second-order derivative y''} can be approximated by differencing twice.

  y''=\code{diff(diff(c(ya,y,yb))/dx)/dx}

  3. Finally, function \code{multiroot.1D} is used to locate the root.

  See the examples
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
  \code{multiroot.1D} makes use of function \code{steady.1D}.

  It is NOT guaranteed that the method will converge to the root.
}
\keyword{optimize}

