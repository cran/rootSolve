\name{multiroot}
\alias{multiroot}
\title{Solves for n roots of n (nonlinear) equations}
\description{Given a vector of n variables, and a set of n (nonlinear) equations in these variables,

estimates the root of the equations, i.e. the variable values where all function values = 0.

Assumes a full Jacobian matrix, uses the Newton-Raphson method}
\usage{multiroot(f, start, maxiter=100, 
  rtol=1e-6, atol=1e-8, ctol=1e-8, 
  useFortran=TRUE, positive=FALSE, ...)}
\arguments{
  \item{f }{function for which the root is sought; it must return a vector with as many values as the length of \code{start}}
  \item{start }{vector containing initial guesses for the unknown x; if \code{start} has a name attribute, the names will be used to label the output vector.}
  \item{maxiter }{maximal number of iterations allowed}
  \item{rtol }{relative error tolerance, either a scalar or a vector, one value for each element in the unknown x }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one value for each element in x}
  \item{ctol }{a scalar. If between two iterations, the maximal change in the variable values is less than this amount, then it is assumed that the root is found}
  \item{useFortran }{logical, if FALSE, then an \R -implementation of the Newton-Raphson method is used - see details}
  \item{positive}{if TRUE, the unknowns y are forced to be non-negative numbers}
  \item{... }{additional arguments passed to function \code{f}}
}
\value{
  a list containing:
  \item{root }{the location (x-values) of the root}
  \item{f.root }{the value of the function evaluated at the \code{root}}
  \item{iter }{the number of iterations used}
  \item{estim.precis }{the estimated precision for \code{root}. 
    It is defined as the mean of the absolute function values (mean(abs(\code{f.root}))) }
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
# example 1  
# 2 simultaneous equations
model <- function(x) c(F1=x[1]^2+ x[2]^2 -1,F2=x[1]^2- x[2]^2 +0.5)

(ss<-multiroot(f=model,start=c(1,1)))

# example 2
# 3 equations, two solutions
model <- function(x) c(F1= x[1] + x[2] + x[3]^2 - 12,
                       F2= x[1]^2 - x[2] + x[3] - 2,
                       F3= 2 * x[1] - x[2]^2 + x[3] - 1 )

# first solution
(ss<-multiroot(model,c(1,1,1),useFortran=FALSE))
(ss<-multiroot(f=model,start=c(1,1,1)))

# second solution; use different start values
(ss<-multiroot(model,c(0,0,0)))
model(ss$root)

# example 3: find a matrix
f2<-function(x)
 {
 X<-matrix(nr=5,x)
 X \%*\% X \%*\% X -matrix(nr=5,data=1:25,byrow=TRUE)
 }
x<-multiroot(f2, start= 1:25 )$root
X<-matrix(nr=5,x)

X\%*\%X\%*\%X
}
\details{\code{start} gives the initial guess for each variable; different initial guesses may return different roots.

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

In addition, the solver will stop if between two iterations, the maximal change in the values of \bold{x} is less than \code{ctol}.

There is no checking whether the requested precision exceeds the capabilities of the machine.}
\seealso{\code{\link{steady}}, \code{\link{steady.1D}}, and \code{\link{steady.band}}, the root solvers 
for a system of ordinary differential equations. Here variables are state-variables and the function
 estimates the rate of change.}
\note{The Fortran implementation of the Newton-Raphson method function (the default) is generally faster than the \R implementation.

\code{multiroot} makes use of function \code{steady}. Technically, \code{multiroot} is just a wrapper around function \code{steady}.

The \R implementation has been included for didactic purposes.

It is NOT guaranteed that the method will converge to the root.}
\keyword{optimize}

