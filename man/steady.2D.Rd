\name{steady.2D}
\alias{steady.2D}
\title{Steady-state solver for 2-Dimensional ordinary differential equations}
\description{Estimates the steady-state condition for a system of ordinary differential equations that
result from 2-Dimensional reaction-transport models that include transport only between adjacent layers }
\usage{steady.2D(y, time=0, func, parms=NULL, nspec=NULL, dimens,...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{time }{time for which steady-state is wanted; the default is time=0}
  \item{func }{either an \R-function that computes the values of the derivatives in the ode system (the model defininition) at time \code{time}, 
    or a character string giving the name of a compiled function in a dynamically loaded shared library.
          If \code{func}  is an \R-function, it must be defined as:
      \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
      in the integration, \code{y} is the current estimate of the variables
      in the ODE system.  If the initial values \code{y} has a names
      attribute, the names will be available inside \code{func}.  \code{parms} is
      a vector or list of parameters; ... (optional) are any other arguments passed to the function. 
  
      The return value of \code{func} should be a list, whose first element is a
      vector containing the derivatives of \code{y} with respect to
      \code{time}, and whose next elements are global values whose steady-state value is also required.}
  \item{parms }{parameters passed to \code{func}}
  \item{nspec }{the number of *species* (components) in the model. }
  \item{dimens}{a 2-valued vector with the dimensionality of the model, i.e. the number of *boxes* in x- and y-direction}
  \item{... }{additional arguments passed to function \code{stodes}}
}
\value{
   A list containing 
  \item{y }{ A vector with the state variable values from the last iteration during estimation of steady-state condition of the system of equations.     
    If \code{y} has a names attribute, it will be used to label the output values. }
  \item{... }{ the number of "global" values returned}
  The output will have the attribute \code{steady}, which returns TRUE, if steady-state has been reached and the attribute
    \code{precis} with the precision attained during each iteration.
}
\note{It is advisable though not mandatory to specify BOTH \code{nspec} and \code{dimens}. In this case, the solver can check whether the input makes sense 
(as nspec*dimens[1]*dimens[2] = length(y))

do NOT use this method for problems that are not 2D}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
#############################################################
# Diffusion in 2-D; imposed boundary conditions
#############################################################
diffusion2D <- function(t,Y,par)
  {
   y    <- matrix(nr=n,nc=n,data=Y)  # vector to 2-D matrix
   dY   <- -r*y        # consumption
   BND   <- rep(1,n)   # boundary concentration 

   #diffusion in X-direction; boundaries=imposed concentration
   Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
   dY   <- dY - (Flux[2:(n+1),]-Flux[1:n,])/dx

   #diffusion in Y-direction
   Flux <- -Dy * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dy
   dY    <- dY - (Flux[,2:(n+1)]-Flux[,1:n])/dy

   return(list(as.vector(dY)))
  }

  # parameters
  dy    <- dx <- 1   # grid size
  Dy    <- Dx <- 1   # diffusion coeff, X- and Y-direction
  r     <- 0.025     # consumption rate

  n  <- 100
  y  <- matrix(nr=n,nc=n,10.)

  ST3 <- steady.2D(y,func=diffusion2D,parms=NULL,pos=TRUE,dimens=c(n,n),
                 lrw=1000000,atol=1e-10,rtol=1e-10,ctol=1e-10)
  y <- matrix(nr=n,nc=n,data=ST3$y)
    filled.contour(y,color.palette=terrain.colors)
    
}
\details{This is the method of choice for 2-dimensional models, that are only subjected to transport between adjacent layers. 

  Based on the dimension of the problem, the method first calculates the sparsity pattern of the Jacobian, under the assumption 
  that transport is onely occurring between adjacent layers. Then \code{stodes} is called to find the steady-state.
  
  As \code{stodes} is used, it will probably be necessary to specify the length of the real work array, \code{lrw}.

Although a reasonable guess of \code{lrw} is made, it is likely that this will be too low. 
  In this case, \code{steady.2D} will return with an error message telling
  the size of the work array actually needed. In the second try then, set \code{lrw} equal to this number.

  See \code{stodes} for the additional options}
\seealso{\code{\link{stodes}} for the additional options 
  \item \code{\link{steady}}, for solving steady-state when the jacobian matrix is full 
  \item \code{\link{steady.1D}}, for solving steady-state for 1-D models 
  \item \code{\link{steady.2D}}, for steady-state estimation of 2-D models 
  \item \code{\link{steady.band}}, for steady-state solution, when the jacobian matrix is banded  }
\keyword{math}

