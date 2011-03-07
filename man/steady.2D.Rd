\name{steady.2D}
\alias{steady.2D}
\title{
  Steady-state solver for 2-Dimensional ordinary differential equations
}
\description{
  Estimates the steady-state condition for a system of ordinary differential
  equations that result from 2-Dimensional partial differential equation
  models that have been converted to ODEs by numerical differencing.

  It is assumed that exchange occurs only between adjacent layers.
}
\usage{
steady.2D(y, time = 0, func, parms = NULL, 
          nspec = NULL, dimens,
          names = NULL, cyclicBnd = NULL, ...)
}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector.
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
  \item{nspec }{the number of *species* (components) in the model.
  }
  \item{dimens }{a 2-valued vector with the dimensionality of the model,
    i.e. the number of *boxes* in x- and y-direction.
  }
  \item{cyclicBnd }{if not \code{NULL} then a number or a 2-valued vector
    with the dimensions where a cyclic boundary is used - \code{1}: x-dimension,
    \code{2}: y-dimension; see details.
  }
  \item{names }{the names of the components; used to label the output, 
    and for plotting.
  }
  \item{... }{additional arguments passed to function \code{stodes}.
  }
}
\value{
  A list containing
  \item{y }{ a vector with the state variable values from the last iteration
    during estimation of steady-state condition of the system of equations.
  }
  \item{... }{ the number of "global" values returned.
  }
  The output will have the attribute \code{steady}, which returns \code{TRUE},
  if steady-state has been reached and the attribute
  \code{precis} with the precision attained during each iteration.
}
\note{
  It is advisable though not mandatory to specify BOTH \code{nspec} and
  \code{dimens}. In this case, the solver can check whether the input makes
  sense (as nspec*dimens[1]*dimens[2] = length(y))

  do NOT use this method for problems that are not 2D.
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{
## =======================================================================
## Diffusion in 2-D; imposed boundary conditions
## =======================================================================
diffusion2D <- function(t, Y, par)   {

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
y  <- matrix(nrow = n, ncol = n, 10.)

ST3 <- steady.2D(y, func = diffusion2D, parms = NULL, pos = TRUE,
                 dimens = c(n, n), lrw = 1000000, 
                 atol = 1e-10, rtol = 1e-10, ctol = 1e-10)

image(ST3, legend = TRUE)

# The hard way of plotting:    
#y <- matrix(nr = n, nc = n, data = ST3$y)
#     filled.contour(y, color.palette = terrain.colors)

## =======================================================================
## Diffusion in 2-D; extra flux on 2 boundaries, cyclic in y
## =======================================================================

diffusion2Db <- function(t, Y, par)  {

   y    <- matrix(nr=nx,nc=ny,data=Y)  # vector to 2-D matrix
   dY   <- -r*y        # consumption

   BNDx   <- rep(1,nx)   # boundary concentration
   BNDy   <- rep(1,ny)   # boundary concentration

   #diffusion in X-direction; boundaries=imposed concentration
   Flux <- -Dx * rbind(y[1,]-BNDy,(y[2:nx,]-y[1:(nx-1),]),BNDy-y[nx,])/dx
   dY   <- dY - (Flux[2:(nx+1),]-Flux[1:nx,])/dx

   #diffusion in Y-direction
   Flux <- -Dy * cbind(y[,1]-BNDx,(y[,2:ny]-y[,1:(ny-1)]),BNDx-y[,ny])/dy
   dY    <- dY - (Flux[,2:(ny+1)]-Flux[,1:ny])/dy

   # extra flux on two sides
   dY[,1] <- dY[,1]+  10
   dY[1,] <- dY[1,]+  10

   # and exchange between sides on y-direction
   dY[,ny] <- dY[,ny]+ (y[,1]-y[,ny])*10

   return(list(as.vector(dY)))
}

# parameters
dy    <- dx <- 1   # grid size
Dy    <- Dx <- 1   # diffusion coeff, X- and Y-direction
r     <- 0.025     # consumption rate

nx  <- 50
ny  <- 100
y  <- matrix(nrow = nx, ncol = ny, 10.)
print(system.time(
  ST3 <- steady.2D(y, func = diffusion2Db, parms = NULL, pos = TRUE,
          dimens = c(nx, ny), verbose = TRUE, lrw = 283800, 
          atol = 1e-10, rtol = 1e-10, ctol = 1e-10, 
          cyclicBnd = 2)       # y-direction: cyclic boundary
))
image(ST3)
#y <- matrix(nr=nx,nc=ny,data=ST3$y)
#    filled.contour(y,color.palette=terrain.colors)


}
\details{
  This is the method of choice for 2-dimensional models, that are only
  subjected to transport between adjacent layers.

  Based on the dimension of the problem, the method first calculates the
  sparsity pattern of the Jacobian, under the assumption 
  that transport is only occurring between adjacent layers.
  Then \code{stodes} is called to find the steady-state.
  
  As \code{stodes} is used, it will probably be necessary to specify the
  length of the real work array, \code{lrw}.

  Although a reasonable guess of \code{lrw} is made, it may occur that this
  will be too low.
  In this case, \code{steady.2D} will return with an error message telling
  the size of the work array actually needed. In the second try then, set
  \code{lrw} equal to this number.

  In some cases, a cyclic boundary condition exists. This is when the first
  boxes in x-or y-direction interact with the last boxes. In this case, there
  will be extra non-zero fringes in the Jacobian which need to be taken
  into account. The occurrence of cyclic boundaries can be
  toggled on by specifying argument \code{cyclicBnd}. For innstance,
  \code{cyclicBnd = 1} indicates that a cyclic boundary is required only for
  the x-direction, whereas \code{cyclicBnd = c(1,2)} imposes a cyclic boundary
  for both x- and y-direction. The default is no cyclic boundaries.

  See \code{stodes} for the additional options.
}
\seealso{

 \code{\link{steady}}, for a general interface to most of the steady-state
   solvers

 \code{\link{steady.band}}, to find the steady-state of ODE models with a
   banded Jacobian

 \code{\link{steady.1D}},
   \code{\link{steady.3D}}, steady-state solvers for 1-Dand 3-D
   partial differential equations.

 \code{\link{stode}}, iterative steady-state solver for ODEs with full
   or banded Jacobian.

 \code{\link{stodes}}, iterative steady-state solver for ODEs with arbitrary
   sparse Jacobian.

 \code{\link{runsteady}}, steady-state solver by dynamically running to
   steady-state

}

\keyword{math}

