\name{steady.3D}
\alias{steady.3D}
\title{
  Steady-state solver for 3-Dimensional ordinary differential equations
}
\description{
  Estimates the steady-state condition for a system of ordinary differential
  equations that result from 3-Dimensional partial differential equation
  models that have been converted to ODEs by numerical differencing.

  It is assumed that exchange occurs only between adjacent layers.
}
\usage{
steady.3D(y, time = 0, func, parms = NULL, nspec = NULL, 
          dimens, names = NULL, method = "stodes", 
          jactype = NULL, cyclicBnd = NULL, times = time, ...)
}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector.
  }
  \item{time, times }{time for which steady-state is wanted; the default is \code{times=0}.
    (note- since version 1.7, 'times' has been added as an alias to 'time').
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
  \item{dimens }{a 3-valued vector with the dimensionality of the model,
    i.e. the number of *boxes* in x-, y- and z- direction.
  }
  \item{method }{the solution method, one of "stodes", or "runsteady".
    When \code{method} = 'runsteady', then solver \code{lsodes}, the sparse solver is used by default, 
    unless argument \code{jactype} is set to \code{"2D"}, in which case \code{lsode} will be used (likely less efficient).
    in which case \code{lsodes} is used and the structure of the jacobian is determined by the solver.
  }
  \item{jactype }{the jacobian type - default is a regular 2-D structure where layers only interact with adjacent layers in both directions. 
    If the structure does not comply with this, the jacobian can be set equal to \code{'sparse'}.
  }
  \item{cyclicBnd }{if not \code{NULL} then a number or a 3-valued vector
    with the dimensions where a cyclic boundary is used - \code{1}: x-dimension,
    \code{2}: y-dimension; \code{3}: z-dimension;see details.
  }
  \item{names }{the names of the components; used to label the output, 
    and for plotting.
  }
  \item{... }{additional arguments passed to function \link{stodes}.
    See note.
  }
}
\value{
  A list containing
  \item{y }{ a vector with the state variable values from the last iteration
    during estimation of steady-state condition of the system of equations.
  }
  \item{... }{ the "global" values returned.
  }
  The output will have the attribute \code{steady}, which returns \code{TRUE},
  if steady-state has been reached and the attribute

  \code{precis} with the precision attained during each iteration.
  Another attribute, called \code{dims} returns a.o. the length of the work 
  array actually required. 
  This can be specified with input argument \code{lrw}. See note and first example.
}
\note{
  It is advisable though not mandatory to specify BOTH \code{nspec} and
  \code{dimens}. In this case, the solver can check whether the input makes
  sense (as nspec*dimens[1]*dimens[2]*dimens[3] = length(y))

  do NOT use this method for problems that are not 3D.

  It is likely that the estimated length of the work array (argument \code{lrw}),
  required for the solver \link{stodes} will be too small. 
  If that is the case, the solver will return with an error
  saying to increase \code{lrw}. The current value of the work array can be 
  found via the \code{attributes} of the output. See first example.
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
## =======================================================================
## Diffusion in 3-D; imposed boundary conditions
## =======================================================================
diffusion3D <- function(t, Y, par)   {
   yy    <- array(dim=c(n,n,n),data=Y)  # vector to 3-D array
   dY   <- -r*yy        # consumption
   BND   <- rep(1,n)   # boundary concentration
   for (i in 1:n) {
     y <- yy[i,,]
     #diffusion in X-direction; boundaries=imposed concentration
     Flux <- -Dy * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dy
     dY[i,,]   <- dY[i,,] - (Flux[2:(n+1),]-Flux[1:n,])/dy

     #diffusion in Y-direction
     Flux <- -Dz * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dz
     dY[i,,]    <- dY[i,,] - (Flux[,2:(n+1)]-Flux[,1:n])/dz
   }
   for (j in 1:n) {
     y <- yy[,j,]
     #diffusion in X-direction; boundaries=imposed concentration
     Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
     dY[,j,]   <- dY[,j,] - (Flux[2:(n+1),]-Flux[1:n,])/dx

     #diffusion in Y-direction
     Flux <- -Dz * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dz
     dY[,j,]    <- dY[,j,] - (Flux[,2:(n+1)]-Flux[,1:n])/dz
   }
   for (k in 1:n) {
     y <- yy[,,k]
     #diffusion in X-direction; boundaries=imposed concentration
     Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
     dY[,,k]   <- dY[,,k] - (Flux[2:(n+1),]-Flux[1:n,])/dx

     #diffusion in Y-direction
     Flux <- -Dy * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dy
     dY[,,k]    <- dY[,,k] - (Flux[,2:(n+1)]-Flux[,1:n])/dy
   }
   return(list(as.vector(dY)))
}

  # parameters
  dy    <- dx <- dz <-1   # grid size
  Dy    <- Dx <- Dz <-1   # diffusion coeff, X- and Y-direction
  r     <- 0.025     # consumption rate

  n  <- 10
  y  <- array(dim=c(n, n, n), data = 10.)

  # stodes is used, so we should specify the size of the work array (lrw)
  # We take a rather large value initially

  print(system.time(
  ST3 <- steady.3D(y, func =diffusion3D, parms = NULL, pos = TRUE,
                   dimens = c(n, n, n), lrw = 100000,
                   atol = 1e-10, rtol = 1e-10, ctol = 1e-10, 
                   verbose = TRUE)
  ))

  # the actual size of lrw is in the attributes()$dims vector.     
  # it is best to set lrw as small as possible 
  attributes(ST3)     

  # image plot
  y <- array(dim=c(n, n, n), data = ST3$y)
  filled.contour(y[ , ,n/2], color.palette = terrain.colors)
    
  # rootSolve's image plot, looping over 3rd dimension
  image(ST3, mfrow = c(4,3))

  # loop over 1st dimension, contours, legends added
  image(ST3, mfrow = c(2, 2), add.contour = TRUE, legend = TRUE,
        dimselect = list(x = c(1, 4, 8, 10)))

}
\details{
  This is the method of choice to find the steady-state for 
  3-dimensional models, that are only subjected to transport between 
  adjacent layers.

  Based on the dimension of the problem, the method first calculates the
  sparsity pattern of the Jacobian, under the assumption 
  that transport is only occurring between adjacent layers.
  Then \code{stodes} is called to find the steady-state.
  
  In some cases, a cyclic boundary condition exists. This is when the first
  boxes in x-, y-, or z-direction interact with the last boxes.
  In this case, there will be extra non-zero fringes in the Jacobian which
  need to be taken into account. The occurrence of cyclic boundaries can be
  toggled on by specifying argument \code{cyclicBnd}. For innstance,
  \code{cyclicBnd = 1} indicates that a cyclic boundary is required only for
  the x-direction, whereas \code{cyclicBnd = c(1,2)} imposes a cyclic boundary
  for both x- and y-direction. The default is no cyclic boundaries.

  As \code{stodes} is used, it will probably be necessary to specify the
  length of the real work array, \code{lrw}.

  Although a reasonable guess of \code{lrw} is made, it may occur that this
  will be too low.
  In this case, \code{steady.3D} will return with an error message telling
  the size of the work array actually needed. In the second try then, set
  \code{lrw} equal to this number.

  As \code{stodes} is used, it will probably be necessary to specify the
  length of the real work array, \code{lrw}.

  Although a reasonable guess of \code{lrw} is made, it may occur that this
  will be too low.
  In this case, \code{steady.2D} will return with an error message telling
  that there was insufficient storage. In the second try then, increase
  \code{lrw}. you may need to experiment to find suitable value. The smaller the better.

  The error message that says to increase \code{lrw} may look like this:
\preformatted{   
 In stodes(y = y, time = time, func = func, parms = parms, nnz = c(nspec,  :
  insufficient storage in nnfc
}      

  See \code{stodes} for the additional options.
}
\seealso{

 \code{\link{steady}}, for a general interface to most of the steady-state
   solvers

 \code{\link{steady.band}}, to find the steady-state of ODE models with a
   banded Jacobian

 \code{\link{steady.1D}}, \code{\link{steady.2D}},
   steady-state solvers for 1-D and 2-D
   partial differential equations.

 \code{\link{stode}}, iterative steady-state solver for ODEs with full
   or banded Jacobian.

 \code{\link{stodes}}, iterative steady-state solver for ODEs with arbitrary
   sparse Jacobian.

 \code{\link{runsteady}}, steady-state solver by dynamically running to
   steady-state

}

\keyword{math}

