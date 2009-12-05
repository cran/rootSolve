\name{plot.steady1D}
\alias{plot.steady1D}
\alias{image.steady2D}
\title{
 Plot Method for steady1D and steady2D Objects
}
\description{
Plot the output of steady-state solver routines.
}
\usage{
\method{plot}{steady1D}(x, which = NULL, grid = NULL, 
  xyswap =FALSE, ask = NULL, ...)
\method{image}{steady2D}(x, which = NULL, add.contour = FALSE, 
  grid = NULL, ask = NULL, ...)
}

\arguments{
  \item{x }{an object of class \code{steady1D}, or \code{steady2D} as 
    returned by the solvers \code{steady.1D} and \code{steady.2D}, and 
    to be plotted.}

  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables. }

  \item{grid }{For 1-D plots of output generated with \code{steady.1D}, 
    a vector of values against which the 1-D steady-state solution 
    has to be plotted. If \code{NULL}, then steady-state solutions 
    are plotted against the index. 
    
    for \code{image} plots of output generated with \code{steady.2D}: 
    the x- and y-grid, as a \code{list}.}

  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.}
  
  \item{xyswap }{if \code{TRUE}, then x-and y-values are swapped and the 
    y-axis is from top to bottom. Useful for drawing vertical profiles.}
    
  \item{add.contour }{if \code{TRUE}, will add contours to the image plot. }

  \item{\dots}{additional graphics arguments passed to
    \code{\link{plot.default}} (for 1D) or \code{image} (for 2D)}
}

\details{
  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow=c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.

  Other graphical parameters can be passed as well. Parameters
  \code{xlab} and \code{ylab} are vectorized, so it is possible to
  assign specific axis labels to individual plots.
}

\seealso{
 \code{\link{steady.1D}}, \code{\link{steady.2D}}  
}
\examples{
## =======================================================================
##  EXAMPLE 1: 1D model, BOD + O2                                
## =======================================================================
## Biochemical Oxygen Demand (BOD) and oxygen (O2) dynamics
## in a river

#==================#
# Model equations  #
#==================#
O2BOD <- function(t,state,pars) {
  BOD <- state[1:N]
  O2  <- state[(N+1):(2*N)]

# BOD dynamics
  FluxBOD <-  v*c(BOD_0,BOD)  # fluxes due to water transport
  FluxO2  <-  v*c(O2_0,O2)
  
  BODrate <- r*BOD*O2/(O2+10)  # 1-st order consumption, Monod in oxygen

#rate of change = flux gradient - consumption  + reaeration (O2)
  dBOD         <- -diff(FluxBOD)/dx  - BODrate
  dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)

  return(list(c(dBOD=dBOD,dO2=dO2)))

 }    # END O2BOD
 
 
#==================#
# Model application#
#==================#
# parameters
dx      <- 100       # grid size, meters
v       <- 1e2       # velocity, m/day
x       <- seq(dx/2,10000,by=dx)  # m, distance from river
N       <- length(x)
r       <- 0.1       # /day, first-order decay of BOD
p       <- 0.1       # /day, air-sea exchange rate
O2sat   <- 300       # mmol/m3 saturated oxygen conc
O2_0    <- 50        # mmol/m3 riverine oxygen conc
BOD_0   <- 1500      # mmol/m3 riverine BOD concentration

# initial guess:
state <- c(rep(200,N),rep(200,N))

# running the model
 out   <- steady.1D (y=state, func=O2BOD, parms=NULL,
                     nspec=2, pos=TRUE, names=c("BOD","O2"))

# output
 plot(out, type="l", lwd=2, ylab=c("mmol/m3","mmol O2/m3"), grid=x)

## =======================================================================
## Diffusion in 2-D; zero-gradient boundary conditions
## =======================================================================

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
  Y  <- matrix(nr=n,nc=n,10.)

  ST <- steady.2D(Y,func=diffusion2D,parms=NULL,pos=TRUE,dimens=c(n,n),
                 lrw=1000000,atol=1e-10,rtol=1e-10,ctol=1e-10)
  grid <- list(x=seq(dx/2, by = dx, length.out=n), 
               y=seq(dy/2, by = dy, length.out=n))
  image(ST, grid = grid)

}
\keyword{ hplot }

