\name{plot.steady1D}
\alias{plot.steady1D}
\alias{image.steady2D}
\alias{subset.steady2D}
\alias{image.steady3D}
\alias{summary.rootSolve}
\title{
 Plot and Summary Method for steady1D, steady2D and steady3D Objects
}
\description{
Plot the output of steady-state solver routines.
}
\usage{
\method{plot}{steady1D}(x, ..., which = NULL, grid = NULL, 
              xyswap = FALSE, ask = NULL, 
              obs = NULL, obspar = list(), vertical = FALSE)
\method{image}{steady2D}(x, which = NULL, add.contour = FALSE, 
              grid = NULL, ask = NULL, 
              method = "image", legend = FALSE, ...)
\method{subset}{steady2D}(x, which = NULL, ...)
\method{image}{steady3D}(x, which = NULL, dimselect = NULL, 
              add.contour = FALSE, grid = NULL, ask = NULL, 
              method = "image", legend = FALSE, ...)
\method{summary}{rootSolve}(object, ...)
}

\arguments{
  \item{x }{an object of class \code{steady1D}, or \code{steady2D} as 
    returned by the solvers \code{steady.1D} and \code{steady.2D}, and 
    to be plotted.
    
    For \code{steady1D} objects, it is allowed to pass several objects 
    after \code{x} (unnamed) - see second example.
    
    }

  \item{which }{the name(s) or the index to the variables that should be
    plotted.  Default = all variables. }

  \item{grid }{For 1-D plots of output generated with \code{steady.1D}, 
    a vector of values against which the 1-D steady-state solution 
    has to be plotted. If \code{NULL}, then steady-state solutions 
    are plotted against the index. 
    
    for \code{image} plots of output generated with \code{steady.2D} or 
    \code{steady.3D}: the x- and y-grid, as a \code{list}.}

  \item{ask }{logical; if \code{TRUE}, the user is \emph{ask}ed before
    each plot, if \code{NULL} the user is only asked if more than one
    page of plots is necessary and the current graphics device is set
    interactive, see \code{\link{par}(ask=.)} and
    \code{\link{dev.interactive}}.}
  
  \item{xyswap }{if \code{TRUE}, then x-and y-values are swapped and the 
    y-axis is from top to bottom. Useful for drawing vertical profiles.}
    
  \item{vertical }{if \code{TRUE}, then 1. x-and y-values are swapped, the
    y-axis is from top to bottom, the x-axis is on top, margin 3 and the
    main title gets the value of the x-axis.
    Useful for drawing vertical profiles; see example 2.}

  \item{obs }{a \code{data.frame} or \code{matrix} with "observed data" that
    will be added as \code{points} to the plots. \code{obs} can also be a 
    \code{list} with multiple data.frames and/or matrices containing 
    observed data. 
    
    The first column of \code{obs} should contain the grid-variables as 
    specified in argument \code{grid}. 
    The other columns contain the observed values
    and they should have names that are known in \code{x}.    
    
    If the first column of \code{obs} consists of factors, or characters (strings),
    then it is assumed that the data are presented in long (database) format,
    where the first three columns contain (name, grid, value).  
    
    If \code{obs} is not \code{NULL} and 
    \code{which} is \code{NULL}, then the variables, common to both 
    \code{obs} and \code{x} will be plotted.
    }

  \item{obspar }{additional graphics arguments passed to \code{points}, for 
    plotting the observed data. If \code{obs} is a \code{list} containing multiple
    observed data sets, then the graphics arguments can be  a vector
    or a list (e.g. for \code{xlim}, \code{ylim}), specifying each data set 
    separately.
  }

  \item{dimselect }{a \code{list} or \code{NULL}. The dimension over which the 
    3-D image loops.  If \code{NULL}, will loop over the 3rd (z) dimension. This is
    similar as setting \code{dimselect = list(z = 1:Nz)} where Nz is the number
    of grid cells in the 3rd dimension; setting 
    \code{dimselect = list(z = seq(1, Nz, by =2))} will loop over the 3rd dimension,
    but every 2nd cell; \code{dimselect = list(x = ...)} or \code{dimselect = list(y = ...)}
    will loop over the x respectively y-dimension. See \code{\link{steady.3D}}
    for some examples. 
    }

  \item{add.contour }{if \code{TRUE}, will add contours to the image plot. }

  \item{method }{the name of the plotting function to use, one of "image", 
    "filled.contour", "contour" or "persp".}

  \item{legend }{if \code{TRUE}, a color legend will be drawn next to 
    the "image", or "persp" plot.}

  \item{object }{object of class \code{rootSolve} whose summary has to be 
    calculated.
  } 

  \item{\dots}{additional arguments passed to the methods.
  
    The graphical arguments are passed to
    \code{\link{plot.default}} (for 1D) or \code{image} (for 2D, 3D)
    
    For \code{plot.steady1D}, the dots may contain other objects of class 
    \code{steady1D}, as returned by \code{steady.1D}, and 
    to be plotted on the same graphs as \code{x}  - see second example.
    \code{x} and and these other objects should be compatible, 
    i.e. the column names should be the same.
    
    For \code{plot.steady1D}, the arguments after \ldots must be matched exactly.     
    }
}

\details{
  The number of panels per page is automatically determined up to 3 x 3
  (\code{par(mfrow=c(3, 3))}). This default can be overwritten by
  specifying user-defined settings for \code{mfrow} or \code{mfcol}.
  Set \code{mfrow} equal to \code{NULL} to avoid the plotting function to 
  change user-defined \code{mfrow} or \code{mfcol} settings

  Other graphical parameters can be passed as well. Parameters
  are vectorized, either according to the number of plots 
  (\code{xlab, ylab}, \code{main, sub}, \code{xlim, ylim}, \code{log},
   \code{asp, ann, axes, frame.plot},\code{panel.first,panel.last},
   \code{cex.lab,cex.axis,cex.main}) or 
  according to the number of lines within one plot (other parameters 
  e.g. \code{col}, \code{lty}, \code{lwd} etc.) so it is possible to
  assign specific axis labels to individual plots, resp. different plotting 
  style. Plotting parameter \code{ylim}, or \code{xlim} can also be a list 
  to assign different axis limits to individual plots.
  
  Similarly, the graphical parameters for observed data, as passed by 
  \code{obspar} can be vectorized, according to the number of observed 
  data sets.
  
  For \code{steady3D} objects, 2-D images are generated by looping over
  one of the axies; by default the third axis. See  \code{\link{steady.3D}}. 
}

\seealso{
 \code{\link{steady.1D}}, \code{\link{steady.2D}}, \code{\link{steady.3D}}  
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
O2BOD <- function(t, state, pars) {
  BOD <- state[1:N]
  O2  <- state[(N+1):(2*N)]

# BOD dynamics
  FluxBOD <-  v * c(BOD_0, BOD)  # fluxes due to water transport
  FluxO2  <-  v * c(O2_0, O2)
  
  BODrate <- r*BOD*O2/(O2+10)  # 1-st order consumption, Monod in oxygen

#rate of change = flux gradient - consumption  + reaeration (O2)
  dBOD         <- -diff(FluxBOD)/dx  - BODrate
  dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)

  return(list(c(dBOD = dBOD, dO2 = dO2)))
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
state <- c(rep(200,N), rep(200,N))

# running the model
out   <- steady.1D (y = state, func = O2BOD, parms = NULL,
                    nspec = 2, pos = TRUE, 
                    names = c("BOD", "O2"))

summary(out)

# output
plot(out, grid = x, type = "l", lwd = 2, 
     ylab = c("mmol/m3", "mmol O2/m3"))

# observations
obs <- matrix (ncol = 2, data = c(seq(0, 10000, 2000),
                                c(1400, 900,400,100,10,10)))

colnames(obs) <- c("Distance", "BOD")

# plot with observations
plot(out, grid = x, type = "l", lwd = 2, ylab = "mmol/m3", obs = obs, 
     pch = 16, cex = 1.5)

# similar but data in "long" format
OUT <- data.frame(name = "BOD", obs)
\dontrun{
plot(out, grid = x, type = "l", lwd = 2, ylab = "mmol/m3", obs = OBS, 
     pch = 16, cex = 1.5)
}

## =======================================================================
##  EXAMPLE 2: 1D model, BOD + O2 - second run                               
## =======================================================================
# new runs with different v
v       <- 50       # velocity, m/day

# running the model a second time
out2   <- steady.1D (y = state, func = O2BOD, parms = NULL,
                     nspec = 2, pos = TRUE, names = c("BOD", "O2"))


v       <- 200       # velocity, m/day

# running the model a second time
out3   <- steady.1D (y = state, func = O2BOD, parms = NULL,
                     nspec = 2, pos = TRUE, names = c("BOD", "O2"))

# output of all three scenarios at once
plot(out, out2, out3, type = "l", lwd = 2, 
     ylab = c("mmol/m3", "mmol O2/m3"), grid = x,
     obs = obs, which = c("BOD", "O2"))
  
# output of all three scenarios at once, and using vertical style
plot(out, out2, out3, type = "l", lwd = 2, vertical = TRUE,
     ylab = "Distance [m]",
     main = c("BOD [mmol/m3]", "O2 [mmol O2/m3]"), grid = x,
     obs = obs, which = c("BOD", "O2"))

# change plot pars
plot(out, out2, out3, type = "l", lwd = 2, 
     ylab = c("mmol/m3", "mmol O2/m3"), 
     grid = x, col = c("blue", "green"), log = "y",  
     obs = obs, obspar = list(pch = 16, col = "red", cex = 2))

## =======================================================================
## EXAMPLE 3: Diffusion in 2-D; zero-gradient boundary conditions
## =======================================================================

diffusion2D <- function(t,Y,par)  {
   y    <- matrix(nr=n,nc=n,data=Y)  # vector to 2-D matrix
   dY   <- -r*y        # consumption
   BND   <- rep(1,n)   # boundary concentration 

   #diffusion in X-direction; boundaries=imposed concentration
   Flux <- -Dx * rbind(y[1,]-BND, (y[2:n,]-y[1:(n-1),]), BND-y[n,])/dx
   dY   <- dY - (Flux[2:(n+1),]-Flux[1:n,])/dx

   #diffusion in Y-direction
   Flux <- -Dy * cbind(y[,1]-BND, (y[,2:n]-y[,1:(n-1)]), BND-y[,n])/dy
   dY    <- dY - (Flux[ ,2:(n+1)]-Flux[ ,1:n])/dy
                                              
   return(list(as.vector(dY)))
}

  # parameters
dy    <- dx <- 1   # grid size
Dy    <- Dx <- 1   # diffusion coeff, X- and Y-direction
r     <- 0.025     # consumption rate

n  <- 100
Y  <- matrix(nrow = n, ncol = n, 10.)

ST <- steady.2D(Y, func = diffusion2D, parms = NULL, pos = TRUE,
                dimens = c(n, n), lrw = 1000000,
                atol = 1e-10, rtol = 1e-10, ctol = 1e-10)
grid <- list(x = seq(dx/2, by = dx, length.out = n), 
             y = seq(dy/2, by = dy, length.out = n))
image(ST, grid = grid)
summary(ST)
}
\keyword{ hplot }

