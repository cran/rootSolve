\name{steady.band}
\alias{steady.band}
\title{
  Steady-state solver for ordinary differential equations; assumes a
  banded jacobian
}
\description{
  Estimates the steady-state condition for a system of ordinary differential
  equations.

  Assumes a banded jacobian matrix.

}
\usage{
steady.band(y, time=0, func, parms=NULL, 
  nspec=NULL, bandup=nspec, banddown=nspec, ...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ODE system, a vector.
    If \code{y} has a name attribute, the names will be used to label the
    output matrix.
  }
  \item{time }{time for which steady-state is wanted; the default is
    \code{time}=0.
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
  \item{bandup }{the number of nonzero bands above the Jacobian diagonal.
  }
  \item{banddown }{the number of nonzero bands below the Jacobian diagonal.
  }
  \item{... }{additional arguments passed to function \code{stode}.
  }
}
\value{
  A list containing 
  \item{y }{ a vector with the state variable values from the last iteration
    during estimation of steady-state condition of the system of equations.
    If \code{y} has a names attribute, it will be used to label the output
    values.
  }
  \item{... }{ the number of "global" values returned.
  }
  The output will have the attribute \code{steady}, which returns \code{TRUE},
  if steady-state has been reached and the attribute
  \code{precis} with the precision attained during each iteration.
}
\author{
  Karline Soetaert <k.soetaert@nioo.knaw.nl>
}
\examples{
#############################################################
# 1000 simultaneous equations, solved 6 times for different 
# values of parameter "decay"
#############################################################
  model <- function (time,OC,parms,decay)
  {
   # model of particles (OC) that sink out of the water and decay
   # exponentially declining sinking rate, maximal 100 m/day
   sink <- 100*exp(-0.01*dist)
  
   # boundary flux; upper boundary=imposed concentration (100)
   Flux <- sink * c(100 ,OC)     
   
   # Rate of change= Flux gradient and first-order consumption
   dOC  <- -diff(Flux)/dx - decay*OC
   list(dOC,maxC=max(OC))
  }
  dx   <- 1                          # thickness of boxes
  dist <- seq(0,1000,by=dx)          # water depth at each modeled box interface

  ss   <- NULL 
  for (decay in seq(from=0.1,to=1.1,by=0.2))
    ss   <- cbind(ss,steady.band(runif(1000),func=model,
    parms=NULL,nspec=1,decay=decay)$y)  

  matplot(ss,1:1000,type="l",lwd=2,main="steady.band", ylim=c(1000,0),
          ylab="water depth, m", xlab="concentration of sinking particles")
  legend("bottomright",legend=seq(from=0.1,to=1.1,by=0.2),lty=1:10,
          title="decay rate",col=1:10,lwd=2)

#############################################################
# 5001 simultaneous equations: solve
# dy/dt = 0 = d2y/dx2 + 1/x*dy/dx + (1-1/(4x^2)y - sqrx(x)*cos(x),
# over the interval [1,6]
# with boundary conditions: y(1)=1, y(6)=-0.5
#############################################################

  derivs <- function(t,y,parms, x,dx,N,y1,y6)
  {

  # Numerical approximation of derivates:
  # d2y/dx2 = (yi+1-2yi+yi-1)/dx^2
   d2y <- (c(y[-1],y6) -2*y + c(y1,y[-N])) /dx/dx

  # dy/dx = (yi+1-yi-1)/(2dx)
   dy  <- (c(y[-1],y6) - c(y1,y[-N])) /2/dx

   res <- d2y+dy/x+(1-1/(4*x*x))*y-sqrt(x)*cos(x)
   return(list(res))
  }

dx     <- 0.001
x      <- seq(1,6,by=dx)
N      <- length(x)
y  <- steady.band(y=rep(1,N),time=0,func=derivs,x=x,dx=dx,
                  N=N,y1=1,y6=-0.5,nspec=1)$y
plot(x,y,type="l",main="5001 nonlinear equations - banded Jacobian")

# add the analytic solution for comparison:
xx     <- seq(1,6,by=0.1)
ana <- 0.0588713*cos(xx)/sqrt(xx)+1/4*sqrt(xx)*cos(xx)+
       0.740071*sin(xx)/sqrt(xx)+1/4*xx^(3/2)*sin(xx)
points(xx,ana)
legend("topright",pch=c(NA,1),lty=c(1,NA),c("numeric","analytic"))

          }
\details{
  This is the method of choice for single-species 1-D models.

  For multi-species 1-D models, this method can only be used if the state
  variables are arranged per box, per species (e.g.
  A[1],B[1],A[2],B[2],A[3],B[3],.... for species A, B).

  Usually a 1-D *model* function will have the species arranged as
  A[1],A[2],A[3],....B[1],B[2],B[3],....
  in this case, use \code{steady.1D}
}
\seealso{

 \code{\link{steady}}, for a general interface to most of the steady-state
   solvers

 \code{\link{steady.1D}}, \code{\link{steady.2D}},
   \code{\link{steady.3D}}, steady-state solvers for 1-D, 2-D and 3-D
   partial differential equations.

 \code{\link{stode}}, iterative steady-state solver for ODEs with full
   or banded Jacobian.

 \code{\link{stodes}}, iterative steady-state solver for ODEs with arbitrary
   sparse Jacobian.

 \code{\link{runsteady}}, steady-state solver by dynamically running to
   steady-state

}

\keyword{math}

