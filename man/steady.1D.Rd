\name{steady.1D}
\alias{steady.1D}
\title{Steady-state solver for multicomponent 1-D ordinary differential equations}
\description{Estimates the steady-state condition for a system of ordinary differential equations that result
from 1-Dimensional reaction-transport models that include transport only between adjacent layers and that model many species.}
\usage{steady.1D(y, time=0, func, parms=NULL, nspec=NULL, dimens=NULL, 
                method = "stode",...)}
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
  \item{nspec }{the number of *species* (components) in the model. If NULL, then \code{dimens} should be specified}
  \item{dimens}{the number of *boxes* in the model. If NULL, then \code{nspec} should be specified}
  \item{method }{the solution method, one of "stode", "stodes" or "runsteady"}
  \item{... }{additional arguments passed to the solver function as defined by \code{method}}
}
\value{
  A list containing 
  \item{y }{ A vector with the state variable values from the last iteration during estimation of steady-state condition of the system of equations.     
    If \code{y} has a names attribute, it will be used to label the output values. }
  \item{... }{ the number of "global" values returned}

  The output will have the attribute \code{steady}, which returns TRUE, if steady-state has been reached and the attribute
    \code{precis} with the precision attained during each iteration.
}
\note{It is advisable though not mandatory to specify BOTH \code{nspec} and \code{dimens}. 
In this case, the solver can check whether the input makes sense 
(i.e. if  nspec*dimens = length(y))}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{##################################################################
######  EXAMPLE 1: BOD + O2                                 ######
##################################################################
# Biochemical Oxygen Demand (BOD) and oxygen (O2) dynamics
# in a river

#==================#
# Model equations  #
#==================#
O2BOD <- function(t,state,pars)

{
  BOD <- state[1:N]
  O2  <- state[(N+1):(2*N)]

# BOD dynamics
  FluxBOD <-  v*c(BOD_0,BOD)  # fluxes due to water transport
  FluxO2  <-  v*c(O2_0,O2)
  
  BODrate <- r*BOD*O2/(O2+10)  # 1-st order consumption, Monod in oxygen

#rate of change = flux gradient - consumption  + reaeration (O2)
  dBOD         <- -diff(FluxBOD)/dx  - BODrate
  dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)

  return(list(c(dBOD=dBOD,dO2=dO2),BODrate=BODrate))

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
print(system.time(
 out   <- steady.1D (y=state,func=O2BOD,parms=NULL, nspec=2,pos=TRUE)))

#==================#
# Plotting output  #
#==================#
mf <- par(mfrow=c(2,2))
plot(x,out$y[(N+1):(2*N)],xlab= "Distance from river", 
     ylab="mmol/m3",main="Oxygen",type="l")

plot(x,out$y[1:N],xlab= "Distance from river",
     ylab="mmol/m3",main="BOD",type="l")

plot(x,out$BODrate,xlab= "Distance from river",
     ylab="mmol/m3/d",main="BOD decay rate",type="l")
par(mfrow=mf)

# same, but now running dynamically to steady-state
print(system.time(
out   <- steady.1D (y=state,func=O2BOD,parms=NULL, nspec=2,
                    time=c(0,1000),method="runsteady")))

##################################################################
######  EXAMPLE 2: Silicate diagenesis                      ######
##################################################################
# Example from the book:
# Soetaert and Herman (2008).
# a practical guide to ecological modelling -
# using R as a simulation platform.
# Springer

#====================#
# Model equations    #
#====================#

SiDIAmodel <- function (time=0,    # time, not used here
                        Conc,      # concentrations: BSi, DSi
                        parms=NULL) # parameter values; not used
{
 BSi<- Conc[1:N]
 DSi<- Conc[(N+1):(2*N)]

# transport           
# diffusive fluxes at upper interface of each layer

# upper concentration imposed (bwDSi), lower: zero gradient
 DSiFlux <- -SedDisp *   IntPor *diff(c(bwDSi ,DSi,DSi[N]))/thick    
 BSiFlux <- -Db      *(1-IntPor)*diff(c(BSi[1],BSi,BSi[N]))/thick 

 BSiFlux[1] <- BSidepo                # upper boundary flux is imposed

# BSi dissolution     
 Dissolution <- rDissSi * BSi*(1.- DSi/EquilSi )^pow 
 Dissolution <- pmax(0,Dissolution)

# Rate of change= Flux gradient, corrected for porosity + dissolution
 dDSi     <- -diff(DSiFlux)/thick/Porosity      +    # transport
              Dissolution * (1-Porosity)/Porosity    # biogeochemistry

 dBSi     <- -diff(BSiFlux)/thick/(1-Porosity)  - Dissolution				

 return(list(c(dBSi=dBSi,dDSi=dDSi), # Rates of changes
        Dissolution=Dissolution,     # Profile of dissolution rates
        DSiSurfFlux =DSiFlux[1],     # DSi sediment-water exchange rate 
        DSIDeepFlux =DSiFlux[N+1],   # DSi deep-water (burial) flux
        BSiDeepFlux =BSiFlux[N+1]))  # BSi deep-water (burial) flux
}

#====================#
# Model run          #
#====================#
# sediment parameters
thick    <- 0.1                       # thickness of sediment layers (cm)
Intdepth <- seq(0,10,by=thick)        # depth at upper interface of layers
Nint     <- length(Intdepth)          # number of interfaces
Depth    <- 0.5*(Intdepth[-Nint] +Intdepth[-1]) # depth at middle of layers
N        <- length(Depth)                       # number of layers

por0    <- 0.9                         # surface porosity (-)
pordeep <- 0.7                         # deep porosity    (-)
porcoef <- 2                           # porosity decay coefficient  (/cm)
# porosity profile, middle of layers
Porosity <- pordeep + (por0-pordeep)*exp(-Depth*porcoef)    
# porosity profile, upper interface 
IntPor   <- pordeep + (por0-pordeep)*exp(-Intdepth*porcoef)  

dB0      <- 1/365           # cm2/day      - bioturbation coefficient
dBcoeff  <- 2
mixdepth <- 5               # cm
Db       <- pmin(dB0,dB0*exp(-(Intdepth-mixdepth)*dBcoeff))

# biogeochemical parameters
SedDisp  <- 0.4             # diffusion coefficient, cm2/d
rDissSi  <- 0.005           # dissolution rate, /day
EquilSi  <- 800             # equilibrium concentration
pow      <- 1
BSidepo  <- 0.2*100          # nmol/cm2/day
bwDSi    <- 150              # mmol/m3

# initial guess of state variables-just random numbers between 0,1
Conc     <- runif(2*N)

# three runs with different deposition rates
BSidepo  <- 0.2*100          # nmol/cm2/day
sol  <- steady.1D (Conc, func=SiDIAmodel, parms=NULL, nspec=2)
CONC <- sol$y

BSidepo  <- 2*100          # nmol/cm2/day
sol2 <- steady.1D (Conc, func=SiDIAmodel, parms=NULL, nspec=2)
CONC <- cbind(CONC,sol2$y)

BSidepo  <- 3*100          # nmol/cm2/day
sol3 <- steady.1D (Conc, func=SiDIAmodel, parms=NULL, nspec=2)
CONC  <- cbind(CONC,sol3$y)

DSi  <- CONC[(N+1):(2*N),]
BSi  <- CONC[1:N,]

#====================#
# plotting output    #
#====================#
par(mfrow=c(2,2))

matplot(DSi,Depth,ylim=c(10,0),xlab="mmolSi/m3 Liquid",
        main="DSi",type="l",lwd=c(1,2,1),col="black")
matplot(BSi,Depth,ylim=c(10,0),xlab="mmolSi/m3 Solid" ,
        main="BSi",type="l",lwd=c(1,2,1),col="black")
legend("right",c("0.2","2","3"),title="mmol/m2/d",
       lwd=c(1,2,1),lty=1:3)
plot(Porosity,Depth,ylim=c(10,0),xlab="-" ,
     main="Porosity",   type="l",lwd=2)
plot(Db,Intdepth,ylim=c(10,0),xlab="cm2/d",
     main="Bioturbation",type="l",lwd=2)
mtext(outer=TRUE,side=3,line=-2,cex=1.5,"SiDIAmodel")}
\details{This is the method of choice for multi-species 1-dimensional models, that are only subjected to transport between adjacent layers 
    \cr More specifically, this method is to be used if the state variables are arranged per species: 
    \cr A[1],A[2],A[3],....B[1],B[2],B[3],.... (for species A, B)) 


Two methods are implemented.
\item The default method rearranges the state variables as A[1],B[1],...A[2],B[2],...A[3],B[3],.... This reformulation leads to a banded Jacobian with (upper and lower) half bandwidth = number of species. 
Then function \code{stode} solves the banded problem.
\item The second method uses function \code{stodes}. Based on the dimension of the problem, the method first calculates the sparsity pattern of the Jacobian, under the assumption 
  that transport is only occurring between adjacent layers. Then \code{stodes} is called to solve the problem.
  \cr As \code{stodes} is used to estimate steady-state, it may be necessary to specify the length of the real work array, \code{lrw}.
  \cr Although a reasonable guess of \code{lrw} is made, it is possible that this will be too low. 
  In this case, \code{steady.1D} will return with an error message telling
  the size of the work array actually needed. In the second try then, set \code{lrw} equal to this number.

  For single-species 1-D models, use \code{\link{steady.band}}. 

 If state variables are arranged as (e.g. A[1],B[1],A[2],B[2],A[3],B[3],... then the model should be solved with  \code{\link{steady.band}}}
\seealso{\code{\link{stode}} and \code{\link{stodes}} for the additional options 
  \item \code{\link{steady}}, for solving steady-state when the jacobian matrix is full 
  \item \code{\link{steady.2D}}, for steady-state estimation of 2-D models 
  \item \code{\link{steady.band}}, for steady-state solution, when the jacobian matrix is banded}
\keyword{math}

