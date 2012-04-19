\name{stode}
\alias{stode}
\title{
  Iterative steady-state solver for ordinary differential equations (ODE)
  and a full or banded Jacobian.
}
\description{
Estimates the steady-state condition for a system of 
  ordinary differential equations (ODE) written in the form:
  \deqn{dy/dt = f(t,y)}   
  
  i.e. finds the values of \code{y} for which f(t,y) = 0.
  
  Uses a newton-raphson method, implemented in Fortran 77.
  
  The system of ODE's is written as an \R function or defined in
  compiled code that has been dynamically loaded.
}
\usage{
stode(y, time = 0, func, parms = NULL, 
      rtol = 1e-6, atol = 1e-8, ctol = 1e-8, 
      jacfunc = NULL, jactype = "fullint", verbose = FALSE, 
      bandup = 1, banddown = 1, positive = FALSE, 
      maxiter = 100, ynames = TRUE, 
      dllname = NULL, initfunc = dllname, initpar = parms, 
      rpar = NULL, ipar = NULL, nout = 0, outnames = NULL, 
      forcings = NULL, initforc = NULL, fcontrol = NULL, ...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ode system, a vector.
    If \code{y} has a name attribute, the names will be used to label the
    output matrix.
  }
  \item{time }{time for which steady-state is wanted; the default is
    \code{time}=0.
  }
  \item{func }{either a user-supplied function that computes the values of the
    derivatives in the ode system (the \emph{model definition}) at time
    \code{time}, or a character string giving the name of a compiled function
    in a dynamically loaded shared library.

    If \code{func}  is a user-supplied function, it must be called as:
    \code{yprime = func(t, y, parms, ...)}.  \code{t} is the time point
    at which the steady-state is wanted, \code{y} is the current estimate of
    the variables in the ode system. If the initial values \code{y} has a names
    attribute, the names will be available inside \code{func}.  \code{parms} is
    a vector of parameters (which may have a names attribute).

    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{time}, and whose next elements (possibly with a
    \code{\link{names}} attribute) are global values that are required as
    output.
    
    The derivatives
    should be specified in the same order as the state variables \code{y}.

    If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{stode()} is called. see Details for more information.
  }
  \item{parms }{other parameters passed to \code{func} and \code{jacfunc}.
  }
  \item{rtol }{relative error tolerance, either a scalar or a vector, one
    value for each \code{y}.
  }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one
    value for each \code{y}.
  }
  \item{ctol }{if between two iterations, the maximal change in \code{y} is
    less than this amount, steady-state is assumed to be reached.
  }
  \item{jacfunc }{if not \code{NULL}, either a user-supplied \R function that
    estimates the Jacobian of the system of differential equations
    dydot(i)/dy(j), or a character string giving the name of a compiled
    function in a dynamically loaded shared library as provided in
    \code{dllname}. In some circumstances, supplying \code{jacfunc}
    can speed up the computations.  The \R calling sequence for
    \code{jacfunc} is identical to that of \code{func}.      

    If the Jacobian is a full matrix, \code{jacfunc} should return a matrix
    dydot/dy, where the ith row contains the derivative of \eqn{dy_i/dt}
    with respect to \eqn{y_j}, or a vector containing the matrix elements
    by columns (the way \R and Fortran store matrices).

    If the Jacobian is banded, \code{jacfunc} should return a matrix containing
    only the nonzero bands of the jacobian, (dydot/dy), rotated row-wise.
  }
  \item{jactype }{the structure of the Jacobian, one of "fullint", "fullusr",
    "bandusr", or "bandint" - either full or banded and estimated internally
    or by the user.
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
  \item{positive }{either a logical or a vector with indices of the state
    variables that have to be non-negative; if \code{TRUE}, all state
    variables \code{y} are forced to be non-negative numbers.
  }
  \item{maxiter }{maximal number of iterations during one call to the solver.\
  }
  \item{ynames }{if FALSE: names of state variables are not passed to function
    \code{func} ; this may speed up the simulation especially for multi-D
    models.
  }
  \item{dllname }{a string giving the name of the shared library (without
    extension) that contains all the compiled function or subroutine
    definitions referred to in \code{func} and \code{jacfunc}.
  }
  \item{initfunc }{if not NULL, the name of the initialisation function
    (which initialises values of parameters), as provided in \file{dllname}.
    See details.
  }
  \item{initpar }{only when \file{dllname} is specified and an initialisation
    function \code{initfunc} is in the dll: the parameters passed to the
    initialiser, to initialise the common blocks (FORTRAN) or global variables
    (C, C++).
  }
  \item{rpar }{only when \file{dllname} is specified: a vector with double
    precision values passed to the dll-functions whose names are specified
    by \code{func} and \code{jacfunc}.
  }
  \item{ipar }{only when \file{dllname} is specified: a vector with integer
    values passed to the dll-functions whose names are specified by \code{func}
    and \code{jacfunc}.
  }
  \item{nout }{only used if \file{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated in
    the compiled function \code{func}, present in the shared library. Note:
    it is not automatically checked whether this is indeed the number of
    output variables calculated in the dll - you have to perform this check
    in the code - see package vignette.
    }
  \item{outnames }{only used if \file{dllname} is specified and
    \code{nout} > 0: the names of output variables calculated in the
    compiled function \code{func}, present in the shared library.
  }
  \item{forcings }{only used if \file{dllname} is specified: a vector with the 
    forcing function values, or a list with the forcing function data sets, 
    each present as a two-columned matrix, with (time,value); interpolation 
    outside the interval [min(\code{times}), max(\code{times})] is done by 
    taking the value at the closest data extreme.

    This feature is here for compatibility with models defined in compiled code
    from package deSolve; see deSolve's package vignette \code{"compiledCode"}.
  }
  \item{initforc }{if not \code{NULL}, the name of the forcing function
    initialisation function, as provided in
    \file{dllname}. It MUST be present if \code{forcings} has been given a
    value.
    See deSolve's package vignette \code{"compiledCode"}.
  }
  \item{fcontrol }{A list of control parameters for the forcing functions.
    See deSolve's package vignette \code{"compiledCode"}.
  }
  \item{... }{additional arguments passed to \code{func} and \code{jacfunc}
    allowing this to be a generic function.
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
  if steady-state has been reached and the attribute \code{precis} with an
  estimate of the precision attained during each iteration, the mean absolute
  rate of change (sum(abs(dy))/n).
}
\author{
  Karline Soetaert <karline.soetaert@nioz.nl>
}
\examples{
## =======================================================================
## Example 1. A simple sediment biogeochemical model
## =======================================================================

model<-function(t, y, pars)
{

with (as.list(c(y, pars)),{

  Min       = r*OM
  oxicmin   = Min*(O2/(O2+ks))
  anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2)

  dOM  = Flux - oxicmin - anoxicmin
  dO2  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)
  dSO4 = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4)
  dHS  = 0.5*anoxicmin   -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)

  list(c(dOM, dO2, dSO4, dHS), SumS = SO4+HS)
})
}

# parameter values
pars <- c(D = 1, Flux = 100, r = 0.1, rox = 1,
          ks = 1, ks2 = 1, BO2 = 100, BSO4 = 10000, BHS = 0)
# initial conditions
y<-c(OM = 1, O2 = 1, SO4 = 1, HS = 1)

# direct iteration  - enforces  positivitiy..
ST <- stode(y = y, func = model, parms = pars, pos = TRUE)

ST

## =======================================================================
## Example 2. 1000 simultaneous equations
## =======================================================================

model <- function (time, OC, parms, decay, ing) {
 # model describing organic Carbon (C) in a sediment, 
 # Upper boundary = imposed flux, lower boundary = zero-gradient
 Flux  <- v * c(OC[1] ,OC) +              # advection
          -Kz*diff(c(OC[1],OC,OC[N]))/dx  # diffusion;
 Flux[1]<- flux     # imposed flux
 
 # Rate of change= Flux gradient and first-order consumption
 dOC   <- -diff(Flux)/dx - decay*OC

 # Fraction of OC in first 5 layers is translocated to mean depth
 dOC[1:5]  <- dOC[1:5] - ing*OC[1:5]
 dOC[N/2]  <- dOC[N/2] + ing*sum(OC[1:5])
 list(dOC)
}

v    <- 0.1    # cm/yr
flux <- 10
dx   <- 0.01
N    <- 1000 
dist <- seq(dx/2,by=dx,len=N)
Kz   <- 1  #bioturbation (diffusion), cm2/yr
print( system.time(
ss   <- stode(runif(N), func = model, parms = NULL, positive = TRUE, 
              decay = 5, ing = 20)))

plot(ss$y[1:N], dist, ylim = rev(range(dist)), type = "l", lwd = 2,
     xlab = "Nonlocal exchange", ylab = "sediment depth",
     main = "stode, full jacobian")

## =======================================================================
## Example 3. Solving a system of linear equations
## =======================================================================

# this example is included to demonstrate how to use the "jactype" option.
# (and that stode is quite efficient).

A <- matrix(nrow = 500, ncol = 500, runif(500*500))
B <- 1:500

# this is how one would solve this in R
print(system.time(X1 <- solve(A, B)))

# to use stode:
# 1. create a function that receives the current estimate of x
# and that returns the difference A\%*\%x-b, as a list:

fun <- function (t, x, p)  # t and p are dummies here...
  list(A\%*\%x-B)

# 2. jfun returns the Jacobian: here this equals "A"
jfun <- function (t, x, p) # all input parameters are dummies
  A

# 3. solve with jactype="fullusr" (a full Jacobian, specified by user)
print (system.time(
  X <- stode(y = 1:500, func = fun, jactype = "fullusr", jacfunc = jfun)
  ))

# the results are the same (within precision)
sum((X1-X$y)^2)
}

\references{
  For a description of the Newton-Raphson method, e.g.

  Press, WH, Teukolsky, SA, Vetterling, WT, Flannery, BP, 1996.
  Numerical Recipes in FORTRAN. The Art of Scientific computing. 2nd edition.
  Cambridge University Press.

  The algorithm uses LINPACK code:

  Dongarra, J.J., J.R. Bunch,  C.B. Moler and G.W. Stewart, 1979.
  LINPACK user's guide, SIAM, Philadelphia.
}

\details{
  The work is done by a Fortran 77 routine that implements the Newton-Raphson
  method. It uses code from LINPACK.

  The form of the \bold{Jacobian} can be specified by \code{jactype} which can
  take the following values:
  \itemize{
    \item jactype = "fullint" : a full jacobian, calculated internally by
      the solver, the default.
    \item jactype = "fullusr" : a full jacobian, specified by user function
      \code{jacfunc}.
    \item jactype = "bandusr" : a banded jacobian, specified by user function
      \code{jacfunc}; the size of the bands specified by \code{bandup} and
      \code{banddown}.
    \item jactype = "bandint" : a banded jacobian, calculated by the solver;
      the size of the bands specified by \code{bandup} and \code{banddown}.
  }
  if \code{jactype}= "fullusr" or "bandusr" then the user must supply a
  subroutine \code{jacfunc}.

  The input parameters \code{rtol}, \code{atol} and \code{ctol} determine
  the \bold{error control} performed by the solver.

  The solver will control the vector
  \bold{e} of estimated local errors in \bold{y}, according to an
  inequality of the form max-norm of ( \bold{e}/\bold{ewt} )
  \eqn{\leq}{<=} 1, where \bold{ewt} is a vector of positive error
  weights.  The values of \code{rtol} and \code{atol} should all be
  non-negative.
  The form of \bold{ewt} is:
  
  \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{y}) + \mathbf{atol}}{\bold{rtol} * abs(\bold{y}) + \bold{atol}}
  
  where multiplication of two vectors is element-by-element.

  In addition, the solver will stop if between two iterations, the maximal
  change in the values of \bold{y} is less than \code{ctol}.
  

  \bold{Models} may be defined in compiled C or Fortran code, as well as in \R.

  If \code{func} or \code{jacfunc} are a string, then they are assumed to be
  compiled code.

  In this case, \code{dllname} must give the name of the shared library
  (without extension) which must be loaded before \code{stode()} is called.
  
  See vignette("rooSolve") for how a model has to be specified in compiled code.
  Also, vignette("compiledCode") from package \pkg{deSolve} contains examples
  of how to do this.

  If \code{func} is a user-supplied \bold{R-function}, it must be called as:
  yprime = func(t, y, parms,...).
  t is the \code{time} at which the steady-state should be estimated,
  \code{y} is the current estimate of the variables in the ode system.
  The return value of \code{func} should be a list, whose first element is
  a vector containing the derivatives of \code{y} with respect to \code{time},
  and whose next elements contains output variables whose values at
  steady-state are also required.

  An example is given below:

  \code{model<-function(t,y,pars)}\cr
  \code{\{}\cr
  \code{with (as.list(c(y,pars)),\{}\cr 
     \code{Min       = r*OM}\cr
     \code{oxicmin   = Min*(O2/(O2+ks))}\cr
     \code{anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2}\cr
     \code{dOM  = Flux - oxicmin - anoxicmin}\cr
     \code{dO2  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)}\cr
     \code{dSO4 = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4)}\cr
     \code{dHS  =  0.5*anoxicmin  -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)}
     
     
     \code{list(c(dOM,dO2,dSO4,dHS),SumS=SO4+HS)}\cr
     \code{\})}\cr
  \code{\}}\cr

  This model can be solved as follows:

\code{pars <- c(D=1,Flux=100,r=0.1,rox =1,}\cr
\code{ks=1,ks2=1,BO2=100,BSO4=10000,BHS = 0)}\cr       

\code{y<-c(OM=1,O2=1,SO4=1,HS=1)}\cr                       
\code{ST <- stode(y=y,func=model,parms=pars,pos=TRUE))}\cr
   

}

\seealso{

 \code{\link{steady}}, for a general interface to most of the steady-state
   solvers

 \code{\link{steady.band}}, to find the steady-state of ODE models with a
   banded Jacobian

 \code{\link{steady.1D}}, \code{\link{steady.2D}},
   \code{\link{steady.3D}} steady-state solvers for 1-D, 2-D and 3-D
   partial differential equations.

 \code{\link{stodes}}, iterative steady-state solver for ODEs with arbitrary
   sparse Jacobian.

 \code{\link{runsteady}}, steady-state solver by dynamically running to
   steady-state

}

\note{
  The implementation of \code{stode} and substantial parts of the help file
  is similar to the implementation of the integration routines (e.g.
  \code{lsode}) from package deSolve.
}

\keyword{math}

