\name{runsteady}
\alias{runsteady}
\title{
  Dynamically runs a system of ordinary differential equations (ODE) to
  steady-state
}
\description{
  Solves the steady-state condition of ordinary differential equations (ODE)
  in the form:
    \deqn{dy/dt = f(t,y)}  
  
  by dynamically running till the summed absolute values of the derivatives
  become smaller than some predefined tolerance.
    
  The \R function \code{runsteady} makes use of the FORTRAN ODE solver DLSODE,
  written by Alan C. Hindmarsh and Andrew H. Sherman
    
  The system of ODE's is written as an \R function or defined in
  compiled code that has been dynamically loaded. The user has to specify
  whether or not the problem is stiff and choose the appropriate solution
  method (e.g. make choices about the type of the Jacobian).
}

\usage{
runsteady(y, times = c(0, Inf), func, parms, 
          stol = 1e-8, rtol = 1e-6, atol = 1e-6,  
          jacfunc = NULL, jactype = "fullint", mf = NULL, 
          verbose = FALSE, tcrit = NULL, hmin = 0, hmax = NULL, 
          hini = 0, ynames = TRUE, maxord = NULL, bandup = NULL, 
          banddown = NULL, maxsteps = 100000, dllname = NULL, 
          initfunc = dllname, initpar = parms, rpar = NULL, 
          ipar = NULL, nout = 0, outnames = NULL, 
          forcings = NULL, initforc = NULL, fcontrol = NULL, ...)
}

\arguments{
  \item{y }{the initial (state) values for the ODE system. If \code{y} has a
    name attribute, the names will be used to label the output matrix.
  }
  \item{times }{The simulation time. This should be a 2-valued vector,
    consisting of the initial time and the end time.
    The last time value should be large enough to make sure that steady-state
    is effectively reached in this period.
    The simulation will stop either when \code{times[2]} has been reached or
    when \code{maxsteps} have been performed.
  }
  \item{func }{either an \R-function that computes the values of the
    derivatives in the ODE system (the \emph{model definition}) at time
    t, or a character string giving the name of a compiled function
    in a dynamically loaded shared library.
  
    If \code{func}  is an \R-function, it must be defined as:
    \code{yprime = func(t, y, parms,...)}.  \code{t} is the current time point
    in the integration, \code{y} is the current estimate of the variables
    in the ODE system.  If the initial values \code{y} has a names
    attribute, the names will be available inside \code{func}.  \code{parms} is
    a vector or list of parameters; ... (optional) are any other arguments
    passed to the function.
  
    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{time}, and whose next elements are global values that are required at
    each point in \code{times}.
    
    The derivatives
    should be specified in the same order as the state variables \code{y}.
  }
  \item{parms }{vector or list of parameters used in \code{func} or
    \code{jacfunc}.
  }
  \item{stol }{steady-state tolerance; it is assumed that steady-state is
    reached if the average of absolute values of the derivatives
    drops below this number.
  }
  \item{rtol }{relative error tolerance of integrator, either a scalar or an
    array as long as \code{y}.  See details.
  }
  \item{atol }{absolute error tolerance of integrator, either a scalar or an
    array as long as \code{y}.  See details.
  }
  \item{jacfunc }{if not \code{NULL}, an \R function that computes
    the jacobian of the system of differential equations
    dydot(i)/dy(j), or a string giving the name of a function or
    subroutine in \file{dllname} that computes the jacobian (see Details
    below for more about this option).  In some circumstances, supplying
    \code{jacfunc} can speed up the computations, if the system is stiff.
    The \R calling sequence for \code{jacfunc} is identical
    to that of \code{func}.

    If the jacobian is a full matrix, \code{jacfunc} should return a matrix
    dydot/dy, where the ith row contains the derivative of \eqn{dy_i/dt} with
    respect to \eqn{y_j}, or a vector containing the matrix elements by
    columns (the way \R and Fortran store matrices).

    If the jacobian is banded, \code{jacfunc} should return a matrix containing
    only the nonzero bands of the jacobian, rotated row-wise. See first example
    of \code{lsode}.
  }
  \item{jactype }{the structure of the jacobian, 
    one of "fullint", "fullusr", "bandusr" or "bandint" - either full or
    banded and estimated internally or by user; overruled if \code{mf}
    is not NULL.
  }
  \item{mf }{the "method flag" passed to function lsode - overrules
    \code{jactype} - provides more options than \code{jactype} - see details.
  }
  \item{verbose }{if \code{TRUE}: full output to the screen, e.g. will output
    the settings of vectors *istate* and *rstate* - see details.
  }
  \item{tcrit }{if not \code{NULL}, then \code{lsode} cannot integrate past
    \code{tcrit}. The Fortran routine \code{lsode} overshoots its targets
    (times points in the vector \code{times}), and interpolates values
    for the desired time points.  If there is a time beyond which integration
    should not proceed (perhaps because of a singularity),
    that should be provided in \code{tcrit}.
    }
  \item{hmin }{an optional minimum value of the integration
    stepsize. In special situations this parameter may speed up computations
    with the cost of precision. Don't use \code{hmin} if you don't know why!
  }
  \item{hmax }{an optional maximum value of the integration stepsize. If not
    specified, \code{hmax} is set to the largest difference in \code{times},
    to avoid that the simulation possibly ignores short-term events. If 0,
    no maximal size is specified.
  }
  \item{hini }{initial step size to be attempted; if 0, the initial step size
    is determined by the solver.
  }
  \item{ynames }{if \code{FALSE}: names of state variables are not passed to
    function \code{func} ; this may speed up the simulation.
  }
  \item{maxord }{the maximum order to be allowed. \code{NULL} uses the default,
    i.e. order 12 if implicit Adams method (meth=1), order 5 if BDF method
    (meth=2). Reduce maxord to save storage space.
  }
  \item{bandup }{number of non-zero bands above the diagonal, in case the
    jacobian is banded.
  }
  \item{banddown }{number of non-zero bands below the diagonal, in case the
    jacobian is banded.
  }
  \item{maxsteps }{maximal number of steps. The simulation will stop either
    when \code{maxsteps} have been performed or when \code{times[2]} has been
    reached.
  }
  \item{dllname }{a string giving the name of the shared library (without
    extension) that contains all the compiled function or subroutine
    definitions refered to in \code{func} and \code{jacfunc}. See package
    vignette.
  }
  \item{initfunc }{if not NULL, the name of the initialisation function (which
    initialises values of parameters), as provided in \file{dllname}. See
    package vignette.
  }
  \item{initpar }{only when \file{dllname} is specified and an initialisation
    function \code{initfunc} is in the dll: the parameters passed to the
    initialiser, to initialise the common blocks (fortran) or global variables
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
  \item{nout }{only used if \code{dllname} is specified and the model is
    defined in compiled code: the number of output variables calculated in
    the compiled function \code{func}, present in the shared library. Note:
    it is not automatically checked whether this is indeed the number of output
    variables calculed in the dll - you have to perform this check in the code
    - See package vignette of deSolve.
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
  if steady-state has been reached, the attribute \code{precis} with the
  precision attained at the last iteration estimated as the mean absolute
  rate of change (sum(abs(dy))/n), the attribute \code{time} with the
  simulation time reached and the attribute \code{steps} with the number of
  steps performed.
      
  The output will also have the attributes \code{istate}, and \code{rstate},
  two vectors with several useful elements of the dynamic simulation.
  See details.
  The first element of istate returns the conditions under which the last call
  to the integrator returned. Normal is \code{istate[1] = 2}.
  If \code{verbose} = \code{TRUE}, the settings of istate and rstate will
  be written to the screen.

}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{
## =======================================================================
## A simple sediment biogeochemical model
## =======================================================================

model<-function(t, y, pars) {

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
y <- c(OM = 1, O2 = 1, SO4 = 1, HS = 1)

# direct iteration
print( system.time(
  ST <- stode(y = y, func = model, parms = pars, pos = TRUE)
))

print( system.time(
  ST2 <- runsteady(y = y, func = model, parms = pars, times = c(0, 1000))
))

rbind("Newton Raphson" = ST$y, "Runsteady" = ST2$y)
}
\references{ Alan C. Hindmarsh, "ODEPACK, A Systematized Collection of ODE
   Solvers," in Scientific Computing, R. S. Stepleman, et al., Eds.
   (North-Holland, Amsterdam, 1983), pp. 55-64.}
\details{The work is done by the Fortran subroutine \code{dlsode},
    whose documentation should be consulted for details (it is included as
    comments in the source file \file{src/lsode.f}).  The implementation is
    based on the November, 2003 version of lsode, from Netlib.
    
  Before using \code{runsteady}, the user has to decide whether or not the
  problem is stiff. 
  
  If the problem is nonstiff, use method flag \code{mf} = 10, which selects
  a nonstiff (Adams) method, no Jacobian used..  

  If the problem is stiff, there are four standard choices
  which can be specified with \code{jactype} or \code{mf}.
  
  The options for \bold{jactype} are
  \itemize{
    \item jactype = "fullint" : a full jacobian, calculated internally by
      \code{lsode}, corresponds to \code{mf}=22.
    \item jactype = "fullusr" : a full jacobian, specified by user function
      \code{jacfunc}, corresponds to \code{mf}=21.
    \item jactype = "bandusr" : a banded jacobian, specified by user function
      \code{jacfunc}; the size of the bands specified by \code{bandup} and
      \code{banddown}, corresponds to \code{mf}=24.
    \item jactype = "bandint" : a banded jacobian, calculated by \code{lsode};
      the size of the bands specified by \code{bandup} and \code{banddown},
      corresponds to \code{mf}=25.
  }
  More options are available when specifying \bold{mf} directly. 

  The legal values of \code{mf} are 10, 11, 12, 13, 14, 15, 20, 21, 22, 23,
  24, 25.

  \code{mf} is a positive two-digit integer, \code{mf} = (10*METH + MITER),
  where
  \itemize{
    \item METH indicates the basic linear multistep method: METH = 1 means
      the implicit Adams method. METH = 2 means the method based on backward
      differentiation formulas (BDF-s).
    \item MITER indicates the corrector iteration method: MITER = 0 means
      functional iteration (no Jacobian matrix is involved).
      MITER = 1 means chord iteration with a user-supplied full (NEQ by NEQ)
      Jacobian.
      MITER = 2 means chord iteration with an internally generated
      (difference quotient) full Jacobian
      (using NEQ extra calls to \code{func} per df/dy value).
      MITER = 3 means chord iteration with an internally generated diagonal
      Jacobian approximation
      (using 1 extra call to \code{func} per df/dy evaluation).
      MITER = 4 means chord iteration with a user-supplied banded Jacobian.
      MITER = 5 means chord iteration with an internally generated banded
      Jacobian (using ML+MU+1 extra calls to \code{func} per df/dy evaluation).
  
      If MITER = 1 or 4, the user must supply a subroutine \code{jacfunc}.
  }

  Inspection of the example below shows how to specify both a banded and full
  jacobian. 
  
  The input parameters \code{rtol}, and \code{atol} determine the \bold{error
  control} performed by the solver. 

  See \code{\link{stode}} for details.
   
  \bold{Models} may be defined in compiled C or Fortran code, as well as in
  an R-function. See function \code{\link{stode}} for details.
  
  The output will have the \bold{attributes *istate*, and *rstate*}, two
  vectors with several useful elements. 

  if \code{verbose} = TRUE, the settings of istate and rstate will be written
  to the screen.
      
  the following elements of \bold{istate} are meaningful:
  \itemize{
    \item el 1 : gives the conditions under which the last call to the
      integrator returned.
      2  if lsode was successful, -1 if excess work done, -2 means excess
      accuracy requested. (Tolerances too small),
      -3 means illegal input detected. (See printed message.),
      -4 means repeated error test failures. (Check all input),
      -5 means repeated convergence failures. (Perhaps bad Jacobian
      supplied or wrong choice of MF or tolerances.),
      -6 means error weight became zero during problem. (Solution
      component i vanished, and atol or atol(i) = 0.)
    \item el 12 : The number of steps taken for the problem so far.
    \item el 13 : The number of evaluations for the problem so far.,
    \item el 14 : The number of Jacobian evaluations  and LU decompositions
      so far.,
    \item el 15 : The method order last used (successfully).,
    \item el 16 : The order to be attempted on the next step.,
    \item el 17 : if el 1 =-4,-5: the largest component in the error vector,
  }

  \bold{rstate} contains the following:
  \itemize{
    \item 1: The step size in t last used (successfully).
    \item 2: The step size to be attempted on the next step.
    \item 3: The current value of the independent variable which the solver
      has actually reached, i.e. the current internal mesh point in t.
    \item 4: A tolerance scale factor, greater than 1.0, computed when a
      request for too much accuracy was detected.
  }
  For more information, see the comments in the original code lsode.f}

\seealso{

 \code{\link{steady}}, for a general interface to most of the steady-state
   solvers

 \code{\link{steady.band}}, to find the steady-state of ODE models with a
   banded Jacobian

 \code{\link{steady.1D}}, \code{\link{steady.2D}},
   \code{\link{steady.3D}} steady-state solvers for 1-D, 2-D and 3-D
   partial differential equations.

 \code{\link{stode}}, iterative steady-state solver for ODEs with full
   or banded Jacobian.

 \code{\link{stodes}}, iterative steady-state solver for ODEs with arbitrary
   sparse Jacobian.

}


\keyword{math}

