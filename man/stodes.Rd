\name{stodes}
\alias{stodes}
\title{Steady-state solver for ordinary differential equations (ODE) with a sparse jacobian.}
\description{Estimates the steady-state condition for a system of
  ordinary differential equations (ODE) in the form:
  \deqn{dy/dt = f(t,y)}

  and where the jacobian matrix df/dy has an arbitrary sparse structure.

  Uses a newton-raphson method, implemented in Fortran.

  The system of ODE's is written as an \R function or defined in
  compiled code that has been dynamically loaded.}
\usage{stodes(y, time=0, func, parms=NULL, 
rtol=1e-6, atol=1e-8, ctol=1e-8, sparsetype="sparseint",verbose=FALSE,
nnz=NULL, inz=NULL, lrw=NULL, ngp=NULL, positive=FALSE, maxiter=100, 
ynames=TRUE, dllname=NULL, initfunc=dllname, initpar=parms, 
rpar=NULL, ipar=NULL, nout=0, outnames = NULL, ...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ode system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{time }{time for which steady-state is wanted; the default is \code{time}=0}
  \item{func }{either a user-supplied function that computes the values of the
    derivatives in the ode system (the \emph{model definition}) at time \code{time}, or a character string
    giving the name of a compiled function in a dynamically loaded
    shared library.

    If \code{func}  is a user-supplied function, it must be called as:
    \code{yprime = func(t, y, parms)}.  \code{t} is the time point
    at which the steady-state is wanted, \code{y} is the current estimate of the variables
    in the ode system.  If the initial values \code{y} has a names
    attribute, the names will be available inside \code{func}.  \code{parms} is
    a vector of parameters (which may have a names attribute).

    The return value of \code{func} should be a list, whose first element is a
    vector containing the derivatives of \code{y} with respect to
    \code{time}, and whose next elements (possibly with a
    \code{\link{names}} attribute) are global values that are required as
    output.
        If \code{func} is a string, then \code{dllname} must give the name
    of the shared library (without extension) which must be loaded
    before \code{stodes()} is called.  see Details for more information.}
  \item{parms }{other parameters passed to \code{func} }
  \item{rtol }{relative error tolerance, either a scalar or a vector, one value for each y }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one value for each y}
  \item{ctol }{if between two iterations, the maximal change in y is less than this amount, steady-state is reached}
  \item{sparsetype }{the sparsity pattern, to date only "sparseint", sparse jacobian, estimated internally by stodes}
  \item{verbose }{if TRUE: full output to the screen, e.g. will output the steady-state settings}
  \item{nnz }{the number of nonzero elements in the sparse Jacobian (if this is unknown, use an estimate); If NULL, a guess will be made, and
    		if not sufficient, \code{stodes} will return with a message indicating the size actually required.
    \cr If a solution is found, the minimal value of \code{nnz} actually required is returned by the solver (1st element of attribute \code{dims})}
  \item{inz }{(row,column) indices to the nonzero elements in the sparse Jacobian. If this is NULL, the 
    sparsity will be determined by \code{stodes}}
  \item{lrw }{the length of the work array of solver; due to the sparsicity, this cannot be readily predicted. If NULL, a guess will be made, and
    		if not sufficient, \code{stodes} will return with a message indicating the size actually required.
    Therefore, some experimentation may be necessary to estimate the value of \code{lrw} 
    \cr If a solution is found, the minimal value of \code{lrw} actually required is returned by the solver (3rd element of attribute \code{dims})}
  \item{ngp }{number of groups of independent state variables. Due to the sparsicity, this cannot be readily predicted. If NULL, a guess will be made, and
		if not sufficient, \code{stodes} will return with a message indicating the size actually required.
    Therefore, some experimentation may be necessary to estimate the value of \code{ngp} 
    \cr If a solution is found, the minimal value of \code{ngp} actually required is returned by the solver (2nd element of attribute \code{dims}    }
  \item{positive }{if TRUE, the state variables are forced to be non-negative numbers}
  \item{maxiter }{maximal number of iterations during one call to the solver}
  \item{ynames }{if FALSE: names of state variables are not passed to function \code{func} ; this may speed up the simulation especially for multi-D models}
  \item{dllname }{a string giving the name of the shared library (without extension) that contains all the compiled function or subroutine definitions referred to in \code{func}.}
  \item{initfunc }{if not NULL, the name of the initialisation function (which initialises values of parameters), as provided in \file{dllname}. See details.}
  \item{initpar }{only when \file{dllname} is specified and an initialisation function \code{initfunc} is in the dll: the parameters passed to the initialiser, to initialise the common blocks (FORTRAN) or global variables (C, C++)}
  \item{rpar }{only when \file{dllname} is specified: a vector with double precision values passed to the dll-functions whose names are specified by \code{func}}
  \item{ipar }{only when \file{dllname} is specified: a vector with integer values passed to the dll-functions whose names are specified by \code{func}}
  \item{nout }{only used if \file{dllname} is specified: the number of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{outnames }{only used if \file{dllname} is specified and \code{nout} > 0: the names of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{... }{additional arguments passed to \code{func} allowing this to be a generic function}
}
\value{
  A list containing 
  \item{y }{ A vector with the state variable values from the last iteration during estimation of steady-state condition of the system of equations.     
    If \code{y} has a names attribute, it will be used to label the output values. }
  \item{... }{ the number of "global" values returned}
  The output will have the attribute \code{steady}, which returns TRUE, if steady-state has been reached and the attribute
    \code{precis} with the precision attained during each iteration.
}
\author{Karline Soetaert <k.soetaert@nioo.knaw.nl>}
\examples{# 1000 simultaneous equations
model <- function (time,OC,parms,decay,ing)
{
 # model describing C in a sediment,
 # Upper boundary = imposed flux, lower boundary = zero-gradient
 Flux  <- v * c(OC[1] ,OC) +              # advection
          -Kz*diff(c(OC[1],OC,OC[N]))/dx  # diffusion;
 Flux[1]<- flux     # imposed flux

 # Rate of change= Flux gradient and first-order consumption
 dOC   <- -diff(Flux)/dx - decay*OC

 # Fraction of OC in first 5 layers is translocated to mean depth
 # (layer N/2)
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
ss   <- stodes(runif(N),func=model,parms=NULL,
               positive=TRUE, decay=5,ing=20,verbose=TRUE)

plot(ss$y[1:N],dist,ylim=rev(range(dist)),type="l",lwd=2,
     xlab="Nonlocal exchange",ylab="sediment depth",
     main="stodes, sparse jacobian")}
\references{For a description of the Newton-Raphson method, e.g.
\cr Press, WH, Teukolsky, SA, Vetterling, WT, Flannery, BP, 1996. Numerical Recipes in FORTRAN. The Art of Scientific computing. 2nd edition. Cambridge University Press.

The algorithm uses linear algebra routines from the Yale sparse matrix package: 
\cr Eisenstat, S.C., Gursky, M.C., Schultz, M.H., Sherman, A.H., 1982. Yale Sparse Matrix Package. 
i. The symmetric codes. Int. J. Num. meth. Eng. 18, 1145-1151.}
\details{The work is done by a Fortran 77 routine that implements the Newton-Raphson method.

\code{stodes} is to be used for problems, where the Jacobian has a sparse structure. 

There are two choices for the sparsity specification, depending on whether \code{inz} is present.
\item If matrix \code{inz} is present, the sparsity is determined by the user. \code{inz} should contain indices (row, column) to the nonzero elements in the Jacobian matrix.
In this case, \code{nnz} will be set equal to the number of rows in \code{inz}
\item If matrix \code{inz} is NOT present, the sparsity is estimated by the solver, based on numerical differences.
In this case, it is advisable to provide an estimate of the number of non-zero elements in the Jacobian (\code{nnz})
This value can be approximate; upon return the number of nonzero elements actually required will be known (1st element of attribute \code{dims})</nnz>

Either way, the Jacobian itself is always generated by the solver (i.e. there is no provision to provide an analytic Jacobian). 
\cr This is done by perturbing simulataneously a combination of state variables that do not affect each other. 
This significantly reduces computing time. The number of groups with independent state variables can be given by \code{ngp}

The input parameters \code{rtol}, \code{atol} and \code{ctol} determine the \bold{error
control} performed by the solver. See help for \code{stode} for details.

\bold{Models} may be defined in compiled C or Fortran code, as well as in R. See help for \code{stode} for details.}
\seealso{\item \code{\link{stode}}, steady-state solver with full or banded jacobian.
  \item \code{\link{steady.band}}, for steady-state estimation, when the jacobian matrix is banded, and where the state variables need NOT to be rearranged
  \item \code{\link{steady.1D}}, for steady-state estimation, when the jacobian matrix is banded, and where the state variables need to be rearranged}
\keyword{math}
