\name{stode}
\alias{stode}
\title{Iterative steady-state solver for ordinary differential equations (ODE) and a full or banded Jacobian.}
\description{Estimates the steady-state condition for a system of 
  ordinary differential equations (ODE) written in the form:
  \deqn{dy/dt = f(t,y)}   
  
  i.e. finds the values of \code{y} for which f(t,y) = 0.
  
  Uses a newton-raphson method, implemented in Fortran 77.
  
  The system of ODE's is written as an \R function or defined in
  compiled code that has been dynamically loaded.}
\usage{stode(y, time=0, func, parms=NULL, 
  rtol=1e-6, atol=1e-8, ctol=1e-8, jacfunc=NULL, 
  jactype="fullint", verbose=FALSE, bandup=1, banddown=1,   
  positive=FALSE, maxiter=100, ynames=TRUE, 
  dllname=NULL, initfunc=dllname, initpar=parms, 
  rpar=NULL, ipar=NULL, nout=0, outnames = NULL,...)}
\arguments{
  \item{y }{the initial guess of (state) values for the ode system, a vector. If \code{y} has a name attribute, the names will be used to label the output matrix.}
  \item{time }{time for which steady-state is wanted; the default is \code{time}=0}
  \item{func }{either a user-supplied function that computes the values of the
    derivatives in the ode system (the \emph{model definition}) at time \code{time}, or a character string
    giving the name of a compiled function in a dynamically loaded
    shared library.

    If \code{func}  is a user-supplied function, it must be called as:
    \code{yprime = func(t, y, parms, ...)}.  \code{t} is the time point
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
    before \code{stode()} is called.  see Details for more information.}
  \item{parms }{other parameters passed to \code{func} and \code{jacfunc}}
  \item{rtol }{relative error tolerance, either a scalar or a vector, one value for each \code{y} }
  \item{atol }{absolute error tolerance, either a scalar or a vector, one value for each \code{y}}
  \item{ctol }{if between two iterations, the maximal change in \code{y} is less than this amount, steady-state is assumed to be reached}
  \item{jacfunc }{if not \code{NULL}, either a user-supplied \R function that estimates
    the Jacobian of the system of differential equations dydot(i)/dy(j), 
    or a character string giving the name of a compiled function in a dynamically loaded shared library as provided in 
    \code{dllname}. In some circumstances, supplying
    \code{jacfunc} can speed up the computations.  The \R calling sequence for
    \code{jacfunc} is identical to that of \code{func}.      
    \cr If the Jacobian is a full matrix, \code{jacfunc} should return a matrix dydot/dy, where the ith
    row contains the derivative of \eqn{dy_i/dt} with respect to \eqn{y_j},
    or a vector containing the matrix elements by columns (the way \R and Fortran store matrices).
    \cr If the Jacobian is banded, \code{jacfunc} should return a matrix containing only the 
    nonzero bands of the jacobian, (dydot/dy), rotated row-wise.}
  \item{jactype }{the structure of the Jacobian, one of "fullint", "fullusr", "bandusr", or "bandint" - either full or banded and estimated internally or by user}
  \item{verbose }{if TRUE: full output to the screen, e.g. will output the steady-state settings}
  \item{bandup }{number of non-zero bands above the diagonal, in case the Jacobian is banded}
  \item{banddown }{number of non-zero bands below the diagonal, in case the jacobian is banded}
  \item{positive }{if TRUE, the state variables \code{y} are forced to be non-negative numbers}
  \item{maxiter }{maximal number of iterations during one call to the solver}
  \item{ynames }{if FALSE: names of state variables are not passed to function \code{func} ; this may speed up the simulation especially for multi-D models}
  \item{dllname }{a string giving the name of the shared library (without extension) that contains all the compiled function or subroutine definitions referred to in \code{func} and \code{jacfunc}.}
  \item{initfunc }{if not NULL, the name of the initialisation function (which initialises values of parameters), as provided in \file{dllname}. See details.}
  \item{initpar }{only when \file{dllname} is specified and an initialisation function \code{initfunc} is in the dll: the parameters passed to the initialiser, to initialise the common blocks (FORTRAN) or global variables (C, C++)}
  \item{rpar }{only when \file{dllname} is specified: a vector with double precision values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{ipar }{only when \file{dllname} is specified: a vector with integer values passed to the dll-functions whose names are specified by \code{func} and \code{jacfunc}}
  \item{nout }{only used if \file{dllname} is specified and the model is defined in compiled code: the number of output variables calculated in the compiled function \code{func}, present in the shared library. Note:
    it is not automatically checked whether this is indeed the number of output variables calculed in the dll - you have to perform this check in the code - see help of daspk or lsoda}
  \item{outnames }{only used if \file{dllname} is specified and \code{nout} > 0: the names of output variables calculated in the compiled function \code{func}, present in the shared library}
  \item{... }{additional arguments passed to \code{func} and \code{jacfunc} allowing this to be a generic function}
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
ss   <- stode(runif(N),func=model,parms=NULL,positive=TRUE, decay=5,ing=20)

plot(ss$y[1:N],dist,ylim=rev(range(dist)),type="l",lwd=2,
xlab="Nonlocal exchange",ylab="sediment depth",main="stode, full jacobian")}
\references{For a description of the Newton-Raphson method, e.g.
\cr Press, WH, Teukolsky, SA, Vetterling, WT, Flannery, BP, 1996. Numerical Recipes in FORTRAN. The Art of Scientific computing. 2nd edition. Cambridge University Press.
\cr The algorithm uses LINPACK code:
\cr Dongarra, J.J., J.R. Bunch,  C.B. Moler and G.W. Stewart, 1979. LINPACK user's guide, SIAM, Philadelphia.}
\details{The work is done by a Fortran 77 routine that implements the Newton-Raphson method. It uses code from LINPACK. 

The form of the \bold{Jacobian} can be specified by \code{jactype} which can take the following values:

    \item jactype = "fullint" : a full jacobian, calculated internally by the solver, the default
    \item jactype = "fullusr" : a full jacobian, specified by user function \code{jacfunc}
    \item jactype = "bandusr" : a banded jacobian, specified by user function \code{jacfunc}; the size of the bands specified by \code{bandup} and \code{banddown}
    \item jactype = "bandint" : a banded jacobian, calculated by the solver; the size of the bands specified by \code{bandup} and \code{banddown}

if \code{jactype}= "fullusr" or "bandusr" then the user must supply a subroutine \code{jacfunc}.

The input parameters \code{rtol}, \code{atol} and \code{ctol} determine the \bold{error
control} performed by the solver. 

The solver will control the vector
  \bold{e} of estimated local errors in \bold{y}, according to an
  inequality of the form max-norm of ( \bold{e}/\bold{ewt} )
  \eqn{\leq}{<=} 1, where \bold{ewt} is a vector of positive error
  weights.  The values of \code{rtol} and \code{atol} should all be
  non-negative.
  The form of \bold{ewt} is:
  
  \deqn{\mathbf{rtol} \times \mathrm{abs}(\mathbf{y}) + \mathbf{atol}}{\bold{rtol} * abs(\bold{y}) + \bold{atol}}
  
  where multiplication of two vectors is element-by-element.

In addition, the solver will stop if between two iterations, the maximal change in the values of \bold{y} is less than \code{ctol}.
  

\bold{Models} may be defined in compiled C or Fortran code, as well as in R. 
\cr If \code{func} or \code{jacfunc} are a string, then they are assumed to be compiled code.
\cr In this case, \code{dllname} must give the name of the shared library (without extension) which must be loaded before lsode() is called. 

\item If \code{func} is a user-supplied \bold{R-function}, it must be called as: yprime = func(t, y, parms,...). 
t is the \code{time} at which the steady-state should be estimated, \code{y} is the current estimate of the variables in the ode system. 
The return value of \code{func} should be a list, whose first element is a vector containing the derivatives of \code{y} with respect to \code{time}, 
and whose next elements contains output variables whose values at steady-state are also required . 
\cr An example is given below:

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
\code{ST <- stode(y=y,fun=model,parms=pars,pos=TRUE))}\cr 
   

  For \bold{code written in C}, the calling sequence for \code{func} must be as follows:
  
  \code{void anoxmod(int *neq, double *t, double *y, double *ydot,}\cr
  \code{double *yout, int *ip)}\cr

  \code{\{}\cr
  \code{double OM, O2, SO4, HS;}\cr
  \code{double Min, oxicmin, anoxicmin;}\cr
  \code{if (ip[0] <1) error("nout should be at least 1");}\cr
  \code{OM  = y[0];}\cr
  \code{O2  = y[1];}\cr
  \code{SO4 = y[2];}\cr
  \code{HS  = y[3];}\cr
  \code{Min       = r*OM;}\cr
  \code{oxicmin   = Min*(O2/(O2+ks));}\cr
  \code{anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2);}\cr

  \code{ydot[0]  = Flux - oxicmin - anoxicmin;}\cr
  \code{ydot[1]  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2);}\cr
  \code{ydot[2]  = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4);}\cr
  \code{ydot[3]  =  0.5*anoxicmin  -rox*HS*(O2/(O2+ks)) + D*(BHS-HS);}\cr

  \code{yout[0] = SO4+HS;}\cr
  \code{\}}\cr

  where \code{*neq} is the number of equations, \code{*t} is the value
  of the independent variable, \code{y} points to a double precision
  array of length \code{*neq} that contains the current value of the
  state variables, and \code{ydot} points to an array that will contain
  the calculated derivatives.
  \cr \code{yout} points to a double precision vector whose first \code{nout} values are other output variables (different from the state variables \code{y}), 
  and the next values are double precision values as passed by parameter \code{rpar} when calling the steady-state solver. The key to the elements of \code{yout} is set in \code{*ip}.
  \cr \code{*ip} points to an integer vector whose length is at least 3; the first element contains the number of output values (which should be equal to \code{nout}), 
  its second element contains the length of \code{*yout}, and the third element contains the length of \code{*ip}; 
  next are integer values, as passed by parameter \code{ipar} when calling the steady-state solver.


  For \bold{code written in Fortran}, the calling sequence for \code{func} must be as in the following
  example:
  
  \code{subroutine model (neq, t, y, ydot, yout, ip)}\cr
  \code{double precision t, y(4), ydot(4), yout(*)}\cr
  \code{double precision OM,O2,SO4,HS}\cr  
  \code{double precision min, oxicmin, anoxicmin}\cr  

  \code{integer neq, ip(*)}\cr
  \code{double precision D, Flux, r, rox, ks, ks2, BO2, BSO4, BHS}\cr  
  \code{common /myparms/D, Flux, r, rox, ks, ks2, BO2, BSO4, BHS}\cr

  \code{IF (ip(1) < 1) call rexit("nout should be at least 1")}\cr
  \code{OM  = y(1)}\cr  
  \code{O2  = y(2)}\cr  
  \code{SO4 = y(3)}\cr  
  \code{HS  = y(4)}\cr  
  \code{Min       = r*OM}\cr  
  \code{oxicmin   = Min*(O2/(O2+ks))}\cr  
  \code{anoxicmin = Min*(1-O2/(O2+ks))* SO4/(SO4+ks2)}\cr  

  \code{ydot(1)  = Flux - oxicmin - anoxicmin}\cr                
  \code{ydot(2)  = -oxicmin      -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)}\cr                
  \code{ydot(3)  = -0.5*anoxicmin  +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4}\cr                
  \code{ydot(4)  =  0.5*anoxicmin  -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)}\cr                        

  \code{yout(1) = SO4+HS}\cr    
  \code{return}\cr
  \code{end}\cr
 
  Note that we start by checking whether enough room is allocated for the output variables, else an error is passed to R (\code{rexit}) and the integration is stopped. 

  In this example, parameters are kept in a common block (called \code{myparms}) in the Fortran code

  In order to put parameters in the common block from the calling R code, an \bold{initialisation subroutine} as specified in \code{initfunc}
  should be defined. This function has as its sole argument a function \code{steadyparms} that fills a double array with double precision
  values.  In the example here, the initialisation subroutine is called \code{myinit}:
  
\code{subroutine myinit(steadyparms)}

\code{external steadyparms}\cr
\code{double precision parms(9)}\cr
\code{common /myparms/parms}

\code{call steadyparms(9, parms)}

\code{return}\cr
\code{end}

   Here \code{myinit} just calls \code{steadyparms} with the dimension of the parameter vector, and the array \code{parms} that will contain the parameter values.

The corresponding C-code is:

\code{void initanox (void (* steadyparms)(int *, double *))}\cr
  \code{\{}\cr
\code{   int N = 9;}\cr
\code{   steadyparms(&N, parms);}\cr
  \code{\}}\cr


  If it is desired to supply a Jacobian to the solver, then the Jacobian
  must be defined in compiled code if the ode system is.  The C function
  call for such a function must be as follows:
  
  \code{void myjac(int *neq, double *t, double *y, int *ml,}\cr
  \code{           int *mu, double *pd, int *nrowpd, double *yout, int *ip)}\cr

The corresponding subroutine call in Fortran is:
      
\code{subroutine myjac (neq, t, y, ml, mu, pd, nrowpd, yout, ip)}\cr
\code{integer neq, ml, mu, nrowpd, ip(*)}\cr
\code{double precision y(*), pd(nrowpd,*), yout(*)}\cr

To run the model using e.g. the Fortran code, the code in anoxmod.f must first be compiled.
\cr This can be done in R itself:  

\code{system("R CMD SHLIB anoxmod.f")}

which will create file anoxmod.dll

After loading the DLL, the model can be solved:

\code{dyn.load("anoxmod.dll")}  

\code{ST2 <- stode(y=y,fun="model",parms=pars,)}\cr
\code{dllname="anoxmod",initfunc="myinit",pos=TRUE,nout=1)}

Examples in both C and Fortran are in the \file{dynload} subdirectory of
the \code{rootSolve} package directory.}
\seealso{\item \code{\link{stodes}}, for steady-state estimation, when the jacobian matrix is sparse.
  \item \code{\link{steady.band}}, for steady-state estimation, when the jacobian matrix is banded, and where the state variables need NOT to be rearranged 
  \item \code{\link{steady.1D}}, for steady-state estimation, when the jacobian matrix is banded, and where the state variables need to be rearranged}
\note{The implementation of \code{stode} and substantial parts of the help file is similar to the implementation of the integration routines (e.g. \code{lsode}) from package deSolve.}
\keyword{math}

