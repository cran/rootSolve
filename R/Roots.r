
################################################################################
##                                                                            ##
##   routines that find the root of a nonlinear function                      ##
##                                                                            ##
## jacobian.full : generates a full jacobian matrix by numerical differencing ##
## jacobian.band : multidiagonal (banded) jacobian matrix by differencing     ##
## multiroot  : root of multiple simultaneous nonlinear equations             ##
##                                                                            ##
## internal functions - not to be included in package:                        ##
## -----------------------------                                              ##
## Perturb     : adds the numerical differencing value to a value             ##
##                                                                            ##
################################################################################

# requires: Inverse (banded)
##############################################################################
##                                                                          ##
## NONLINEAR INVERSE MODELLING                                              ## 
##    VARIOUS UTILITIES                                                     ##
##                                                                          ##
##############################################################################
                                                                                                       


################################################################################
## Perturb     : adds the numerical differencing value to a value             ##
################################################################################

# internal function #
perturb <- function (value)   # value to be perturbed

#------------------------------------------------------------------------
# Calculates the numerical difference value  - internal function
#------------------------------------------------------------------------

{
    # A small, positive value, not too small
    pmax(abs(value) * 1e-8,1e-8)

}


################################################################################
## uniroot.all: multiple roots of one nonlinear equation                         ##
################################################################################

uniroot.all <- function(
        f,                    # the function for which the root is sought.
        interval,             # a vector containing the end-points of the interval to be searched for the root.
        lower= min(interval), # the lower end point of the interval to be searched.
        upper= max(interval), # the upper end point of the interval to be searched.
        tol= .Machine$double.eps^0.2,      # the desired accuracy (convergence tolerance).
        maxiter= 1000,       # the maximum number of iterations.
        n = 100,             #number of subintervals in which root is sought
        ... )               #additional named or unnamed arguments to be passed to f (but beware of partial matching to other arguments).
{
# error checking as in uniroot...
if (!missing(interval) && length(interval) != 2) 
        stop("'interval' must be a vector of length 2")
if (!is.numeric(lower) || !is.numeric(upper) || lower >= 
        upper) 
        stop("lower < upper  is not fulfilled")
 
# subdivide interval in n subintervals and estimate the function values
xseq <- seq(lower,upper,len=n+1)
mod  <- f(xseq,...)
# some function values may already be 0
Equi <- xseq[which(mod==0)]

ss   <- mod[1:n]*mod[2:(n+1)]  # interval where functionvalues change sign
ii   <- which(ss<0)

for (i in ii) Equi <- c(Equi,uniroot(f,lower=xseq[i],upper=xseq[i+1],...)$root)

return(Equi)
}


################################################################################
## multiroot: root of multiple simultaneous nonlinear equations               ##
################################################################################

multiroot <- function(f,              # function for which the root is sought
                      start,          # vector containing initial guesses for the root
                      maxiter=100,    # maximal number of iterations
                      rtol=1e-6,        # relative tolerance  
                      atol=1e-8,        # absolute tolerance 
                      ctol=1e-8,        # minimal change in dy 
                      useFortran=TRUE,
                      ...)            # additional arguments passed to function 'f'
 {
    N        <- length(start)  
    if (!is.numeric(start)) 
        stop("start conditions should be numeric")
    if (!is.numeric(maxiter))
        stop("`maxiter' must be numeric")
    if (as.integer(maxiter) < 1)
        stop ("maxiter must be >=1")    
    if (!is.numeric(rtol))
        stop("`rtol' must be numeric")
    if (!is.numeric(atol))
        stop("`atol' must be numeric")
    if (!is.numeric(ctol))
        stop("`ctol' must be numeric")
    if (length(atol) > 1 && length(atol) != N)
        stop("`atol' must either be a scalar, or as long as `start'")
    if (length(rtol) > 1 && length(rtol) != N)
        stop("`rtol' must either be a scalar, or as long as `y'")
    if (length(ctol) > 1)
        stop("`ctol' must be a scalar")

if (useFortran)
{ 
  Fun <- function (time=0,x,parms=NULL,...)
     list(f(x,...))
 
  x <- steady(y=start,time=0,func=Fun,parms=NULL,atol=atol,
              rtol=rtol,ctol=ctol,jacfunc="fullint",maxiter=maxiter)
  precis <- attr(x,"precis")
  attributes(x)<-NULL

  x        <- unlist(x)
  reffx    <-f(x,...)
  i        <- length(precis)
    
} else {  # simple R-implementation

  precis   <- NULL

  x        <- start
  jacob    <- matrix(nrow=N,ncol=N,data=0)  
  reffx    <- f(x,...)     # function value, 
  
  if (length (reffx) != N) 
        stop("'f', function must return as many function values as elements in start")


  for (i in 1:maxiter)
   {
    refx   <- x     
    pp     <- max(abs(reffx)) # check convergence... 
    precis <- c(precis,pp)
    ewt    <- rtol*abs(x)+atol
    if (max(abs(reffx/ewt))<1) break

    # estimate jacobian
    delt   <- perturb(x)

       for (j in 1:N)
          {
           x[j] <- x[j]+delt[j]   # Perturb the state variables one by one
           fx   <-  f(x,...)      # recalculate function value  

           # impact of the current state var on rate of change of all state vars
           jacob [,j] <- (fx-reffx)/delt[j]             

           x[j] <- refx[j]   # restore
           }

   # new estimate 
    relchange <- as.numeric(solve(jacob,-1*reffx))
    if (max(abs(relchange)) < ctol) break
    x  <- x + relchange
    reffx  <- f(x,...)     # function value, 
    
   }
 
 }
  names(x) <- names(start)
  return(list(root=x,f.root=reffx,iter=i,estim.precis=precis[length(precis)]))
}  # end multiroot


################################################################################
## gradient  : generates a full jacobian matrix by numerical differencing     ##
################################################################################

gradient<- function(f,    # function returning a set of function values, as a vector
                    x,    # variables
                    ...)  # additional arguments passed to function "f"...)              
{

#------------------------------------------------------------------------
# Given a set of variables, a func, 
# Estimates the gradient matrix by successively perturbing the variables
#------------------------------------------------------------------------
# Reference value of variables and function value
       if (!is.numeric(x)) stop("x-values should be numeric")
       
       refx <- x
       reff <- f(x,...)

       Nx <- length(x)
       Nf <- length(reff)
# Perturb the state variables one by one
       delt   <- perturb(x)
       jacob  <- matrix(nrow=Nf,ncol=Nx,data=0)
       for (j in 1:Nx)
          {

           x[j] <- x[j]+delt[j]

           # recalculate model function value
           newf  <- f(x,...)

           # impact of the current variable on function values
           jacob [,j] <- (newf-reff)/delt[j]             

           x[j] <- refx[j]   # restore
           }
        colnames(jacob) <- names(x)   
        return(jacob)

} ## END gradient

################################################################################
## jacobian.full  : generates a full jacobian matrix by numerical differencing    ##
################################################################################

jacobian.full<- function(y,             # (state) variables
                         func,          # function that calculates rate of change
                         dy=NULL,       # reference rate of change
                         time=0,        # time passed to function 'func'
                         parms=NULL,    # parameter values passed to function 'func'
                         ...)           # other arguments passed to function 'func'
{

#------------------------------------------------------------------------
# Given a set of state variables, a model, 
# Estimates the Jacobian matrix by successively perturbing the states
#------------------------------------------------------------------------
# Reference value of state variable and of rate of change (dC)
       if (!is.numeric(y)) stop("y-values should be numeric")
       N    <- length(y)
       refy <- y
       ifelse (is.null(dy), refdy <- unlist( func(time,y,parms,...))[1:N], refdy <- dy)
       if (! is.numeric(refdy)) stop("dy-values should either be NULL or numeric")
       if (length(refdy) != N)  stop("function should return at least one value for each y")
       ynames  <- attr(y,"names")
 
# Perturb the state variables one by one
       delt   <- perturb(y)

       ny <-length(y)
       jacob  <- matrix(nrow=ny,ncol=ny,data=0)
       for (j in 1:ny)
          {

           y[j] <- y[j]+delt[j]

           # recalculate model rate of change    
           dy  <-  unlist( func(time,y,parms,...))[1:N]

           # impact of the current state var on rate of change of all state vars
           jacob [,j] <- (dy-refdy)/delt[j]             

           y[j] <- refy[j]   # restore
           }
        colnames (jacob) <- ynames          
           
        return(jacob)

} ## END jacobian.full


################################################################################
## jacobian.band  : multidiagonal (banded) jacobian matrix by differencing        ##
################################################################################

jacobian.band<- function(y,         # (state) variables
                     func,          # function that calculates rate of change
                     bandup=1,      # number of nonzero bands above
                     banddown=1,    # number of nonzero bands below diagonal
                     dy=NULL,       # reference rate of change 
                     time=0,        # time passed to function 'func'
                     parms=NULL,    # parameter values passed to function 'func'
                     ...)           # other arguments passed to function 'func'

{

#------------------------------------------------------------------------
# Given a set of state variables, a model, 
# Estimates the Jacobian matrix by successively perturbing the states
# Assumes banded jacobian matrix; 
# bandup, banddwon: number of bands above and below diagonal
#------------------------------------------------------------------------
# Reference value of state

   if (!is.numeric(y)) stop("y-values should be numeric")   
   ny <- length(y)
   if(is.null(dy)) dy <-  unlist( func(time,y,parms,...))[1:ny]
   if (! is.numeric(dy)) stop("dy-values should either be NULL or numeric")
   if (length(dy) != ny) stop("function should return at least one value for each y")
   
   ynames  <- attr(y,"names")

   refy  <- y
   refdy <- dy

# Perturb the state variables, assume banded structure (only 3*nspec iterations)

   nband  <- bandup+banddown+1
   jacob  <- matrix(nrow=nband,ncol=ny,data=0)
   delt   <- perturb(y)        # perturbation factors
   for ( j in 1:nband)
     {
      kpert    <- seq(j,ny,nband)                 # list of state var to perturb
      y[kpert] <- y[kpert] + delt[kpert]   # perturbed state var

      # new rate of change
      dy <-  unlist( func(time,y,parms,...))[1:ny]
      for (k in kpert)
        { 
          iseq <- seq(max(k-bandup,1),min(k+banddown,ny))
          # impact of the selected states on the rate of change of all states
          jacob [iseq-k+bandup+1,k] <- (dy[iseq] -refdy[iseq]      )/delt[k]
        }
      y<-refy
    } 
    colnames (jacob) <- ynames          
    return(jacob )   # jacobian matrix, banded format 

} ## END jacobian.band

