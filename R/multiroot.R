## =============================================================================
## multiroot.1D: root of multiple nonlinear equations
## resulting by discretizing (partial) differential equations
## =============================================================================

multiroot.1D <- function (f, start, maxiter=100,
       rtol=1e-6, atol=1e-8, ctol=1e-8, nspec = NULL,
       dimens = NULL, verbose=FALSE, positive=FALSE, names=NULL, ...)  {
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

    Fun <- function (time=0,x,parms=NULL)
      list(f(x,...))

    x <- steady.1D(y=start,time=0,func=Fun,parms=NULL,atol=atol,
              rtol=rtol,ctol=ctol,positive=positive,method="stode",
              nspec=nspec,dimens=dimens,names=names,cyclicBnd=NULL)
    precis <- attr(x,"precis")
    attributes(x)<-NULL

    x        <- unlist(x)
    reffx    <-f(x,...)
    i        <- length(precis)

  names(x) <- names(start)

  return(list(root=x,f.root=reffx,iter=i,estim.precis=precis[length(precis)]))

}

## =============================================================================
## multiroot: root of multiple simultaneous nonlinear equations
## =============================================================================

multiroot <- function(f, start, maxiter=100, rtol=1e-6, atol=1e-8, ctol=1e-8,
                      useFortran=TRUE,positive=FALSE,
                      jacfunc=NULL, jactype="fullint", verbose=FALSE,
                      bandup=1, banddown=1, ...)  {

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

  if (useFortran) {
    Fun <- function (time=0,x,parms=NULL)
      list(f(x,...))
    JacFunc <- jacfunc
    if (!is.null (jacfunc))
      JacFunc <-  function (time=0,x,parms=NULL)
         jacfunc(x,...)
    method <- "stode"
    if (jactype=="sparse") {
      method <- "stodes"
      if (! is.null (jacfunc)) stop("jacfunc can not be used when jactype='sparse'")
      x <- stodes(y=start,time=0,func=Fun,parms=NULL,atol=atol,positive=positive,
              rtol=rtol,ctol=ctol,maxiter=maxiter, verbose=verbose)
    }else
      x <- steady(y=start,time=0,func=Fun,parms=NULL,atol=atol,positive=positive,
              rtol=rtol,ctol=ctol,maxiter=maxiter,method=method,
              jacfunc=JacFunc,jactype=jactype,verbose=verbose,bandup=bandup,
              banddown=banddown)
    precis <- attr(x,"precis")
    attributes(x)<-NULL

    x        <- unlist(x)
    reffx    <-f(x,...)
    i        <- length(precis)
    
  } else {  # simple R-implementation - ignores the Jacobian settings

    precis   <- NULL

    x        <- start
    jacob    <- matrix(nrow=N,ncol=N,data=0)
    reffx    <- f(x,...)     # function value,
  
    if (length (reffx) != N)
      stop("'f', function must return as many function values as elements in start")


    for (i in 1:maxiter) {
      refx   <- x
      pp     <- mean(abs(reffx)) # check convergence...
      precis <- c(precis,pp)
      ewt    <- rtol*abs(x)+atol
      if (max(abs(reffx/ewt))<1) break

      # estimate jacobian
      delt   <- perturb(x)

      for (j in 1:N) {
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
    
    } # end for
  } # end fortran/R
  names(x) <- names(start)

  return(list(root=x,f.root=reffx,iter=i,estim.precis=precis[length(precis)]))
}  # end multiroot

