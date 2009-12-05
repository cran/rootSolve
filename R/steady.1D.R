## =============================================================================
## steady.1D -- solves the steady-state condition of
## ordinary differential equation systems resulting from
## multi-component 1-D PDE models
## has a similar calling sequence as ode.1D from package deSolve
## =============================================================================

steady.1D    <- function (y, time=0, func, parms=NULL, nspec = NULL,
                       dimens = NULL, names = NULL, method="stode",
                       cyclicBnd = NULL, ...) {
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
    stop ("cannot run steady.1D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens) && is.null(nspec)) 
    stop ("cannot run steady.1D: either nspec or dimens should be specified")
  N     <- length(y)
  if (is.null(nspec)  )
    nspec <- N/dimens
  else if (N%%nspec != 0) 
    stop("cannot run steady.1D: nspec is not an integer fraction of number of state variables")
  if (! is.null(names) && length(names) != nspec)
    stop("length of names should equal nspec")
  if (! is.null (cyclicBnd)) method <- "stodes"
  
  if (nspec == 1 & method == "stode") {
    out <- steady.band(y, time, func, parms, nspec, ...)
    if (! is.null(names)) {
      out[[1]] <- matrix(ncol=nspec,data=out[[1]])
      colnames(out[[1]]) <- names
    }
    return(out)

  }
  if (method=="stodes") {
    dimens <- N/nspec
    Bnd <- 0
    if (! is.null(cyclicBnd)) {
    if (max(cyclicBnd) > 1 )
      stop ("cannot run steady.1D: cyclicBnd should be NULL or a value not exceeding 1")
    Bnd <-1
  }

    out <- stodes(y=y,time=time,func=func,parms=parms,
                  nnz=c(nspec,dimens,Bnd),sparsetype="1D",...)
  } else if (is.character(func)) {
    ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec

    if (method!="stode")
      stop ("cannot run steady.1D: method should be one of stode or stodes if func is a DLL")
    out <- stode(y=y[ii],time=time,func=func,parms=parms,
                jactype="1Dint",bandup=nspec,banddown=N/nspec,...)                    
                
    out[[1]][ii] <- out[[1]]
  } else {

  # internal function #

    bmodel <- function (time,state,pars,model,...) {
      Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
      c(list(Modconc[[1]][ii]),Modconc[-1])        # ii: reorder rate of change
    }

    ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec
    ij    <- as.vector(t(matrix(nrow=nspec,1:N)))   # from ordering per spec -> per slice
    
    bmod  <- function(time,state,pars,...) {
      bmodel(time,state,pars,func,...)
    }

    if (method=="stode")
      out <- stode(y[ii],time,func=bmod,parms=parms,
                bandup=nspec,banddown=nspec,jactype="bandint",...) else
      out <- runsteady(y[ii],times=time,func=bmod,parms=parms,
                bandup=nspec,banddown=nspec,jactype="bandint",...)

    out[[1]][ii] <- out[[1]]
    if (! is.null(names)) {
      out[[1]] <- matrix(ncol=nspec,data=out[[1]])
      colnames(out[[1]]) <- names
    }
  }
  class(out) <- c("steady1D","list")    # a steady-state 
  attr(out,"dimens") <- dimens
  attr(out, "nspec") <- nspec
  attr(out,"ynames") <- names
  return(out)
}
