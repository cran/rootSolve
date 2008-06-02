### steady, steady.1D, steady.band -- solves the steady-state condition of 
### ordinary differential equation systems 
### steady is designed for arbitrary sparse systems with a full Jacobian
### steady.1D is designed for solving multi-component 1-D reaction-transport models 
### steady.band is designed for solving single-component 1-D reaction-transport models
### these functions have similar calling sequence as integration routines from
### package deSolve

steady  <- function (y,
                     time=0,
                     func,
                     parms=NULL,
                     method="stode",
                     ...)
{
 if (method=="stode") stode(y,time,func,parms=parms,...) else
 if (method=="stodes") stodes(y,time,func,parms=parms,...) else
 if (method=="runsteady") runsteady(y,time,func,parms=parms,...) 
  
}

steady.1D    <- function (y,
                       time=0,
                       func,
                       parms=NULL,
                       nspec = NULL,
                       dimens = NULL,
                       method="stode",
                       ...)
{
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
     stop ("cannot run steady.1D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens) && is.null(nspec)) 
     stop ("cannot run steady.1D: either nspec or dimens should be specified")
  N     <- length(y)
  if (is.null(nspec)  ) nspec = N/dimens
  else if (N%%nspec != 0) 
     stop("cannot run steady.1D: nspec is not an integer fraction of number of state variables")
   
  if (method=="stodes")
  {
    dimens <- N/nspec
    out <- stodes(y=y,time=time,func=func,parms=parms,
                  nnz=c(nspec,dimens),sparsetype="1D",...)                    
   } else if (is.character(func))
  {
  ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec

  if (method!="stode") stop ("cannot run steady.1D: method should be one of stode or stodes if func is a DLL")
    out <- stode(y=y[ii],time=time,func=func,parms=parms,
                jactype="1Dint",bandup=nspec,banddown=N/nspec,...)                    
                
  out[[1]][ii] <- out[[1]] 
  } else {

  # internal function #
  bmodel <- function (time,state,pars,model,...)
  {
     Modconc <-  model(time,state[ij],pars,...)   # ij: reorder state variables
     c(list(Modconc[[1]][ii]),Modconc[-1])        # ii: reorder rate of change     
  }

  ii    <- as.vector(t(matrix(ncol=nspec,1:N)))   # from ordering per slice -> per spec
  ij    <- as.vector(t(matrix(nrow=nspec,1:N)))   # from ordering per spec -> per slice
    
  bmod  <- function(time,state,pars,...) {bmodel(time,state,pars,func,...)}
  if (method=="stode")
    out <- stode(y[ii],time,func=bmod,parms=parms,
                 bandup=nspec,banddown=nspec,jactype="bandint",...) else
    out <- runsteady(y[ii],time,func=bmod,parms=parms,
                 bandup=nspec,banddown=nspec,jactype="bandint",...) 
                 
                 
  out[[1]][ii] <- out[[1]]
  }
  return(out)
}

steady.2D    <- function (y,
                       time=0,
                       func,
                       parms=NULL,
                       nspec=NULL,
                       dimens,
                       ...)
{
  if (any(!is.na(pmatch(names(list(...)), "jacfunc")))) 
     stop ("cannot run steady.2D with jacfunc specified - remove jacfunc from call list")
  if (is.null(dimens)) 
     stop ("cannot run steady.2D: dimens should be specified")
  if (length(dimens)!=2) 
     stop ("cannot run steady.2D: dimens should contain 2 values")
  N     <- length(y)
  if (N%%prod(dimens) != 0) 
     stop("cannot run steady.2D: dimensions are not an integer fraction of number of state variables")
  if (is.null(nspec)) 
     nspec = N/prod(dimens)
  else if (nspec * prod(dimens) != N) 
     stop("cannot run steady.2D: dimens[1]*dimens[2]*nspec is not equal to number of state variables")
  

  out <- stodes(y=y,time=time,func=func,parms=parms,
                nnz=c(nspec,dimens),sparsetype="2D",...)                    
  return(out)
}

steady.band  <- function (y,
                          time=0,
                          func,
                          parms=NULL,
                          nspec=NULL,
                          bandup=nspec,
                          banddown=nspec,
                          ...)
{
if (is.null(bandup)  ) stop ("cannot run steady.band: bandup is not specified")
if (is.null(banddown)) stop ("cannot run steady.band: banddown is not specified")

  stode(y,time,func,parms=parms,
        bandup=bandup,banddown=banddown,jactype="bandint",...) 
}

### stode -- solves for the root (steady-state) of
###       ordinary differential equation systems defined in func
###       'y' contains the initial guesses for the state variables
###       'parms' is a vector of parameters for func.  They should not
###       change during the iterations. `rtol', and `atol'
###       are, respectively, the relative tolerance parameter, and the
###       absolute tolerance parameter.  `atol' may be scaler or vector.
###       `rtol' is a scaler or a vector.
###
###       The return value is a vector with the last value of the steady-state.
###       iteration. If attribute "steady" is true, this is 
###       the steady-state condition 
###
###      'func' may be a string instead of an R function.  If
###       so, then if jacfunc is not NULL, it must be a character string
###       as well.  In these cases, 'func' is the name
###       of a function to be found in the dll named 'dllname' 
###       (without extension). 'jacfunc' points to the name of the jacobian.
###       initfunc is the name of the function that is the initializer for the problem.
### the implementation of this function is very similar to the implementation of 
### the integration routines in the deSolve package.
   
stode         <- function(y,              # state variables
                          time=0,           # time at which output is wanted  
                          func,             # function that returns rate of change
                          parms=NULL,       # other parameters passed to func and jacfunc                        
                          rtol=1e-6,        # relative tolerance  
                          atol=1e-8,        # absolute tolerance 
                          ctol=1e-8,        # minimal change in dy 
                          jacfunc=NULL,     # jacobian 
                          jactype = "fullint",  # jacobian
                          verbose=FALSE,    # 
                          bandup=1,         # upper band
                          banddown=1,       # lower band
                          positive = FALSE,
                          maxiter=100,    # maximal number of steps during one call to the solver
                          ynames=TRUE,      # if false: names of state variables are not passed to function func
                          dllname=NULL,     # 
                          initfunc=dllname, # 
                          initpar=parms,    # to initialise common block/global variables                          
                          rpar=NULL,           #  
                          ipar=NULL,          #
                          nout=0,           # only if dllname is present: the number of output variables
                          outnames = NULL,  #
                          ...)              # accessory parameters passed to ??
{
### check input
    if (!is.numeric(y))
        stop("`y' must be numeric")
    n <- length(y)
    if (! is.null(time)&&!is.numeric(time))
        stop("`time' must be NULL or numeric")
    if (!is.function(func) && !is.character(func))
        stop("`func' must be a function or character vector")
    if (is.character(func) && (is.null(dllname) || !is.character(dllname)))
        stop("You need to specify the name of the dll or shared library where func can be found (without extension)")
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
    if (!is.null(jacfunc) && !(is.function(jacfunc) || is.character(jacfunc)))
        stop("`jacfunc' must be a function or character vector")
    if (length(atol) > 1 && length(atol) != n)
        stop("`atol' must either be a scalar, or as long as `y'")
    if (length(rtol) > 1 && length(rtol) != n)
        stop("`rtol' must either be a scalar, or as long as `y'")
    if (length(ctol) > 1)
        stop("`ctol' must be a scalar")
    itol <- 1    # length atol and rtol ==1
    if (length(atol)==n && length(rtol)==n) itol <- 4 else
    if (length(atol)==n && length(rtol)!=n) itol <- 2 else 
    if (length(atol)!=n && length(rtol)==n) itol <- 3 
    
### Jacobian, method flag
       if (jactype == "fullint" ) imp <- 22 # full jacobian, calculated internally
  else if (jactype == "fullusr" ) imp <- 21 # full jacobian, specified by user function
  else if (jactype == "bandusr" ) imp <- 24 # banded jacobian, specified by user function
  else if (jactype == "bandint" ) imp <- 25 # banded jacobian, specified internally
  else if (jactype == "1Dint"   ) imp <- 0  # banded jacobian, specified+rearranged internally  
  else stop("jactype must be one of fullint, fullusr, bandusr, or bandint")

  if (imp == 0)
  {
   nspec<-bandup
   ndim <-banddown
   banddown <- nspec
  } else
  {
   nspec<-0
   ndim <-0
  }
  # check if jacfunc is specified if it is needed. 
  if (imp %in% c(21,24) && is.null(jacfunc)) 
  stop ("stode: cannot estimate steady-state: *jacfunc* NOT specified; either specify *jacfunc* or change *jactype*")

  if (imp %in% c(24,25)) nabd <- 1+2*bandup+banddown  else nabd <- n

### print to screen...
  if (verbose)
  {
   print("Steady-state settings")
   if (is.character(func)) print(paste("model function a DLL: ",func))
   if (is.character(jacfunc)) print(paste ("jacobian specified as a DLL: ",jacfunc))
   print("jacobian method")
   df   <- c("method flag, mf","jsv", "meth","miter","itask")
   if (imp==22)txt<-"full jacobian, calculated internally" else 
   if (imp==21)txt<-"full jacobian, specified by user function" else 
   if (imp==24)txt<-"banded jacobian, specified by user function" else 
   if (imp==0 )txt<-"banded jacobian, 1-D, specified internally" else 
               txt<-"banded jacobian, calculated internally"    
   print(data.frame("implicit method", value=imp,message=txt))
  }

### model and jacobian function
    Ynames <- attr(y,"names")
    JacFunc <- NULL    
    ModelInit <- NULL
    if(!is.null(dllname))
    {
        if (is.loaded(initfunc, PACKAGE = dllname,
            type = "") || is.loaded(initfunc, PACKAGE = dllname,
            type = "Fortran")) 
            ModelInit <- getNativeSymbolInfo(initfunc, PACKAGE = dllname)$address
     }

    ## If func is a character vector, then
    ## copy its value to funcname 
    ## check to make sure it describes
    ## a function in a loaded dll
    if (is.character(func)) {
        funcname <- func
      ## get the pointer and put it in func        
       if(is.loaded(funcname, PACKAGE = dllname)) {
       Func <- getNativeSymbolInfo(funcname, PACKAGE = dllname)$address
        } else stop(paste("cannot calculate steady-state: dyn function not loaded: ",funcname))

      ## is there a jacobian?
       if (!is.null(jacfunc)) {
           if (!is.character(jacfunc))
              stop("If 'func' is dynloaded, so must 'jacfunc' be")
           jacfuncname <- jacfunc
           if(is.loaded(jacfuncname, PACKAGE = dllname))
           {JacFunc <- getNativeSymbolInfo(jacfuncname, PACKAGE = dllname)$address
           } else stop(paste("cannot calculate steady-state: jacobian function not loaded ",jacfunc))
        }

      ## If we go this route, the number of "global" results is in nout        
      Nglobal <- nout
      rho     <- NULL
      if (is.null(outnames))
         { Nmtot   <- NULL} else
      if (length(outnames) == nout) 
         { Nmtot   <- outnames} else
      if (length(outnames) > nout) 
         Nmtot <- outnames[1:nout] else
         Nmtot <- c(outnames,(length(outnames)+1):nout)

    }
    else {
      rho <- environment(func)
      # func and jac are overruled, either including ynames, or not
      # This allows to pass the "..." arguments and the parameters
        
      if(ynames)
        {
         Func    <- function(time,state) 
         { attr(state,"names") <- Ynames 
           func   (time,state,parms,...)[1]}   
         
         Func2   <- function(time,state) 
         { attr(state,"names") <- Ynames
           func   (time,state,parms,...)}    
         
         JacFunc <- function(time,state) 
         { attr(state,"names") <- Ynames
           jacfunc(time,state,parms,...)}    
        } else {                            # no ynames...
         Func    <- function(time,state) 
           func   (time,state,parms,...)[1] 
        
         Func2   <- function(time,state) 
           func   (time,state,parms,...)    
         
         JacFunc <- function(time,state) 
           jacfunc(time,state,parms,...)    
        }

      ## Call func once to figure out whether and how many "global"
      ## results it wants to return and some other safety checks
        
        tmp <- eval(Func2(time, y), rho) 
        if (!is.list(tmp))
            stop("Model function must return a list\n")
        if (length(tmp[[1]]) != length(y))
            stop(paste("The number of derivatives returned by func() (",
                length(tmp[[1]]), "must equal the length of the initial conditions vector (",
                length(y), ")", sep = ""))

      # use "unlist" here because some output variables are vectors/arrays
        Nglobal <- if (length(tmp) > 1)   
            length(unlist(tmp[-1]))  else 0
        Nmtot <- attr(unlist(tmp[-1]),"names")

      if (imp %in% c(21,24))
      {tmp <- eval(JacFunc(time, y), rho) 
       if (!is.matrix(tmp)) stop("Jacobian function must return a matrix\n")
       dd <- dim(tmp)
       if((imp ==24 && dd != c(bandup+banddown+1,n)) ||
          (imp ==21 && dd != c(n,n))) stop("Jacobian dimension not ok") 
      } 

    }
    
### calling solver

    storage.mode(y) <- "double"
    storage.mode(rtol) <- storage.mode(atol) <- storage.mode(ctol) <- "double"
    if (is.null(ipar)) ipar<-0
    if (is.null(rpar)) rpar<-0
                  
    out <- .Call("call_dsteady", y, as.double(time), Func, as.double(initpar),    
        ctol, atol, rtol, as.integer(itol), rho,  JacFunc, ModelInit, as.integer(verbose),
        as.integer(imp),as.integer(bandup),as.integer(banddown),as.integer(maxiter), 
        as.integer(positive),as.integer(Nglobal),as.integer(nabd),
        as.integer(nspec),as.integer(ndim),
        as.double (rpar), as.integer(ipar),PACKAGE = "rootSolve")

### saving results    
    precis <- attr(out, "precis")
    steady <- attr(out, "steady")

    attributes(out)<-NULL
    if (Nglobal > 0) {
       if (!is.character(func)) {         # if a DLL: already done...    
            y <- out                      # state variables of this time step
            if(ynames)  attr(y,"names")  <-  Ynames
            out2 <- Func2(time, y)[-1]      
            out <- c(list(y=out), out2)                     
        } else out <- list(y=out[1:n],var=out[(n+1):(n+Nglobal)])
    } else out <- list(y=out)

    attr(out, "precis") <- precis
    attr(out, "steady") <- (steady==1   )
    if (verbose)
    {
      print("precision at each steady state step")
      print(precis)    }
    return(out)

}


                                                                                                                                                                                                            