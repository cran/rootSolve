
### stodes -- sparse solver for the root (steady-state) of
###       ordinary differential equation systems
###       Sparse Jacobian
###

stodes        <- function(y,              # state variables
                          time=0,           # time at which output is wanted
                          func,             # function that returns rate of change
                          parms=NULL,       # other parameters passed to func
                          rtol=1e-6,        # relative tolerance
                          atol=1e-8,        # absolute tolerance
                          ctol=1e-8,        # minimal change in dy
                          sparsetype="sparseint", # sparsity type
                          verbose=FALSE,    #
                          nnz=NULL,
                          inz=NULL,
                          lrw=NULL,ngp=NULL,
                          positive = FALSE,
                          maxiter=100,    # maximal number of steps during one call to the solver
                          ynames=TRUE,      # if false: names of state variables are not passed to function func
                          dllname=NULL,     #
                          initfunc=dllname, #
                          initpar=parms,    # to initialise common block/global variables
                          rpar=NULL,           #
                          ipar=NULL,          #
                          nout=0,           # only if dllname is present: the number of output variables
                          outnames=NULL,
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

    Type <- 1            # sparsity to be determined numerically
    ian <- 0
    jan <- 0
    if (is.null(ngp))   ngp = n+1
        
    if(sparsetype=="sparseint")
    {
    if (! is.null(inz))  # sparsity is imposed; create ian, jan 
    {Type <- 0
     nnz <- nrow(inz)
     jan <- numeric(nnz)
     ian <- numeric(n+1)
     iw <- 1
     ian[1] <- 1
   # column indices should be sorted...
     rr  <- inz[,2]
     if (min(rr[2:nnz]-rr[1:(nnz-1)])<0) stop ("cannot proceed: row indices in nnz should be sorted")
     for(i in 1:n)
     {
      ii <- which (rr==i)
      il <- length(ii)
      i1 <- ian[i]
      i2 <- ian[i]+il-1
      ian[i+1] <- i2+1
      if (il>0) jan[i1:i2] <- inz[ii,1]
      }
        
     } else if (is.null(nnz))   nnz = n*n
    } else if (sparsetype == "1D")   {
      Type   <- 2
      nspec  <- nnz[1] 
      nnz    <- c(n*(2+nspec)-2*nspec,nnz)
      ngp    <- 3*nspec+1
    } else if (sparsetype =="2D")    {
      Type   <- 3
      nspec  <- nnz[1] 
      dimens <- nnz[2:3]
      nnz   <- c(n*(4+nspec)-2*nspec*(sum(dimens)),nnz)
      ngp    < 4*nspec+1 

    } else stop("cannot run stodes: sparsetype not known ")

    if (is.null(lrw))   lrw = 3*n+4*nnz[1]



### print to screen...
  if (verbose)
  {
   print("Steady-state settings")
   if (is.character(func)) print(paste("model function a DLL: ",func))
   if (sparsetype=="sparseint")txt<-"sparse jacobian, calculated internally" else 
   if (sparsetype=="1D")    txt<-"sparse 1-D jacobian, calculated internally" else 
   if (sparsetype=="2D")    txt<-"sparse 2-D jacobian, calculated internally" 
   print(data.frame(sparseType = sparsetype, message=txt))

  }

### model and jacobian function
    Ynames <- attr(y,"names")
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

        } else {                            # no ynames...
         Func    <- function(time,state)
           func   (time,state,parms,...)[1]

         Func2   <- function(time,state)
           func   (time,state,parms,...)

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

    }

### calling solver
    imp <- 22
    storage.mode(y) <- "double"
    storage.mode(rtol) <- storage.mode(atol) <- storage.mode(ctol) <- "double"

    out <- .Call("call_stsparse", y, as.double(time), Func, as.double(initpar),
        ctol, atol, rtol, as.integer(itol), rho,  ModelInit, as.integer(verbose),
        as.integer(imp),as.integer(nnz),as.integer(lrw),as.integer(ngp),as.integer(maxiter),
        as.integer(positive),as.integer(Nglobal),
        as.double (rpar), as.integer(ipar), as.integer(Type),
        as.integer(ian),as.integer(jan), PACKAGE = "rootSolve")

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
    attr(out, "steady") <- (steady[1]==1   )
    attr(out, "dims"  ) <- steady[2:4]

    if (verbose)
    {
      print("precision at each steady state step")
      print(precis)    
      print("")
      print("--------------------")
      print(" Memory requirements")
      print("--------------------")      
      nf <- c(" nnz","ngp","nsp") 
      df <- c( " the number of nonzero elements",
               " the number of independent groups of state variables ",
               " the length of the work array actually required."              )

       print(data.frame(par=nf,mess=df, val=steady[2:4]))
      
      }
      
    return(out)

}

