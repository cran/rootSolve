## =============================================================================
## steady, -- solves the steady-state condition of
## ordinary differential equation systems
## has similar calling sequence as integration routines from package deSolve
## =============================================================================

steady  <- function (y, time=0, func, parms=NULL, method="stode", ...)  {

  if (method=="stode")
    stode(y,time,func,parms=parms,...)  else
  if (method=="stodes")
    stodes(y,time,func,parms=parms,...) else
  if (method=="runsteady")
    runsteady(y,times=time,func,parms=parms,...)

}


                                                                                                                                                                                                            