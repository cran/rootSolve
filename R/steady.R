## =============================================================================
## steady, -- solves the steady-state condition of
## ordinary differential equation systems
## has similar calling sequence as integration routines from package deSolve
## =============================================================================

steady  <- function (y, time=NULL, func, parms=NULL, method="stode", ...)  {

  if (!method %in% c("stode", "stodes","runsteady"))
    stop (" 'method' should be one of 'stode', 'stodes', 'runsteady'")   

  if (is.null(time)) {
    if (method %in% c("stode", "stodes")) 
      time <- 0
    else
      time <- c(0,Inf)  
  }  

  if (method=="stode")
    stode(y,time,func,parms=parms,...)  else
  if (method=="stodes")
    stodes(y,time,func,parms=parms,...) else
  if (method=="runsteady")
    runsteady(y,times=time,func,parms=parms,...)

}


                                                                                                                                                                                                            