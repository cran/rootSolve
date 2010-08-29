### ============================================================================
### first some common functions
### ============================================================================

## =============================================================================
## Set the mfrow parameters and whether to "ask" for opening a new device
## =============================================================================

setplotpar <- function(nmdots, dots, nv, ask) {
    if (!any(match(nmdots, c("mfrow", "mfcol"), nomatch = 0))) {
      nc <- min(ceiling(sqrt(nv)),3)
      nr <- min(ceiling(nv/nc),3)
      mfrow <- c(nr, nc)
    }
    else if ("mfcol" %in% nmdots)
        mfrow <- rev(dots$mfcol)
    else mfrow <- dots$mfrow

    if (! is.null(mfrow)) {
      mf <- par(mfrow=mfrow)
    }

   ## interactively wait if there are remaining figures
    if (is.null(ask))
      ask <- prod(par("mfrow")) < length(which) && dev.interactive()

    return(ask)
}

## =============================================================================
## find a variable
## =============================================================================

selectstvar <- function (which,var) {

    if (!is.numeric(which)) {
        ln <- length(which)
        # keep ordering...
        Select <- NULL
        for ( i in 1:ln) {
          ss <- which(which[i]==var)
          if (length(ss) ==0)
            stop("variable", which[i], "not in var")
          Select <- c(Select,ss)
        }        
    }
    else {
        Select <- which   # "Select now refers to the column number
        if (max(Select) > length(var))
            stop("index in 'which' too large")
        if (min(Select) < 1)
            stop("index in 'which' should be > 0")
    }
  return(Select)
}

### ============================================================================
### S3 methods
### ============================================================================

plot.steady1D <- function (x, which = NULL, grid = NULL, xyswap=FALSE, 
  ask = NULL, ...) {

# if x is vector, check if there is more than one species...  
    X <- x$y
    if (is.vector(X)) {
      nspec <- attributes(x)$nspec
      if (length(X)%%nspec != 0) 
        stop("length of 'x' should be a multiple of 'nspec' if x is a vector")
      x <- matrix(nc = nspec, data = X)
    } else x <- X   # only state variables
   
      
    if (is.null(which)) which <- 1:ncol(x)
    var <- colnames(x)
    if(is.null(var)) var <- 1:ncol(x)
    which <- selectstvar(which,var)
    
    if (is.null(grid)) 
       grid <- 1:nrow(x)
    if (length(grid) != nrow(x)) 
      stop("length of grid (x-axis) should be = number of rows in 'x$y'")
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  ""   else dots$ylab 

    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        dots$main <- Main[i]
            
        if (! xyswap) {            
          dots$xlab <- xxlab[i]
          dots$ylab <- yylab[i]
          do.call("plot", c(alist(grid, x[, i]), dots))
        } else {
          if (is.null(labs)){
            dots$ylab <- xxlab[i]
            dots$xlab <- yylab[i]
          } else {
            dots$xlab <- xxlab[i]
            dots$ylab <- yylab[i]
            dots$ylim <- rev(range(grid))    # y-axis reversed
          }
          do.call("plot", c(alist(x[, i], grid), dots))
        } 
    }
}

### ============================================================================
# to drape a color over a persp plot.
### ============================================================================


drapecol <- function (A, col = colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100), NAcol = "white") 
{
    nr <- nrow(A)
    nc <- ncol(A)
    ncol <- length(col)
    AA <- 0.25 * (A[1:(nr - 1), 1:(nc - 1)] + A[1:(nr - 1), 2:nc] + 
        A[2:nr, 1:(nc - 1)] + A[2:nr, 2:nc])
    Ar <- range(AA, na.rm = TRUE)
    rn <- Ar[2] - Ar[1]
    ifelse(rn != 0, drape <- col[1 + trunc((AA - Ar[1])/rn * 
        (ncol - 1))], drape <- rep(col[1], ncol))
    drape[is.na(drape)] <- NAcol
    return(drape)
}


### ============================================================================

image.steady2D <- function (x, which = NULL, 
    add.contour = FALSE, grid = NULL, ask = NULL, method="image", ...) {

# if x is vector, check if there is more than one species...  
    X <- x$y
    out <- list()
    nspec <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    if (is.vector(X)) {
      if (length(X) - nspec*prod(dimens) != 0) 
        stop("length of 'x' should be = 'nspec' * prod(dimens) if x is a vector")
      x <- matrix(nc = nspec, data = X)
      
      for ( i in 1:nspec){
        istart <- (i-1)*prod(dimens) 
        out[[i]] <- matrix(nr=dimens[1], nc=dimens[2], data =
          X[(istart+1):(istart+prod(dimens))])
      }
    } else 
        out <- X   # only state variables
      
    if (is.null(which)) which <- 1:nspec
    var <- 1:nspec
    which <- selectstvar(which,var)
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  "y"   else dots$ylab 

    if (method=="persp")
      dotscol <- dots$col 

    else if (method == "filled.contour")
    dots$color.palette <- if (is.null(dots$color.palette)) 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))else dots$color.palette
    else
    dots$col <- if (is.null(dots$col)) 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col
    
    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
        dots$main <- Main[i]
            
        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]
        List <- alist(z=out[[i]])
        if (! is.null(grid)) {
          List$x <- grid[[1]]
          List$y <- grid[[2]]
        }
        if (method=="persp")
           if(is.null(dotscol))  
             dots$col <- drapecol(out[[i]],
               colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
              "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100))
           else
              dots$col<-drapecol(out[[i]],dotscol)
        
        do.call(method, c(List, dots)) 
        if (add.contour) do.call("contour", c(List, add=TRUE))
          
    }
}

### ============================================================================

image.steady3D <- function (x, which = NULL, 
    add.contour = FALSE, grid = NULL, ask = NULL, method="image", ...) {

# if x is vector, check if there is more than one species...  
    X <- x$y
    out    <- list()
    nspec  <- attributes(x)$nspec
    dimens <- attributes(x)$dimens
    Nz <- dimens[3]
    if (is.vector(X)) {
      if (length(X) - nspec*prod(dimens) != 0) 
        stop("length of 'x' should be = 'nspec' * prod(dimens) if x is a vector")
      x <- matrix(nc = nspec, data = X)
      
      for ( i in 1:nspec){
        istart <- (i-1)*prod(dimens) 
        out[[i]] <- array(dim = dimens, data =
          X[(istart+1):(istart+prod(dimens))])
      }
    } else 
        out <- X   # only state variables
      
    if (is.null(which)) which <- 1:nspec
    var <- 1:nspec
    which <- selectstvar(which,var)
    
    np <- length(which)

    dots <- list(...)
    nmdots <- names(dots)

    # number of figures in a row and 
    # interactively wait if there are remaining figures
   
    ask <- setplotpar(nmdots, dots, np, ask)
    if (ask) {
        oask <- devAskNewPage(TRUE)
        on.exit(devAskNewPage(oask))
    }
    
    Main <-  if (is.null(dots$main)) var else rep(dots$main, length.out =np)

    labs <- (is.null(dots$xlab) && is.null(dots$ylab))
    xxlab <- if (is.null(dots$xlab))  "x"  else dots$xlab
    yylab <- if (is.null(dots$ylab))  "y"   else dots$ylab 

    if (method=="persp")
      dotscol <- dots$col 

    else if (method == "filled.contour")
    dots$color.palette <- if (is.null(dots$color.palette)) 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))else dots$color.palette
    else
    dots$col <- if (is.null(dots$col)) 
      colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
             "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100) else dots$col
    
    ## allow individual xlab and ylab (vectorized)
    xxlab <- rep(xxlab, length.out = np)
    yylab <- rep(yylab, length.out = np)

    for (i in which) {
            
        dots$xlab <- xxlab[i]
        dots$ylab <- yylab[i]
        List <- list()
        if (! is.null(grid)) {
          List$x <- grid[[1]]
          List$y <- grid[[2]]
        }
        for (z in 1:Nz){
          dots$main <- paste("var",Main[i],"z = ", z)
          zdat <- out[[i]][,,z]
          List$z <- zdat
          if (method=="persp")
            if(is.null(dotscol))  
              dots$col <- drapecol(zdat,
                colorRampPalette(c("#00007F", "blue", "#007FFF", "cyan",
               "#7FFF7F", "yellow", "#FF7F00", "red", "#7F0000"))(100))
            else
              dots$col<-drapecol(zdat,dotscol)
        
          do.call(method, c(List, dots)) 
          if (add.contour) do.call("contour", c(List, add=TRUE))
       }   
    }
}

