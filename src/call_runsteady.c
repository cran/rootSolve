#include <time.h>
#include <string.h>
#include "steady.h"

/* definition of the call  to the fortran functions - in file lsode.f*/

void F77_NAME(dlsode)(void (*)(int *, double *, double *, double *, double *, int *),
		     int *, double *, double *, double *,
		     int *, double *, double *, int *, int *,
		     int *, double *,int *,int *, int *,
		     void (*)(int *, double *, double *, int *,
			      int *, double *, int *, double *, int *),
		     int *, double *, int *);

/* interface between fortran function call and R function 
   Fortran code calls lsode_derivs(N, t, y, ydot, yout, iout) 
   R code called as lsode_deriv_func(time, y) and returns ydot 
   Note: passing of parameter values and "..." is done in R-function runsteady*/

static void lsode_derivs (int *neq, double *t, double *y, 
                          double *ydot, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                              REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(lsode_deriv_func,Time,Y));   incr_N_Protect();
  PROTECT(ans = eval(R_fcall, lsode_envir));           incr_N_Protect();

  for (i = 0; i < *neq; i++)   ydot[i] = REAL(VECTOR_ELT(ans,0))[i];

  my_unprotect(2);
}

/* interface between fortran call to jacobian and R function */

static void lsode_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;

                             REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++) REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(lsode_jac_func,Time,Y));    incr_N_Protect();
  PROTECT(ans = eval(R_fcall, lsode_envir));          incr_N_Protect();

  for (i = 0; i < *neq * *nrowpd; i++)  pd[i] = REAL(ans)[i];

  my_unprotect(2);
}


/* give name to data types */
typedef void deriv_func(int *, double *, double *,double *, double *, int *);
typedef void jac_func  (int *, double *, double *, int *,
		                    int *, double *, int *, double *, int *);
typedef void init_func (void (*)(int *, double *));

/* MAIN C-FUNCTION, CALLED FROM R-code */

SEXP call_lsode(SEXP y, SEXP times, SEXP func, SEXP parms, SEXP stol, SEXP rtol,
		SEXP atol, SEXP rho, SEXP tcrit, SEXP jacfunc, SEXP initfunc,
		SEXP verbose, SEXP iTask, SEXP rWork, SEXP iWork, SEXP jT, 
    SEXP nOut, SEXP lRw, SEXP lIw, SEXP Rpar, SEXP Ipar)

{
/******************************************************************************/
/******                   DECLARATION SECTION                            ******/
/******************************************************************************/

/* These R-structures will be allocated and returned to R*/
  SEXP yout, RWORK, ISTATE;    

  int  i, j, k, latol, lrtol, lrw, liw, maxit;
  double *xytmp, *rwork, tin, tout, *Atol, *Rtol, Stol, *out, *dy, ss, sumder;
  int neq, itol, itask, istate, iopt, *iwork, jt, mflag, nout, ntot, is;
  int *ipar, lrpar, lipar, isDll, isOut;
  
  deriv_func *derivs;
  jac_func   *jac;
  init_func  *initializer;
    
/******************************************************************************/
/******                         STATEMENTS                               ******/
/******************************************************************************/

/*                      #### initialisation ####                              */    
  init_N_Protect();


  jt  = INTEGER(jT)[0];         /* method flag */
  neq = LENGTH(y);              /* number of equations */ 
  
  mflag = INTEGER(verbose)[0];
  nout   = INTEGER(nOut)[0];    /* number of output variables */
  tin  = REAL(times)[0];        /* start and end time */
  tout = REAL(times)[1];
  Stol = REAL(stol)[0];         /* steady-state tolerance */ 

/* The output:
    out and ipar are used to pass output variables (number set by nout)
    followed by other input (e.g. forcing functions) provided 
    by R-arguments rpar, ipar
    ipar[0]: number of output variables, ipar[1]: length of rpar, 
    ipar[2]: length of ipar */
  
  if (inherits(func, "NativeSymbol"))  /* function is a dll */
  {
   isDll = 1;
   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar);    /* length of ipar */
   isOut = 1;
   ntot  = neq + nout;          /* length of yout */

  } else                              /* function is not a dll */
  {
   isDll = 0;
   lipar = 1;
   lrpar = 1; 
   isOut = 0;   
   ntot  = neq ;        

  }
 
   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;              /* first 3 elements of ipar are special */
    ipar[1] = lrpar;
    ipar[2] = lipar;
    /* other elements of ipar are set in R-function lsodx via argument *ipar* */
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];

    /* first nout elements of rpar reserved for output variables 
      other elements are set in R-function lsodx via argument *rpar* */
    for (j = 0; j < nout; j++)        out[j] = 0.;  
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }
   
/* copies of all variables that will be changed in the FORTRAN subroutine */

  xytmp = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < neq; j++) xytmp[j] = REAL(y)[j];
 
  latol = LENGTH(atol);
  Atol = (double *) R_alloc((int) latol, sizeof(double));

  lrtol = LENGTH(rtol);
  Rtol = (double *) R_alloc((int) lrtol, sizeof(double));

  liw = INTEGER (lIw)[0];
  iwork = (int *) R_alloc(liw, sizeof(int));
     for (j=0; j<LENGTH(iWork); j++) iwork[j] = INTEGER(iWork)[j];
  maxit = iwork[5];
  
  lrw = INTEGER(lRw)[0];
  rwork = (double *) R_alloc(lrw, sizeof(double));
     for (j=0; j<length(rWork); j++) rwork[j] = REAL(rWork)[j];

/* initialise global R-variables... */

  PROTECT(Time = NEW_NUMERIC(1));                  incr_N_Protect();
  PROTECT(Y = allocVector(REALSXP,(neq)));         incr_N_Protect();
  PROTECT(yout = allocVector(REALSXP,ntot));       incr_N_Protect();
  PROTECT(st_gparms = parms);                      incr_N_Protect();

 /* The initialisation routine */
  if (!isNull(initfunc))
	  {
	  initializer = (init_func *) R_ExternalPtrAddr(initfunc);
	  initializer(Initstparms);
	  }

/* pointers to functions derivs, jac, passed to FORTRAN */
  dy = (double *) R_alloc(neq, sizeof(double));
  for (j = 0; j < neq; j++) dy[j] = 0.; 

  if (isDll ==1) 
    { /* DLL address passed to fortran */
      derivs = (deriv_func *) R_ExternalPtrAddr(func);  
      /* no need to communicate with R - but output variables set here */
	  
    } else {
      /* interface function between fortran and R passed to Fortran*/ 
      derivs = (deriv_func *) lsode_derivs; 
      /* needed to communicate with R */
      lsode_deriv_func = func;
      lsode_envir = rho;

    }

  if (!isNull(jacfunc)) 
    {
      if (isDll ==1)
	    {
	     jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	    } else  {
	     lsode_jac_func = jacfunc;
	     jac = lsode_jac;
	    }
    }


/* tolerance specifications */
  if (latol == 1 && lrtol == 1 ) itol = 1;
  if (latol  > 1 && lrtol == 1 ) itol = 2;
  if (latol == 1 && lrtol  > 1 ) itol = 3;
  if (latol  > 1 && lrtol  > 1 ) itol = 4;

  itask = INTEGER(iTask)[0];   
  istate = 1;

  iopt = 0;
  ss = 0.;
  is = 0 ;
  for (i = 5; i < 8 ; i++) ss = ss+rwork[i];
  for (i = 5; i < 10; i++) is = is+iwork[i];
  if (ss >0 || is > 0) iopt = 1; /* non-standard input */

/*                     ####   main time loop   ####                           */    
  for (i = 0; i < maxit; i++)
	{  /* one step */
    F77_CALL(dlsode) (derivs, &neq, xytmp, &tin, &tout,
			   &itol, NUMERIC_POINTER(rtol), NUMERIC_POINTER(atol), 
         &itask, &istate, &iopt, rwork,
			   &lrw, iwork, &liw, jac, &jt, out, ipar); 
     /* check steady-state */
    sumder = 0. ;
    derivs (&neq, &tin, xytmp, dy, out, ipar) ;
    for (j = 0; j < neq; j++) sumder = sumder+fabs(dy[j]); 

    if (sumder < Stol)  break; 
    
    /* errors? */
	  if (istate == -2)
	    {
	      for (j = 0; j < lrtol; j++) Rtol[j] *= 10.0;
	      for (j = 0; j < latol; j++) Atol[j] *= 10.0;
	      warning("Excessive precision requested.  `rtol' and `atol' have been scaled upwards by the factor %g\n",10.0);
	      istate = 3;
	    }
	  if (istate == -1) 
     {
      warning("an excessive amount of work (> maxsteps ) was done, but integration was successful - increase maxsteps");
     }
    if (istate == -3)
	   {
	    unprotect_all();
	    error("Illegal input to lsode\n");
	   }
	} 

  /* If here: hopefully steady-state is reached */  
 
	  for (j = 0; j < neq; j++)
	    REAL(yout)[j] = xytmp[j];
    if (isOut == 1) 
      {
      for (j = 0; j < nout; j++)  REAL(yout)[j + neq + 1] = out[j]; 
	    }
/*                    ####  an error occurred   ####                          */    
  if (istate < 0 ) warning("Returning early.  Results are accurate, as far as they go\n");

  PROTECT(ISTATE = allocVector(INTSXP, 22));incr_N_Protect();
  for (k = 0;k<22;k++) INTEGER(ISTATE)[k+1] = iwork[k];
  INTEGER(ISTATE)[0] = istate;
       
  PROTECT(RWORK = allocVector(REALSXP, 6));incr_N_Protect();
  for (k = 0;k<5;k++) REAL(RWORK)[k] = rwork[k+10];

  REAL(RWORK)[5] = sumder;

  setAttrib(yout, install("rstate"), RWORK);    
  setAttrib(yout, install("istate"), ISTATE);    

/*                       ####   termination   ####                            */    
  unprotect_all();
  return(yout);
}

