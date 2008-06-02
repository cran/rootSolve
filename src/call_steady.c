/* Patterned on code call_lsoda.c from package odesolve */

#include <time.h>
#include <string.h>
#include "steady.h"   
                           
void F77_NAME(dsteady)(void (*)(int *, double *, double *, double *, double*, int*),
		     int *, int *, double *, double *, double *, double *,
		     int *, int *, int *, 
		     int *, double *,double *,double *, int*, 
		     void (*)(int *, double *, double *, int *,
			            int *, double *, int *, double*, int*),
		     int *, int *, double *, double *, double *, int *, double *, int *, 
         double *, int *);

typedef void deriv_func(int *, double *, double *,double *,double *, int *);
deriv_func *derivb; 


static void steady_derivs (int *neq, double *t, double *y, double *ydot, double *yout, int *iout)
{

  int i;
  SEXP R_fcall, ans;     

  REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(steady_deriv_func,Time,Y)) ;incr_N_Protect();
  PROTECT(ans = eval(R_fcall, steady_envir))         ;incr_N_Protect();

  for (i = 0; i < *neq; i++)	ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
  my_unprotect(2);      

}

/* assumes 1-D multi-species model; rearrange state vars and rates of change */

static void steady_derivs2(int *neq, double *t, double *y, double *ydot, double *yout, int *iout)
{
  int i, j, ii;
     ii = 0;
     for (i = 0; i < ndim; i++)         /* rearrange states */
      { for (j = 0; j < nspec; j++)
        /*y2[j* ndim+i] = y[i* nspec+j];   */
        y2[j* ndim+i] = y[ii++];   
        }
        
     derivb (neq, t, y2, dy2, yout, iout) ;
     ii = 0;
     for (i = 0; i < ndim; i++)         /* rearrange rates of change */
      { for (j = 0; j < nspec; j++)
       /*  ydot[i* nspec+j]=dy2[j* ndim+i];   */
         ydot[ii++]=dy2[j* ndim+i];            
         }
        
     }


static void steady_jac (int *neq, double *t, double *y, int *ml,
		    int *mu, double *pd, int *nrowpd, double *RPAR, int *IPAR)
{
  int i, j;
  SEXP R_fcall, ans;

  REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(steady_jac_func,Time,Y));  incr_N_Protect();
  PROTECT(ans = eval(R_fcall, steady_envir));        incr_N_Protect();

  for (i = 0; i < *neq; i++)
    for (j = 0; j < *nrowpd; j++)
    {
      pd[i * (*nrowpd) + j] = REAL(ans)[i * (*neq) + j];
    }
  my_unprotect(2);
}

typedef void jac_func(int *, double *, double *, int *,
		                  int *, double *, int *, double *, int *);
typedef void init_func(void (*)(int *, double *));

SEXP call_dsteady(SEXP y, SEXP time, SEXP func, SEXP parms, SEXP chtol, 
		SEXP atol, SEXP rtol, SEXP itol, SEXP rho, SEXP jacfunc, SEXP initfunc, 
		SEXP verbose, SEXP mf, SEXP BU, SEXP BD, SEXP nIter, SEXP Pos, 
    SEXP nOut,SEXP nAbd, SEXP nSpec, SEXP nDim, SEXP Rpar, SEXP Ipar)
{
  SEXP   yout, RWORK, IWORK;
  int    j, k, ny, isOut, maxit, isSteady;
  double *svar, *dy, *beta, *alpha, tin, *Atol, *Rtol, Chtol, *out;
  double *copyvar, *delt, *precis, *ewt ;
  int    neq, bu, bd, jt, niter, mflag, nout, ntot, nabd, pos, *indx, Itol;
  int    *ipar, lrpar, lipar, len, isDll, rearrange;
    
  deriv_func *derivs;
  jac_func   *jac;
  init_func  *initializer;

  init_N_Protect();
  jt = INTEGER(mf)[0];        

  bu = INTEGER(BU)[0];        
  bd = INTEGER(BD)[0];  
  nabd = INTEGER(nAbd)[0];      
  ny = LENGTH(y);
  Itol = INTEGER(itol)[0];
  maxit = INTEGER(nIter)[0];  
  pos = INTEGER(Pos)[0];  
  
  neq = ny; 
  mflag = INTEGER(verbose)[0];
  nout  = INTEGER(nOut)[0];
  rearrange = 0;
  if (jt ==0)   /* state variables and rate of changes need rearranging*/
  {
   jt =25;
   rearrange = 1;
  } 

  if (inherits(func, "NativeSymbol"))  /* function is a dll */
  {
   isDll = 1;
   if (nout > 0) isOut = 1; 
   ntot  = neq + nout;          /* length of yout */

   lrpar = nout + LENGTH(Rpar); /* length of rpar; LENGTH(Rpar) is always >0 */
   lipar = 3 + LENGTH(Ipar);    /* length of ipar */
  } else                        /* function is not a dll */
  {
   isDll = 0;
   isOut = 0;
   ntot = neq;
   lipar = 1;
   lrpar = 1; 
  }

   out   = (double *) R_alloc(lrpar, sizeof(double));
   ipar  = (int *)    R_alloc(lipar, sizeof(int));

   if (isDll ==1)
   {
    ipar[0] = nout;
    ipar[1] = lrpar;
    ipar[2] = lipar;
    for (j = 0; j < LENGTH(Ipar);j++) ipar[j+3] = INTEGER(Ipar)[j];
   
    for (j = 0; j < nout; j++) out[j] = 0.;  
    for (j = 0; j < LENGTH(Rpar);j++) out[nout+j] = REAL(Rpar)[j];
   }

  /* initialise global variables... */

  PROTECT(Time = NEW_NUMERIC(1))                   ;incr_N_Protect(); 
  PROTECT(Y = allocVector(REALSXP, neq))           ;incr_N_Protect();        

  /* copies of all variables that will be changed in the FORTRAN subroutine */

  indx = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) indx[j] = 0;
 
  svar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) svar[j] = REAL(y)[j];

  if (jt > 23)
  {
   delt = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) delt[j] = 0; 
  } else delt = (double *) R_alloc(1, sizeof(double));
  
  dy = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) dy[j] = 0; 

  ewt = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) ewt[j] = 0; 

  beta = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) beta[j] = 0; 

  copyvar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) copyvar[j] = 0; 

  alpha = (double *) R_alloc(neq*nabd, sizeof(double));
    for (j = 0; j < neq*nabd; j++) alpha[j] = 0; 
  
  len = LENGTH(atol);  
  Atol = (double *) R_alloc(len, sizeof(double));
    for (j = 0; j < len; j++) Atol[j] = REAL(atol)[j];

  len = LENGTH(rtol);  
  Rtol = (double *) R_alloc(len, sizeof(double));
    for (j = 0; j < len; j++) Rtol[j] = REAL(rtol)[j];

  Chtol = REAL(chtol)[0];

  precis =(double *) R_alloc(maxit, sizeof(double));
    for (j = 0; j < maxit; j++) precis[j] = 0;
  
  PROTECT(yout = allocVector(REALSXP,ntot))    ; incr_N_Protect();

  PROTECT(st_gparms = parms)                   ; incr_N_Protect();  

 /* The initialisation routine */
  if (!isNull(initfunc))
    	{
	     initializer = (init_func *) R_ExternalPtrAddr(initfunc);
	     initializer(Initstparms); 	}

 /* pointers to functions derivs and jac, passed to the FORTRAN subroutine */

  if (isDll==1) 
    { 
    if (rearrange == 0)
      {
      derivs = (deriv_func *) R_ExternalPtrAddr(func);
      } else {
       nspec=INTEGER(nSpec)[0];
       ndim =INTEGER(nDim)[0];
       derivs = (deriv_func *) steady_derivs2; 
       derivb = (deriv_func *) R_ExternalPtrAddr(func);

       y2 = (double *) R_alloc(neq, sizeof(double));
       dy2 = (double *) R_alloc(neq, sizeof(double));   
      }
    } else /* not a DLL */
    {  derivs = (deriv_func *) steady_derivs;  
      PROTECT(steady_deriv_func = func); incr_N_Protect();
      PROTECT(steady_envir = rho);incr_N_Protect();
    } 
    
   if (!isNull(jacfunc))
    {
      if (inherits(jacfunc,"NativeSymbol"))
     	{
	    jac = (jac_func *) R_ExternalPtrAddr(jacfunc);
	    } else {
	    steady_jac_func = jacfunc;
	    jac = steady_jac;
	    }
    }

    tin = REAL(time)[0];

      
	  F77_CALL(dsteady) (derivs, &neq, &nabd, &tin, svar, beta, alpha,
			   &jt, &bu, &bd, &maxit, &Chtol, Atol, Rtol, &Itol, jac, &pos, &isSteady, 
         delt, copyvar, ewt, indx, precis, &niter, out, ipar);

	  for (j = 0; j < ny; j++)
	    REAL(yout)[j] = svar[j];
   
	  if (isOut == 1) 
    {
        derivs (&neq, &tin, svar, dy, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[j + ny] = out[j]; 
    }
 

  PROTECT(RWORK = allocVector(REALSXP, niter));incr_N_Protect();
  for (k = 0;k<niter;k++) REAL(RWORK)[k] = precis[k];

  setAttrib(yout, install("precis"), RWORK);    

  PROTECT(IWORK = allocVector(INTSXP, 1));incr_N_Protect();
                          INTEGER(IWORK)[0] = isSteady;
  
  setAttrib(yout, install("steady"), IWORK);    
       
  unprotect_all();
  return(yout);
}

