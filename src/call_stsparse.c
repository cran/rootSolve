/* SPARSE steady-state solver */

#include <time.h>
#include <string.h>
#include "steady.h"   
                           
void F77_NAME(dsparse)(void (*)(int *, double *, double *, double *, double*, int*),
		     int *, int *, int *, double *, double *, double *, double *, double *,
		     double *, double *, double *, int*, int*, int*, int*, int*, 
         int*, int*, int*, int*, 
         int *, double *, double *, double *, int *,
		     int *, int *, int *, int *, double *, int*,
         int *, double *, int *, int*);

static void stsparse_derivs (int *neq, double *t, double *y, double *ydot, 
                            double *yout, int *iout)
{
  int i;
  SEXP R_fcall, ans;     

  REAL(Time)[0] = *t;
  for (i = 0; i < *neq; i++)  REAL(Y)[i] = y[i];

  PROTECT(R_fcall = lang3(stsparse_deriv_func,Time,Y)) ;incr_N_Protect();
  PROTECT(ans = eval(R_fcall, stsparse_envir))         ;incr_N_Protect();

  for (i = 0; i < *neq; i++)	ydot[i] = REAL(VECTOR_ELT(ans,0))[i];
  my_unprotect(2);      

}



typedef void deriv_func(int *, double *, double *,double *,double *, int *);
typedef void init_func(void (*)(int *, double *));

SEXP call_stsparse(SEXP y, SEXP time, SEXP func, SEXP parms, SEXP chtol, 
		SEXP atol, SEXP rtol, SEXP itol, SEXP rho, SEXP initfunc, 
		SEXP verbose, SEXP mf, SEXP NNZ, SEXP NSP, SEXP NGP, SEXP nIter, SEXP Posit,
    SEXP Pos, SEXP nOut, SEXP Rpar, SEXP Ipar, SEXP Type, SEXP Ian, SEXP Jan)
{
  SEXP   yout, RWORK, IWORK;
  int    j, k, ny, isOut, maxit, isSteady;
  double *svar, *dsvar, *beta, *alpha, tin, *Atol, *Rtol, Chtol, *out;
  double *x, *precis, *ewt, *rsp ;
  int    neq, nnz, nsp, ngp, jt, niter, mflag, nout, ntot, posit, *pos, ipos, Itol, type;
  int    *R, *C, *IC, *ian, *jan, *igp, *jgp, *isp, *dims;
  int    *ipar, lrpar, lipar, len, isDll ;
    
  deriv_func *derivs;
  init_func  *initializer;

  init_N_Protect();

  jt    = INTEGER(mf)[0];        
  nnz   = INTEGER(NNZ)[0];        
  nsp   = INTEGER(NSP)[0];  
  ngp   = INTEGER(NGP)[0];  
  ny    = LENGTH(y);
  Itol  = INTEGER(itol)[0];
  maxit = INTEGER(nIter)[0];  
  type  = INTEGER(Type)[0];

  posit = INTEGER(Posit)[0];   /* positivity of state variables: either specified at once, or via a vector..*/
  ipos = LENGTH(Pos);
  pos = (int *) R_alloc(ipos, sizeof(int));
    for (j = 0; j < ipos; j++) pos[j] = INTEGER(Pos)[j];


  neq   = ny; 
  mflag = INTEGER(verbose)[0];
  nout  = INTEGER(nOut)[0];

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

  R = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) R[j] = 0;
 
  C = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) C[j] = 0;

  dims = (int *) R_alloc(3, sizeof(int));
    for (j = 0; j < 3; j++) dims[j] = 0;

  IC = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) IC[j] = 0;

  svar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) svar[j] = REAL(y)[j];

  dsvar = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) dsvar[j] = 0; 

  beta = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) beta[j] = 0; 

  x = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) x[j] = 0; 

  alpha = (double *) R_alloc(nnz, sizeof(double));
    for (j = 0; j < nnz; j++) alpha[j] = 0; 

  ewt = (double *) R_alloc(neq, sizeof(double));
    for (j = 0; j < ny; j++) ewt[j] = 0; 

  rsp = (double *) R_alloc(nsp, sizeof(double));
    for (j = 0; j < nsp; j++) rsp[j] = 0.;

  ian = (int *) R_alloc(neq+1, sizeof(int));
   if (type == 0) {for (j = 0; j < neq; j++) ian[j] = INTEGER(Ian)[j];} 
   else {for (j = 0; j < neq; j++) ian[j] = 0;}

  jan = (int *) R_alloc(nnz, sizeof(int));
   if (type == 0) 
   {for (j = 0; j < nnz; j++) jan[j] = INTEGER(Jan)[j];} 
   else {for (j = 0; j < nnz; j++) jan[j] = 0;}
   
  /* 1-D or 2-D problem */
  if (type == 2)          {
    dims[0] = INTEGER(NNZ)[1]; /* number components*/ 
    dims[1] = INTEGER(NNZ)[2]; /* dimension x*/ 
  } else if (type == 3)   {
    dims[0] = INTEGER(NNZ)[1]; /* number components*/ 
    dims[1] = INTEGER(NNZ)[2]; /* dimension x*/ 
    dims[2] = INTEGER(NNZ)[3]; /* dimension y*/     
  }


  igp = (int *) R_alloc(ngp+1, sizeof(int));
    for (j = 0; j < ngp+1; j++) igp[j] = 0;
  
  jgp = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < neq; j++) jgp[j] = 0;

  isp = (int *) R_alloc(2*nsp, sizeof(int));
    for (j = 0; j < 2*nsp; j++) isp[j] = 0;
    
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

  if (inherits(func, "NativeSymbol")) 
    {
      derivs = (deriv_func *) R_ExternalPtrAddr(func);

    } else {  derivs = (deriv_func *) stsparse_derivs;  
      PROTECT(stsparse_deriv_func = func); incr_N_Protect();
      PROTECT(stsparse_envir = rho);incr_N_Protect();
    }
    

    
    tin = REAL(time)[0];
      
	  F77_CALL(dsparse) (derivs, &neq, &nnz, &nsp, &tin, svar, dsvar, beta, x,
         alpha, ewt, rsp, ian, jan, igp, jgp, &ngp, R, C, IC, isp,
			   &maxit,  &Chtol, Atol, Rtol, &Itol, &posit, pos, &ipos, &isSteady,
         precis, &niter, dims, out, ipar, &type);

	  for (j = 0; j < ny; j++)
	    REAL(yout)[j] = svar[j];
   
	  if (isOut == 1) 
    {
        derivs (&neq, &tin, svar, dsvar, out, ipar) ;
	      for (j = 0; j < nout; j++)
	       REAL(yout)[j + ny] = out[j]; 
    }
 
  PROTECT(RWORK = allocVector(REALSXP, niter));incr_N_Protect();
  for (k = 0;k<niter;k++) REAL(RWORK)[k] = precis[k];

  setAttrib(yout, install("precis"), RWORK);    

  PROTECT(IWORK = allocVector(INTSXP, 4));incr_N_Protect();
                          INTEGER(IWORK)[0]   = isSteady;
    for (k = 0; k<3; k++) INTEGER(IWORK)[k+1] = dims[k];
  
  setAttrib(yout, install("steady"), IWORK);    
       
  unprotect_all();
  return(yout);
}
