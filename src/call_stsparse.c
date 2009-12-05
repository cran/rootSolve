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

static void C_stsparse_derivs (int *neq, double *t, double *y, double *ydot, 
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

SEXP call_stsparse(SEXP y, SEXP time, SEXP func, SEXP parms, SEXP chtol, 
		SEXP atol, SEXP rtol, SEXP itol, SEXP rho, SEXP initfunc, 
		SEXP verbose, SEXP mf, SEXP NNZ, SEXP NSP, SEXP NGP, SEXP nIter, SEXP Posit,
    SEXP Pos, SEXP nOut, SEXP Rpar, SEXP Ipar, SEXP Type, SEXP Ian, SEXP Jan)
{
  SEXP   yout, RWORK, IWORK;
  int    j, k, ny, maxit, isSteady;
  double *svar, *dsvar, *beta, *alpha, tin, *Atol, *Rtol, Chtol;
  double *x, *precis, *ewt, *rsp ;
  int    neq, nnz, nsp, ngp, jt, niter, mflag, posit, *pos, ipos, Itol, type;
  int    *R, *C, *IC, *ian, *jan, *igp, *jgp, *isp, *dims;
  int    len, isDll ;
    
  C_deriv_func_type *derivs;
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

  if (inherits(func, "NativeSymbol"))  /* function is a dll */
     isDll = 1;
  else
     isDll = 0;
   if (nout > 0) isOut = 1; 

  /* initialise output ... */
  initOut(isDll, neq, nOut, Rpar, Ipar);

  /* initialise global variables... */
            
  PROTECT(Time = NEW_NUMERIC(1))                   ;incr_N_Protect(); 
  PROTECT(Y = allocVector(REALSXP, neq))           ;incr_N_Protect();        

  /* copies of all variables that will be changed in the FORTRAN subroutine */

  R = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) R[j] = 0;
 
  C = (int *) R_alloc(neq, sizeof(int));
    for (j = 0; j < ny; j++) C[j] = 0;

  dims = (int *) R_alloc(7, sizeof(int));   /* 7 is maximal amount */
    for (j = 0; j < 7; j++) dims[j] = 0;

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
   
  /* 1-D, 2-D, 3-D problem:  */
  if (type == 2)        /* 1=ncomp,2:dim(x), 3: cyclic(x)*/
    for (j = 0; j<3 ; j++) dims[j] = INTEGER(NNZ)[j+1];
  else if (type == 3)   /* 1=ncomp,2-3:dim(x,y), 4-5: cyclic(x,y)*/
    for (j = 0; j<5 ; j++) dims[j] = INTEGER(NNZ)[j+1];
  else if (type == 4)   /* 1=ncomp,2-4:dim(x,y,z), 5-7: cyclic(x,y,z)*/
    for (j = 0; j<7 ; j++) dims[j] = INTEGER(NNZ)[j+1];

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

 /* The initialisation routine */
  initParms(initfunc, parms);

 /* pointers to functions derivs and jac, passed to the FORTRAN subroutine */

  if (isDll)
    {
      derivs = (C_deriv_func_type *) R_ExternalPtrAddr(func);

    } else {  derivs = (C_deriv_func_type *) C_stsparse_derivs;  
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
