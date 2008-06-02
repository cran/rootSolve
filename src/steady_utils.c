/* Define some global variables and functions that operate on some of them */
/* Patterned on code odesolve_utils.c from package odesolve */
#include <R.h>
#include <Rdefines.h>

/* keep track of SEXPs PROTECTed, and UNPROTECTing them in case of a stop. */
long int N_Protected;

void init_N_Protect(void) { N_Protected = 0; }

void incr_N_Protect(void) { N_Protected++; }

void unprotect_all(void) { UNPROTECT((int) N_Protected); }

void my_unprotect(int n)
{
    UNPROTECT(n);
    N_Protected -= n;
}

/* Globals :*/ 

SEXP steady_deriv_func;
SEXP steady_jac_func;
SEXP steady_envir;
SEXP st_gparms;

SEXP stsparse_deriv_func;
SEXP stsparse_jac_func;
SEXP stsparse_envir;

/* runsteady */
SEXP lsode_deriv_func;
SEXP lsode_jac_func;
SEXP lsode_envir;

void Initstparms(int *N, double *parms)
{
  int i, Nparms;

  Nparms = LENGTH(st_gparms);
  if ((*N) != Nparms)
    {
      PROBLEM "Confusion over the length of parms"
      ERROR;
    } 
  else
    {
      for (i = 0; i < *N; i++) parms[i] = REAL(st_gparms)[i];
    }
}
  


