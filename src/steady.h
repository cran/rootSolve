#include <R.h>
#include <Rdefines.h>
/* global variables */
SEXP Time, Y ;
extern SEXP st_gparms;

int ndim, nspec;
double *y2,*dy2;

/* steady_utils.c globals */
extern SEXP steady_deriv_func;
extern SEXP steady_jac_func;
extern SEXP steady_envir;

extern SEXP stsparse_deriv_func;
extern SEXP stsparse_jac_func;
extern SEXP stsparse_envir;

/* runsteady */
extern SEXP lsode_deriv_func;
extern SEXP lsode_jac_func;
extern SEXP lsode_envir;

/* steady_utils.c utilities */
void init_N_Protect(void);
void incr_N_Protect(void);
void unprotect_all(void);
void my_unprotect(int);

void Initstparms(int *, double *);


