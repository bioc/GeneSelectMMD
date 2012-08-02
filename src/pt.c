#include <R.h>
#include <Rmath.h>

//    wrapper for the R function pt()
double F77_SUB(rpt)(double *x, double *n, int *lower_tail, int *log_p) 
      {
      return pt(*x, *n, *lower_tail, *log_p);
      }


