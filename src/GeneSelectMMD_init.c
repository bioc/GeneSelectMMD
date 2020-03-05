#include <R_ext/RS.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Fortran calls */
extern void F77_NAME(llkhfun)(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(maxposfun)(void *, void *, void *, void *, void *);
extern void F77_NAME(myttest)(void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(paraestloop)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(wifun)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

static const R_FortranMethodDef FortranEntries[] = {
    {"llkhfun",     (DL_FUNC) &F77_NAME(llkhfun),      9},
    {"maxposfun",   (DL_FUNC) &F77_NAME(maxposfun),    5},
    {"myttest",     (DL_FUNC) &F77_NAME(myttest),      7},
    {"paraestloop", (DL_FUNC) &F77_NAME(paraestloop), 20},
    {"wifun",       (DL_FUNC) &F77_NAME(wifun),       11},
    {NULL, NULL, 0}
};

void R_init_GeneSelectMMD(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
