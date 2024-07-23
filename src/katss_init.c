
#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* This file was automatically generated using:
tools::package_native_routine_registration_skeleton(".")
*/

/* .Call calls */
extern SEXP count_kmers_R(void *, void *);
extern SEXP enrichments_R(void *, void *, void *, void *);
extern SEXP ikke_R(void *, void *, void *, void *, void *, void *);
extern SEXP seqseq_R(void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"count_kmers_R", (DL_FUNC) &count_kmers_R, 2},
    {"enrichments_R", (DL_FUNC) &enrichments_R, 6},
    {"ikke_R",        (DL_FUNC) &ikke_R,        7},
    {"seqseq_R",      (DL_FUNC) &seqseq_R,      3},
    {NULL, NULL, 0}
};

void R_init_rkatss(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
