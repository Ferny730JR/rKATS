/* This file was automatically generated using:
tools::package_native_routine_registration_skeleton(".")
*/

#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME: 
   Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP count_kmers_R(void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP enrichments_R(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern SEXP ikke_R(void *, void *, void *, void *, void *, void *, void *);
extern SEXP seqseq_R(void *, void *, void *);

static const R_CallMethodDef CallEntries[] = {
    {"count_kmers_R", (DL_FUNC) &count_kmers_R,  9},
    {"enrichments_R", (DL_FUNC) &enrichments_R, 10},
    {"ikke_R",        (DL_FUNC) &ikke_R,         7},
    {"seqseq_R",      (DL_FUNC) &seqseq_R,       3},
    {NULL, NULL, 0}
};

void R_init_rkats(DllInfo *dll)
{
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
