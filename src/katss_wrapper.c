#include <R.h>
#include <Rinternals.h>

#include "mtKmerCounter/include/counter.h"
#include "mtKmerCounter/include/hash_functions.h"
#include "mtKmerCounter/include/seqseq.h"
#include "mtKmerCounter/include/enrichments.h"
#include "mtKmerCounter/source/katss_core.h"
#include "mtKmerCounter/source/helpers/memory_utils.h"

// Function to convert R inputs to C and call count_kmers
SEXP count_kmers_R(SEXP filename, SEXP kmer) {
	const char *c_filename = CHAR(STRING_ELT(filename, 0));
	unsigned int c_kmer = INTEGER(kmer)[0];

	KatssCounter *result = katss_count_kmers(c_filename, c_kmer);
	if(result == NULL) {
		return R_NilValue;
	}
	uint64_t capacity = (uint64_t)result->capacity+1;

	SEXP counts_r = PROTECT(allocVector(REALSXP, capacity));
	SEXP kmer_strings = PROTECT(allocVector(STRSXP, capacity));
	double *counts_r_ptr = REAL(counts_r);

	// copy the contents of counts into the R vector
	for(uint32_t i=0; i < capacity; i++) {
		/* Set k-mer counts */
		if(c_kmer<=12) {
			counts_r_ptr[i] = (double)result->table.small[i];
		} else {
			counts_r_ptr[i] = (double)result->table.medium[i];
		}
		/* Set k-mer string */
		char kseq[c_kmer + 1U];
		katss_unhash(kseq, i, c_kmer, 1);
		SET_STRING_ELT(kmer_strings, i, mkChar(kseq));
	}

	/* Free KatssCounter */
	katss_free_counter(result);

	/* Create data.frame from k-mer strings and counts */
	SEXP df = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(df, 0, kmer_strings);
	SET_VECTOR_ELT(df, 1, counts_r);

	/* Set column names for the data.frame */
	SEXP col_names = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(col_names, 0, mkChar("kmer"));
	SET_STRING_ELT(col_names, 1, mkChar("count"));
	setAttrib(df, R_NamesSymbol, col_names);

	/* Set class attribute to "data.frame" */
	SEXP row_names = PROTECT(allocVector(INTSXP, capacity));
	for(uint32_t i=0; i < capacity; i++) {
		INTEGER(row_names)[i] = i+1;
	}
	setAttrib(df, R_RowNamesSymbol, row_names);
	setAttrib(df, R_ClassSymbol, mkString("data.frame"));

	UNPROTECT(5);
	return df;
}


/* Function to convert R inputs to C and call katss_enrichments */
SEXP enrichments_R(SEXP test_file, SEXP ctrl_file, SEXP kmer, Rboolean normalize) {
	const char *test_filename = CHAR(STRING_ELT(test_file, 0));
	const char *ctrl_filename = CHAR(STRING_ELT(ctrl_file, 0));
	unsigned int c_kmer = INTEGER(kmer)[0];
	bool c_normalize = normalize == TRUE;

	KatssEnrichments *result = katss_enrichments(test_filename, ctrl_filename, c_kmer, c_normalize);
	if(result == NULL) {
		return R_NilValue;
	}
	uint64_t capacity = (uint64_t)result->num_enrichments;

	SEXP counts_r = PROTECT(allocVector(REALSXP, capacity));
	SEXP kmer_strings = PROTECT(allocVector(STRSXP, capacity));
	double *counts_r_ptr = REAL(counts_r);

	// copy the contents of counts into the R vector
	for(uint32_t i=0; i < capacity; i++) {
		/* Set k-mer enrichment */
		counts_r_ptr[i] = result->enrichments[i].enrichment;

		/* Set k-mer string */
		char kseq[c_kmer + 1U];
		katss_unhash(kseq, result->enrichments[i].key, c_kmer, 1);
		SET_STRING_ELT(kmer_strings, i, mkChar(kseq));
	}

	/* Free enrichments */
	katss_free_enrichments(result);

	/* Create data.frame from k-mer strings and counts */
	SEXP df = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(df, 0, kmer_strings);
	SET_VECTOR_ELT(df, 1, counts_r);

	/* Set column names for the data.frame */
	SEXP col_names = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(col_names, 0, mkChar("kmer"));
	SET_STRING_ELT(col_names, 1, mkChar("rval"));
	setAttrib(df, R_NamesSymbol, col_names);

	/* Set class attribute to "data.frame" */
	SEXP row_names = PROTECT(allocVector(INTSXP, capacity));
	for(uint32_t i=0; i < capacity; i++) {
		INTEGER(row_names)[i] = i+1;
	}
	setAttrib(df, R_RowNamesSymbol, row_names);
	setAttrib(df, R_ClassSymbol, mkString("data.frame"));

	UNPROTECT(5);
	return df;
}


/* ikkr_R: C wrapper to perform iterative k-mer knockout enrichments in R */
SEXP ikke_R(SEXP test_file, SEXP ctrl_file, SEXP kmer, SEXP iterations, SEXP normalize, SEXP threads) {
	const char *test_filename = CHAR(STRING_ELT(test_file, 0));
	const char *ctrl_filename = CHAR(STRING_ELT(ctrl_file, 0));
	unsigned int c_kmer = INTEGER(kmer)[0];
	uint64_t c_iterations = REAL(iterations)[0];
	bool c_normalize = asLogical(normalize) == TRUE;
	int c_threads = INTEGER(threads)[0];

	KatssEnrichments *result;
	/* If single threaded, use regular ikke*/
	if(c_threads < 2)
		result = katss_ikke(test_filename, ctrl_filename, c_kmer, c_iterations, c_normalize);
	
	/* Else if multithreaded, use mt version */
	else
		result = katss_ikke_mt(test_filename, ctrl_filename, c_kmer, c_iterations, c_normalize, c_threads);

	/* If IKKE failed, return NULL */
	if(result == NULL) {
		return R_NilValue;
	}

	/* Prepare variable to begin copying to R */
	uint64_t capacity = (uint64_t)result->num_enrichments;
	SEXP enrichments_r = PROTECT(allocVector(REALSXP, capacity));
	SEXP kmer_strings = PROTECT(allocVector(STRSXP, capacity));
	double *enrichments_r_ptr = REAL(enrichments_r);

	/* copy the enrichments and k-mer sequences into the R vectors */
	for(uint32_t i=0; i < capacity; i++) {
		/* Set k-mer enrichment */
		enrichments_r_ptr[i] = result->enrichments[i].enrichment;

		/* Set k-mer string */
		char kseq[c_kmer + 1U];
		katss_unhash(kseq, result->enrichments[i].key, c_kmer, 1);
		SET_STRING_ELT(kmer_strings, i, mkChar(kseq));
	}

	/* Free enrichments */
	katss_free_enrichments(result);

	/* Create data.frame from k-mer strings and counts */
	SEXP df = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(df, 0, kmer_strings);
	SET_VECTOR_ELT(df, 1, enrichments_r);

	/* Set column names for the data.frame */
	SEXP col_names = PROTECT(allocVector(STRSXP, 2));
	SET_STRING_ELT(col_names, 0, mkChar("kmer"));
	SET_STRING_ELT(col_names, 1, mkChar("rval"));
	setAttrib(df, R_NamesSymbol, col_names);

	/* Set class attribute to "data.frame" */
	SEXP row_names = PROTECT(allocVector(INTSXP, capacity));
	for(uint32_t i=0; i < capacity; i++) {
		INTEGER(row_names)[i] = i+1;
	}
	setAttrib(df, R_RowNamesSymbol, row_names);
	setAttrib(df, R_ClassSymbol, mkString("data.frame"));

	UNPROTECT(5);
	return df;
}


static inline double
seqseq_single_R(const char *seq, const char *search)
{
	const char *result = seqseq(seq, search);
	if(result == NULL)
		return (double)(0);
	return (double)(result - seq + 1);
}

static inline SEXP
seqseq_multi_R(const char *seq, const char *search)
{
	size_t search_len = strlen(search);
	size_t indcs_len = 1024;
	size_t num_found = 0;
	double *indeces = s_malloc(indcs_len * sizeof(double));

	/* Begin searching for all occurrences of search seq */
	register char *found = (char *)seq;
	while((found = seqseq(found, search)) != NULL) {
		indeces[num_found++] = (double)(found - seq + 1);
		found += search_len;
		if(num_found >= indcs_len) {
			indcs_len *= 2;
			indeces = s_realloc(indeces, indcs_len * sizeof(double));
		}
	}

	/* If none were found, return 0 */
	if(num_found == 0)
		return ScalarReal((double)(0));

	/* Store all matches found in vector */
	SEXP all_matches = PROTECT(allocVector(REALSXP, num_found));
	for(size_t i = 0; i < num_found; i++)
		REAL(all_matches)[i] = indeces[i];

	/* Free resources and return matches */
	free(indeces);
	UNPROTECT(1);
	return all_matches;
}

/* Search for a sequence within a sequence */
SEXP seqseq_R(SEXP seq, SEXP search, SEXP all_matches) {
	const char *seq_ptr = CHAR(STRING_ELT(seq, 0));
	const char *search_ptr = CHAR(STRING_ELT(search, 0));

	/* Find all matches of search in seq and return vector */
	if(asLogical(all_matches) == TRUE) {
		return seqseq_multi_R(seq_ptr, search_ptr);

	/* Find the first occurrence, e.g if search is in the seq */
	} else {
		return ScalarReal(seqseq_single_R(seq_ptr, search_ptr));
	}
}
