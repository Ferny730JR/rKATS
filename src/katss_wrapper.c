#include <R.h>
#include <Rinternals.h>

#include "katss.h"

#include "memory_utils.h"
#include "counter.h"
#include "hash_functions.h"
#include "seqseq.h"
#include "enrichments.h"

#define ALGO_COUNT 0
#define ALGO_RVALS 1
#define ALGO_IKKES 2

SEXP
katssdata_to_df(KatssData *data, KatssOptions *opts, int algo)
{
	uint64_t capacity = data->num_kmers;

	/* Dataframe that will store data */
	SEXP df = NULL;
	SEXP col_names = NULL;
	SEXP row_names = NULL;

	/* Vectors that will store columns */
	SEXP kmers = NULL;
	SEXP rvals = NULL;
	SEXP stdev = NULL;
	SEXP pvals = NULL;

	/* Create vectors that will store data in R */
	kmers = PROTECT(allocVector(STRSXP, capacity));
	rvals = PROTECT(allocVector(REALSXP, capacity));
	if(opts->bootstrap_iters > 0)
		stdev = PROTECT(allocVector(REALSXP, capacity));
	if(opts->bootstrap_iters > 0 && algo != ALGO_COUNT)
		pvals = PROTECT(allocVector(REALSXP, capacity));

	/* Create pointers to SEXP vectors to write into them */
	double *rvals_p = REAL(rvals);
	double *stdev_p = stdev == NULL ? NULL : REAL(stdev);
	double *pvals_p = pvals == NULL ? NULL : REAL(pvals);

	/* Move data into vectors */
	char *kseq = s_malloc(opts->kmer + 1); // where kmer str will be stored
	for(uint64_t i=0; i < data->num_kmers; i++) {
		/* Set the kmer string in dataframe */
		katss_unhash(kseq, (uint32_t)i, opts->kmer, true);
		SET_STRING_ELT(kmers, i, mkChar(kseq));

		/* Write values */
		if(algo == ALGO_COUNT && opts->bootstrap_iters == 0)
			rvals_p[i] = data->kmers[i].count;
		else
			rvals_p[i] = data->kmers[i].rval;
		if(opts->bootstrap_iters > 0)
			stdev_p[i] = data->kmers[i].stdev;
		if(opts->bootstrap_iters > 0 && algo != ALGO_COUNT)
			pvals_p[i] = data->kmers[i].pval;
	}
	free(kseq);

	/* Free resources allocated to kdata since we are done copying values */
	katss_free_kdata(data);

	/* Create data.frame from values alongside col & row names */
	if(opts->bootstrap_iters > 0 && algo != ALGO_COUNT) {
		df = PROTECT(allocVector(VECSXP, 4));
		col_names = PROTECT(allocVector(STRSXP, 4));
	} else if(opts->bootstrap_iters > 0) {
		df = PROTECT(allocVector(VECSXP, 3));
		col_names = PROTECT(allocVector(STRSXP, 3));
	} else {
		df = PROTECT(allocVector(VECSXP, 2));
		col_names = PROTECT(allocVector(STRSXP, 2));
	}
	row_names = PROTECT(allocVector(INTSXP, capacity));

	/* Set column and row names */
	SET_STRING_ELT(col_names, 0, mkChar("kmer"));
	if(algo == ALGO_COUNT)
		SET_STRING_ELT(col_names, 1, mkChar("count"));
	else 
		SET_STRING_ELT(col_names, 1, mkChar("rval"));
	if(opts->bootstrap_iters > 0 && algo != ALGO_COUNT)
		SET_STRING_ELT(col_names, 3, mkChar("pval"));
	if(opts->bootstrap_iters > 0)
		SET_STRING_ELT(col_names, 2, mkChar("stdev"));
	for(uint64_t i=0; i < capacity; i++) {
		INTEGER(row_names)[i] = i+1;
	}

	/* Set columns of dataframe from computes values */
	SET_VECTOR_ELT(df, 0, kmers);
	SET_VECTOR_ELT(df, 1, rvals);
	if(opts->bootstrap_iters > 0)
		SET_VECTOR_ELT(df, 2, stdev);
	if(opts->bootstrap_iters > 0 && algo != ALGO_COUNT)
		SET_VECTOR_ELT(df, 3, pvals);
	
	/* Set attributes to ensure R recognizes names and data frame */
	setAttrib(df, R_NamesSymbol, col_names);
	setAttrib(df, R_RowNamesSymbol, row_names);
	setAttrib(df, R_ClassSymbol, mkString("data.frame"));

	/* Release protected vectors */
	if(opts->bootstrap_iters > 0 && algo != ALGO_COUNT)
		UNPROTECT(7);
	else if(opts->bootstrap_iters > 0)
		UNPROTECT(6);
	else
		UNPROTECT(5);
	
	/* Return data.frame */
	return df;
}
// Function to convert R inputs to C and call count_kmers
SEXP
count_kmers_R(SEXP filename, SEXP kmer, SEXP klet, SEXP sort, SEXP iters, 
              SEXP sample, SEXP algo, SEXP seed, SEXP threads)
{
	const char *c_filename = CHAR(STRING_ELT(filename, 0));

	/* Initialize the options */
	KatssOptions opts;
	katss_init_options(&opts);

	/* Modify the options based on user input */
	opts.kmer = INTEGER(kmer)[0];
	opts.probs_ntprec = INTEGER(klet)[0];
	opts.sort_enrichments = INTEGER(sort)[0];
	opts.bootstrap_iters = INTEGER(iters)[0];
	opts.bootstrap_sample = INTEGER(sample)[0];
	opts.seed = INTEGER(seed)[0];
	opts.threads = INTEGER(threads)[0];
	opts.enable_warnings = true;
	if(INTEGER(algo)[0] == 1) {
		opts.probs_algo = KATSS_PROBS_NONE;
	} else {
		opts.probs_algo = KATSS_PROBS_USHUFFLE;
	}

	/* Compute the results */
	KatssData *result = katss_count(c_filename, &opts);
	if(result == NULL)
		return R_NilValue;

	/* Turn result into an R data.frame and return it */
	return katssdata_to_df(result, &opts, ALGO_COUNT);
	// return R_NilValue;
}


/* Function to convert R inputs to C and call katss_enrichments */
SEXP
enrichments_R(SEXP test_file, SEXP ctrl_file, SEXP kmer, SEXP probabilistic, 
              SEXP normalize, SEXP verbose)
{
	const char *test_filename = CHAR(STRING_ELT(test_file, 0));
	unsigned int c_kmer = INTEGER(kmer)[0];
	bool c_normalize = asLogical(normalize) == TRUE;
	bool c_verbose = asLogical(verbose);

	/* Compute the actual enrichments */
	KatssEnrichments *result;
	if(c_verbose)
		Rprintf("Begin calculating enrichments...\n");

	if(asLogical(probabilistic) == TRUE)
		result = katss_prob_enrichments(test_filename, c_kmer, c_normalize);
	else {
		const char *ctrl_filename = CHAR(STRING_ELT(ctrl_file, 0));
		result = katss_enrichments(test_filename, ctrl_filename, c_kmer, c_normalize);
	}

	/* If enrichments failed for whatever reason, return NULL */
	if(result == NULL) {
		Rprintf("Begin calculating enrichments... Failure\n");
		return R_NilValue;
	}
	if(c_verbose)
		Rprintf("Begin calculating enrichments... Success\n"
		        "Creating data.frame...\n");

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

	if(c_verbose)
		Rprintf("Creating data.frame... Success\n");
	return df;
}


/* ikkr_R: C wrapper to perform iterative k-mer knockout enrichments in R */
SEXP
ikke_R(SEXP test_file, SEXP ctrl_file, SEXP kmer, SEXP iterations, SEXP probabilistic,
       SEXP normalize, SEXP threads)
{
	unsigned int c_kmer = INTEGER(kmer)[0];
	uint64_t c_iterations = REAL(iterations)[0];
	bool c_normalize = asLogical(normalize) == TRUE;
	int c_threads = INTEGER(threads)[0];

	KatssEnrichments *result;
	const char *test_filename = CHAR(STRING_ELT(test_file, 0));
	if(asLogical(probabilistic) == TRUE) {
		/* If single threaded, use regular probalistic ikke */
		if(c_threads < 2)
			result = katss_prob_ikke(test_filename, c_kmer, c_iterations, c_normalize);
		
		/* If multithreaded, use mt version */
		else
			result = katss_prob_ikke_mt(test_filename, c_kmer, c_iterations, c_normalize, c_threads);
	} else {
		const char *ctrl_filename = CHAR(STRING_ELT(ctrl_file, 0));
		/* If single threaded, use regular ikke*/
		if(c_threads < 2)
			result = katss_ikke_(test_filename, ctrl_filename, c_kmer, c_iterations, c_normalize);
		
		/* Else if multithreaded, use mt version */
		else
			result = katss_ikke_mt(test_filename, ctrl_filename, c_kmer, c_iterations, c_normalize, c_threads);
	}
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
