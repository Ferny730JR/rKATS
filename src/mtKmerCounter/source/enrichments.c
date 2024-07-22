#include <stdbool.h>
#include <float.h>
#include <math.h>

#include "katss_core.h"
#include "enrichments.h"
#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"

static int compare(const void *a, const void *b);

KatssEnrichments *
katss_compute_enrichments(KatssCounter *test, KatssCounter *control, bool normalize)
{
	/* Test and control counts must be the same size */
	if(test->kmer != control->kmer) {
		return NULL;
	}

	/* Allocate enrichments struct */
	uint32_t num_enrichments = test->capacity;
	KatssEnrichments *enrichments = s_malloc(sizeof *enrichments);
	enrichments->enrichments = s_malloc(num_enrichments * sizeof(KatssEnrichment));
	enrichments->num_enrichments = ((uint64_t)num_enrichments)+1;

	/* Compute enrichments */
	for(uint32_t i=0; i<=test->capacity; i++) {
		double test_count, control_count;
		katss_get_from_hash(test, KATSS_DOUBLE, &test_count, i);
		katss_get_from_hash(control, KATSS_DOUBLE, &control_count, i);

		enrichments->enrichments[i].key = i; // Set key
		if(test_count == 0.0 || control_count == 0.0) { // Determine if enrichment is valid
			enrichments->enrichments[i].enrichment = NAN;
			continue;
		}

		double r_val = (test_count/test->total)/(control_count/control->total);
		if(normalize) {
			r_val = log2(r_val);
		}

		enrichments->enrichments[i].enrichment = r_val;
	}
	qsort(enrichments->enrichments, num_enrichments, sizeof(KatssEnrichment), compare);
	return enrichments;
}


KatssEnrichments *
katss_enrichments(const char *test_file, const char *control_file, unsigned int kmer, bool normalize)
{
	KatssCounter *test_counts = katss_count_kmers(test_file, kmer);
	if(test_counts == NULL) {
		return NULL;
	}

	KatssCounter *control_counts = katss_count_kmers(control_file, kmer);
	if(control_counts == NULL) {
		return NULL;
	}

	KatssEnrichments *enrichments = katss_compute_enrichments(test_counts, control_counts, normalize);
	katss_free_counter(control_counts);
	katss_free_counter(test_counts);
	return enrichments;
}


KatssEnrichments *
katss_ikke(const char *test_file, const char *control_file, unsigned int kmer, uint64_t iterations, bool normalize)
{
	/* Get the counts for the test_file */
	KatssCounter *test_counts = katss_count_kmers(test_file, kmer);
	if(test_counts == NULL) {
		return NULL;
	}

	/* Get the counts for the control file */
	KatssCounter *control_counts = katss_count_kmers(control_file, kmer);
	if(control_counts == NULL) {
		return NULL;
	}

	KatssEnrichments *enrichments = s_malloc(sizeof(KatssEnrichments));
	enrichments->enrichments = s_malloc(iterations * sizeof *enrichments->enrichments);
	enrichments->num_enrichments = iterations > test_counts->capacity ? ((uint64_t)test_counts->capacity)+1 : iterations;

	/* Get the first top kmer */
	enrichments->enrichments[0] = katss_top_enrichment(test_counts, control_counts, normalize);

	/* Subsequent iterations begin uncounting */
	for(uint32_t i=1; i<iterations; i++) {
		char kseq[17];
		katss_unhash(kseq, enrichments->enrichments[i-1].key, test_counts->kmer, true);
		katss_uncount_kmer(test_counts, test_file, kseq);
		katss_uncount_kmer(control_counts, control_file, kseq);
		enrichments->enrichments[i] = katss_top_enrichment(test_counts, control_counts, normalize);
	}

	return enrichments;
}


KatssEnrichments *
katss_ikke_mt(const char *test_file, const char *control_file, unsigned int kmer, 
              uint64_t iterations, bool normalize, int threads)
{
	/* Get the counts for the test_file */
	KatssCounter *test_counts = katss_count_kmers(test_file, kmer);
	if(test_counts == NULL) {
		return NULL;
	}

	/* Get the counts for the control file */
	KatssCounter *control_counts = katss_count_kmers(control_file, kmer);
	if(control_counts == NULL) {
		return NULL;
	}

	KatssEnrichments *enrichments = s_malloc(sizeof *enrichments);
	enrichments->enrichments = s_malloc(iterations * sizeof *enrichments->enrichments);
	enrichments->num_enrichments = iterations > test_counts->capacity ? ((uint64_t)test_counts->capacity)+1 : iterations;

	/* Get the first top kmer */
	enrichments->enrichments[0] = katss_top_enrichment(test_counts, control_counts, normalize);

	/* Subsequent iterations begin uncounting */
	for(uint32_t i=1; i<iterations; i++) {
		char kseq[17];
		katss_unhash(kseq, enrichments->enrichments[i-1].key, test_counts->kmer, true);
		katss_uncount_kmer_mt(test_counts, test_file, kseq, threads);
		katss_uncount_kmer_mt(control_counts, control_file, kseq, threads);
		enrichments->enrichments[i] = katss_top_enrichment(test_counts, control_counts, normalize);
	}

	return enrichments;
}


KatssEnrichment
katss_top_enrichment(KatssCounter *test, KatssCounter *control, bool normalize)
{
	KatssEnrichment top_kmer = {.enrichment = -DBL_MAX};
	uint32_t top_enrichment_hash = 0;
	double top_enrichment = -DBL_MAX;

	/* Sanity check, make sure total_count is greater than 0 */
	if(!control->total || !test->total) {
		return top_kmer;
	}

	for(uint32_t i=0; i<control->capacity; i++) {
		/* Get frequencies of input and bound */
		double test_frq, control_frq;
		katss_get_from_hash(test, KATSS_DOUBLE, &test_frq, i);
		katss_get_from_hash(control, KATSS_DOUBLE, &control_frq, i);

		if(test_frq == 0 || control_frq == 0) {
			continue;
		}

		control_frq /= control->total;
		test_frq /= test->total;

		/* get enrichment of current kmer */
		double cur_enrichment = test_frq/control_frq;
		if(normalize) {
			cur_enrichment = log2(cur_enrichment);
		}

		if(cur_enrichment > top_enrichment) {
			top_enrichment = cur_enrichment;
			top_enrichment_hash = i;
		}
	}

	/* No top kmer was found (probs because something terrible happened) return empty struct */
	if(top_enrichment == top_kmer.enrichment) {
		return top_kmer;
	}

	/* Fill struct with info */
	top_kmer.enrichment = top_enrichment;
	top_kmer.key = top_enrichment_hash;

	/* Everything (probably) worked! Hurray, now return. */
	return top_kmer;
}


void
katss_free_enrichments(KatssEnrichments *enrichments)
{
	free(enrichments->enrichments);
	free(enrichments);
}


static int
compare(const void *a, const void *b)
{
	const KatssEnrichment *p1 = (KatssEnrichment *)a;
	const KatssEnrichment *p2 = (KatssEnrichment *)b;

	if(isnan(p1->enrichment) && isnan(p2->enrichment)) {
		return 0;
	} else if(isnan(p1->enrichment)) {
		return 1;
	} else if(isnan(p2->enrichment)) {
		return -1;
	}

	if(p1->enrichment < p2->enrichment) {
		return 1;
	} else if(p1->enrichment > p2->enrichment) {
		return -1;
	} else {
		return 0;
	}
}
