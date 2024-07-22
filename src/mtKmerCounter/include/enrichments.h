#ifndef KATSS_ENRICHMENTS_H
#define KATSS_ENRICHMENTS_H

#include <stdint.h>
#include <stdbool.h>
#include "counter.h"

#ifdef __cplusplus
extern "C" {
#endif 

typedef struct KatssEnrichment {
	double enrichment;
	uint32_t key;
} KatssEnrichment;

typedef struct KatssEnrichments {
	KatssEnrichment *enrichments;
	uint64_t num_enrichments;
} KatssEnrichments;

KatssEnrichments *katss_compute_enrichments(KatssCounter *test, KatssCounter *control, bool normalize);
KatssEnrichments *katss_enrichments(const char *test_file, const char *control_file, unsigned int kmer, bool normalize);
KatssEnrichments *katss_ikke_mt(const char *test_file, const char *control_file, unsigned int kmer, 
                                uint64_t iterations, bool normalize, int threads);
KatssEnrichments *katss_ikke(const char *test_file, const char *control_file, unsigned int kmer, uint64_t iterations, bool normalize);
KatssEnrichment katss_top_enrichment(KatssCounter *test, KatssCounter *control, bool normalize);
void katss_free_enrichments(KatssEnrichments *enrichments);

#ifdef __cplusplus
}
#endif

#endif // KATSS_ENRICHMENTS_H
