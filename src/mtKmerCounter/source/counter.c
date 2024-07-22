#include <zlib.h>
#include <stdbool.h>
#include <errno.h>
#include <string.h>

#include "counter.h"
#include "hash_functions.h"
#include "memory_utils.h"

#define BUFFER_SIZE 65536U

/*============ Counting Function Declarations ============*/
static KatssCounter *
count_file(const char *filename, unsigned int kmer, const char filetype);

/*============= Helper Function Declarations =============*/
static char
determine_filetype(const char *filename);

static bool
is_nucleotide(char character);

/*============= Actual Functions Declarations =============*/
KatssCounter *
katss_count_kmers(const char *filename, unsigned int kmer)
{
	char filetype = determine_filetype(filename);
	if(filetype == 'e') { /* Error determining filetype */
		error_message("Unable to read sequence from file.\nCurrent supported file types are:"
		              " FASTA, FASTQ, and file containing sequences per line.");
		return NULL;
	} else if (filetype == 'N') { /* Error opening file */
		return NULL;
	}

	KatssCounter *counter = count_file(filename, kmer, filetype);
	return counter;
}


static KatssCounter *
count_file(const char *filename, unsigned int kmer, const char filetype)
{
	/* Open file and prepare counter & hasher */
	gzFile read_file = gzopen(filename, "r");
	if(read_file == NULL)
		return NULL;

	KatssHasher *hasher = katss_init_hasher(kmer, filetype);
	if(hasher == NULL) {
		gzclose(read_file);
		return NULL;
	}

	KatssCounter *counter = katss_init_counter(kmer);
	if(counter == NULL) {
		gzclose(read_file);
		free(hasher);
		return NULL;
	}

	/* Prepare file reading & hash int */
	char buffer[BUFFER_SIZE+1] = { 0 };
	size_t still_reading;
	uint32_t hash_value;

	do {
		still_reading = gzread(read_file, buffer, BUFFER_SIZE);
		buffer[still_reading] = '\0';

		katss_set_seq(hasher, buffer, filetype);
		while(katss_get_fh(hasher, &hash_value, filetype)) {
			katss_increment(counter, hash_value);
		}
	} while(still_reading == BUFFER_SIZE);
	gzclose(read_file);
	free(hasher);

	return counter;
}

/*==============================================================
|  Helper Functions                                            |
==============================================================*/
static char
determine_filetype(const char *file)
{
	/* Open the gzFile, return 'e' upon error */
	gzFile reads_file = gzopen(file, "r");
	if(reads_file == NULL) {
		error_message("katss: %s: %s", file, strerror(errno));
		gzclose(reads_file);
		return 'N';
	}

	char buffer[BUFFER_SIZE];
	int lines_read = 0;
	int fastq_score_lines = 0;
	int sequence_lines = 0;

	while (gzgets(reads_file, buffer, BUFFER_SIZE) != NULL && lines_read < 10) {
		lines_read++;
		char first_char = buffer[0];

		/* Check if the first line starts with '@' for FASTQ */
		if (first_char == '@' && lines_read % 4 == 1) {
			fastq_score_lines++;

		/* Check if the third line starts with '+' for FASTQ */
		} else if (first_char == '+' && lines_read % 4 == 3) {
			fastq_score_lines++;

		/* Check if the line starts with '>' or ';' for FASTA */
		} else if (first_char == '>' || first_char == ';') {
			gzclose(reads_file);
			return 'a';
		} else {
			// Check for nucleotide characters
			int num_total = 0, num = 0;
			for(int i = 0; buffer[i] != '\0'; i++) {
				if(is_nucleotide(buffer[i])) {
					num++;
				}
				num_total++;
			}
			if((double)num/num_total > 0.9) {
				sequence_lines++;
			}
		}
	}

    gzclose(reads_file);

    if (fastq_score_lines >= 2) {
        return 'q'; // fastq file
    } else if (sequence_lines == 10) {
        return 'r'; // raw sequences file
    } else {
        return 'e'; // unsupported file type
    }
}


static bool
is_nucleotide(char character)
{
	switch(character) {
		case 'A':   return true;
		case 'a':   return true;
		case 'C':   return true;
		case 'c':   return true;
		case 'G':   return true;
		case 'g':   return true;
		case 'T':   return true;
		case 't':   return true;
		case 'U':   return true;
		case 'u':   return true;
		default:    return false;
	}
}
