#' K-mer counting
#'
#' Count the k-mers in a fastq, fasta, or raw sequences file.
#'
#' @param filename Name of the file which you want to count k-mers from
#' The file has to be of either: raw sequences, fasta, or fastq format. Other
#' file types are currently unsupported and will not work properly if used.
#' @param kmer Length of the k-mer you want to count. Currently, only k-mers up
#' to length 16 are supported.
#'
#' @return Dataframe containing the counts for all k-mers
#'
#' @export
#' @useDynLib rkatss, .registration = TRUE
#'
#' @examples
#' test1 <- do.call(paste0, replicate(50, sample(c("A","C","G","T"), 1000, TRUE), FALSE))
#' tf <- tempfile()
#'
#' # Count di-nucleotides in file
#' count_kmers(tf, kmer = 2)
#'
#' # Count mono-nucleotide in file
#' count_kmers(tf, kmer = 1)
#' unlink(tf)
count_kmers <- function(filename, kmer) {
  result <- .Call("count_kmers_R", path.expand(as.character(filename)), as.integer(kmer))
  if(!is.null(result)) {
    return(result)
  }
}


#' Calculate enrichments
#'
#' @param testfile Test sequences
#' @param ctrlfile Control sequences
#' @param kmer     Length of kmer
#' @param normalize Get log2 of enrichments
#'
#' @return dataframe containing enrichments
#' @useDynLib rkatss, .registration = TRUE
#' @export
#'
#' @examples
#' print("Hello, world!")
enrichments <- function(testfile, ctrlfile, kmer, normalize = FALSE) {
  result <- .Call("enrichments_R", path.expand(as.character(testfile)), path.expand(as.character(ctrlfile)), as.integer(kmer), normalize)
  if(!is.null(result)) {
    return(result)
  }
}


#' Title Iterative K-mer Knockout Enrichments
#'
#' @param testfile Test sequences file. Can be in FASTQ, FASTA, or Sequences format
#' @param ctrlfile Control sequences file
#' @param kmer Length of k-mer
#' @param iterations Number of iterations to perform
#' @param normalize  Normalize enrichments to log2
#' @param threads    Number of threads to use
#'
#' @return data.frame containing the enrichments
#' @useDynLib rkatss, .registration = TRUE
#' @export
#'
#' @examples
ikke <- function(testfile, ctrlfile, kmer, iterations = 10, normalize = FALSE,
                 threads = 1) {
  testfile <- path.expand(as.character(testfile))
  ctrlfile <- path.expand(as.character(ctrlfile))
  kmer <- as.integer(kmer)
  iterations <- as.numeric(iterations)
  threads <- as.integer(threads)

  result <- .Call("ikke_R", testfile, ctrlfile, kmer, iterations, normalize,
                  threads)
  if(!is.null(result)) {
    return(result)
  }
}


#' Search for a sequence within a sequence
#'
#' This function returns the index in which a sub-sequence is found within a
#' larger sequence, or 0 if not found. Searching is case insensitive and 'T' and
#' 'U' are equivalent. The maximum length sequence to search for is currently
#' capped at 255.
#'
#' @param sequence Big sequence
#' @param search   Small sequence to search for
#'
#' @return Index in which search was found, 0 if not found
#' @useDynLib rkatss, .registration = TRUE
#' @export
#'
#' @examples
#' ## Find "AAA" in a random sequence
#' seq <-paste(sample(c("A","C","G","T"), 1000, TRUE), collapse="")
#' seqseq(seq, "AAA")
#'
#' ## Searching is case-insensitive
#' seq <- "accgtaagggtgccttac"
#' find2 <- "GGGT"
#' seqseq(seq, "GGGT")
#' seqseq(seq, "gCCtT")
#'
#' ## T and U are both interchangeable
#' seqseq(seq, "GCCUU")
#'
#' ## 0 is returned when search is not in sequence
#' seqseq(seq, "AAA")
seqseq <- function(sequence, search, all.matches = FALSE) {
  return(.Call("seqseq_R", as.character(sequence), as.character(search), all.matches))
}
