#' K-mer counting
#'
#' Count the k-mers in a fastq, fasta, or raw sequences file.
#'
#' @param file Name of the file which you want to count k-mers from
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
count_kmers <- function(file, kmer = 3) {
  if(!is.character(file))
    stop("file must be a character string")
  if(!is.numeric(kmer) && kmer %% 1 != 0) {
    stop("kmer must be an integer")
  }

  file <- path.expand(as.character(file))
  return(.Call("count_kmers_R", file, as.integer(kmer)))
}


#' Calculate k-mer enrichments
#'
#' @param testfile Test sequences
#' @param ctrlfile Control sequences
#' @param kmer     Length of kmer
#' @param normalize Get log2 of enrichments
#' @param probabilistic Compute probabilistic enrichments
#' @param verbose Verbose output for calculations
#'
#' @return dataframe containing enrichments
#' @useDynLib rkatss, .registration = TRUE
#' @export
#'
#' @examples
#' print("Hello, world!")
enrichments <- function(testfile, ctrlfile = NULL, kmer = 3, probabilistic = FALSE,
                        normalize = FALSE, verbose = FALSE) {
  if(!is.character(testfile))
    stop("testfile must be a character string")
  if(!is.character(ctrlfile) && !is.null(ctrlfile))
    stop("ctrlfile must be a character string or NULL")
  if(!is.numeric(kmer) || kmer %% 1 != 0)
    stop("kmer must be an integer")
  if(!is.logical(probabilistic))
    stop("probabilistic must be either TRUE or FALSE")
  if(!is.logical(normalize))
    stop("normalize must be either TRUE or FALSE ")
  if(!is.logical(verbose))
    stop("verbose must be either TRUE or FALSE")

  if(probabilistic && is.character(ctrlfile))
    warning("Ignoring ctrlfile argument")
  if(!probabilistic && is.null(ctrlfile))
    stop("ctrlfile is required when computing non-probabilistic enrichments")
  if(16 >= kmer && kmer>12) {
    menu_title = paste(convert_bytes(4^kmer * 176), "Are you sure you want to proceed?")
    if(utils::menu(c("Yes", "No! Fix your program!"), title = menu_title) == 2)
      return(NULL)
  }

  # Done with argument checks, expand filepaths if necessary
  testfile <- path.expand(as.character(testfile))
  ctrlfile <- path.expand(as.character(ctrlfile))

  # Do calculations in C and return
  return(.Call("enrichments_R", testfile, ctrlfile, as.integer(kmer),
               probabilistic, normalize, verbose))
}


#' Title Iterative K-mer Knockout Enrichments
#'
#' @param testfile Test sequences file. Can be in FASTQ, FASTA, or raw sequences
#' format. Raw sequences format is a file containing only "A", "C", "G", and "T"
#' /"U" characters, in every sequence separated by newline.
#' @param ctrlfile Control sequences file. Can be in FASTQ, FASTA, or raw sequences
#' format. Raw sequences format is a file containing only "A", "C", "G", and "T"
#' /"U" characters, in every sequence separated by newline.
#' @param kmer Length of k-mer.
#' @param iterations Number of iterations to perform
#' @param normalize  Normalize enrichments to log2
#' @param threads    Number of threads to use. Specifying less than 1 thread
#' sets the number of threads as 1.
#' @param probabilistic Calculate probabilistic enrichments.
#'
#' @return data.frame containing the enrichments
#' @useDynLib rkatss, .registration = TRUE
#' @export
#'
#' @examples
#' print("ikke example")
ikke <- function(testfile, ctrlfile = NULL, kmer = 3, iterations = 10,
                 probabilistic = FALSE, normalize = FALSE, threads = 1) {
  if(!is.character(testfile))
    stop("testfile must be a character string")
  testfile <- path.expand(testfile)

  if(!is.character(ctrlfile) && !is.null(ctrlfile))
    stop("ctrlfile must be a character string")
  ctrlfile <- path.expand(ctrlfile)

  if(!is.numeric(kmer) && kmer != as.integer(kmer))
    stop("kmer must be an integer")
  kmer <- as.integer(kmer)

  if(!is.numeric(iterations) && iterations %% 1 != 0)
    stop("iterations must be an integer")
  iterations <- as.numeric(iterations)

  if(!is.logical(probabilistic))
    stop("probabilistic must be either TRUE or FALSE")

  if(!is.numeric(threads) && kmer != as.integer(threads))
    stop("threads must be an integer")
  threads <- as.integer(threads)

  if(probabilistic && !is.null(ctrlfile))
    warning("Ignoring ctrlfile argument")
  if(!probabilistic && is.null(ctrlfile))
    stop("ctrlfile is required when using non-probabilistic ikke")

  # Arguments seem correct, begin function call
  result <- .Call("ikke_R", testfile, ctrlfile, kmer, iterations, probabilistic,
                  normalize, threads)
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
#' @param sequence    Big sequence
#' @param search      Small sequence to search for
#' @param all.matches Find all occurrences of search in sequence
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
  if(!is.character(sequence))
    stop("sequence must be a character string")
  if(!is.character(search))
    stop("search must be a character string")
  if(!is.logical(all.matches))
    stop("all.matches must be either TRUE or FALSE")

  return(.Call("seqseq_R", as.character(sequence), as.character(search), all.matches))
}


convert_bytes <- function(total_bytes) {
  # Conversion factors
  bytes_in_gb <- 1024^3
  bytes_in_mb <- 1024^2
  bytes_in_kb <- 1024

  # Calculating GB, MB, and KB
  gb <- floor(total_bytes / bytes_in_gb)
  remainder_after_gb <- total_bytes %% bytes_in_gb

  mb <- floor(remainder_after_gb / bytes_in_mb)
  remainder_after_mb <- remainder_after_gb %% bytes_in_mb

  kb <- ceiling(remainder_after_mb / bytes_in_kb)

  # Prepare output data frame
  result <- data.frame(gb = NA, mb = NA, kb = kb)
  if (gb > 0) {
    result$gb <- gb
  }
  if (mb > 0) {
    result$mb <- mb
  }

  # Creating the warning message with the sizes
  warning_msg <- paste("Dataframe being output will be ",
                       if (!is.na(result$gb)) paste(result$gb, "GB") else NULL,
                       if (!is.na(result$gb) && !is.na(result$mb)) " " else NULL,
                       if (!is.na(result$mb)) paste(result$mb, "MB") else NULL,
                       if ((!is.na(result$gb) || !is.na(result$mb)) && result$kb > 0) " " else NULL,
                       if (result$kb > 0) paste("",result$kb, "KB") else NULL,
                       ".", sep = "")

  return(warning_msg)
}
