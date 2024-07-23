################################################################################
#                               PUBLIC FUNCTIONS                               #
################################################################################


#' Get the positional weight matrix (PWM) from enriched k-mers
#'
#' @description
#' get_pwms is a simple method to compute the PWMs from k-mer enrichments. To
#' accomplish this, k-mers are first aligned to the most enriched k-mer. K-mers
#' that do not align to the most enriched k-mer are then used to create another
#' PWM. Weight is determined as the enrichment of the k-mer, non significant
#' positions are trimmed, and the PWM is normalized per position.
#'
#' @param df Dataframe containing k-mers and associated enrichments. Meant to
#' complement enrichments/ikke function, so use their output.
#' @param num_pwms Number of PWMs to generate from `df`
#' @param limit Number of rows to use from `df`
#' @param limit_per_logo Number of rows to use per logo
#' @param type Sequence type. Can be one of either 'DNA' or 'RNA'
#'
#' @return List of PWMs
#' @export
#'
#' @examples
#' # Load sample data
#' data(rbfox2_enrichments)
#'
#' # Get the PWM for the enrichments
#' pwm <- get_pwms(rbfox2_enrichments)
#'
#' # Get the primary and secondary PWMs
#' pwm <- get_pwms(rbfox2_enrichments, num_pwms = 2)
#'
#' # Limit the total number of enrichments to use in PWM
#' pwm <- get_pwms(rbfox2_enrichments, limit = 20)
#'
#' # Use 'U' instead of 'T' for RNA-based PWM
#' pwm <- get_pwms(rbfox2_enrichments, type = "RNA")
get_pwms <- function(df, num_pwms=1, limit=Inf, limit_per_logo=Inf, type=c("DNA", "RNA")) {
  # Ensure the dataframe has the correct structure
  if (!all(c("kmer", "rval") %in% colnames(df))) {
    stop("Dataframe must contain 'kmer' and 'rval' columns")
  }
  if(!typeof(df$kmer) == "character") {
    stop("'kmer' column must be of type \"character\"")
  }
  if(!typeof(df$rval) == "double") {
    stop("'rval' column must be of type \"double\"")
  }

  # Trim the dataframe
  df <- df[1:min(limit,nrow(df)),]

  # Get the primary PWM
  df_aligned <- align_kmers(df, limit = limit_per_logo)
  logo_list <- list(create_sequence_logo(df_aligned, type = type[1]))

  # Get the alternative PWMs
  cur_logo <- 1
  while(cur_logo < num_pwms) {
    df_aligned$kmer <- gsub("-","",df_aligned$kmer)
    df <- subset(df, !df$kmer %in% df_aligned$kmer)
    if(nrow(df) == 0) {
      return(logo_list)
    }
    df_aligned <- align_kmers(df, limit = limit_per_logo)
    logo_list[[length(logo_list)+1]] <- create_sequence_logo(df_aligned, type=type[1])
    cur_logo <- cur_logo + 1
  }

  return(logo_list)
}

#' Plot sequence logo from a PWM
#'
#' @description
#' Plots PWMs generated from \code{\link{get_pwms}}. This is a simple interface
#' over the ggseqlogo package. For more complex plots or greater customization,
#' see \link[ggseqlogo]{ggseqlogo} and \link[ggseqlogo]{geom_logo}
#'
#' @param pwm PWM you want to plot
#' @param ncol Number of columns
#' @param method Use bits of probability method for sequence logo.
#' @param title Title of the plot
#' @param subtitle Title for each PWM. Should be same length as the number of PWMs.
#' @param name Set name of y axis.
#'
#' @return ggplot object of PWM plot
#' @export
#' @import ggplot2
#' @import ggseqlogo
#'
#' @examples
#' # Load sample data
#' data(rbfox2_pwms)
#'
#' # Plot single logo
#' rbfox2_plot <- plot_logo(rbfox2_pwms[[1]])
#' print(rbfox2_plot)
#'
#' # Plot single logo using probability
#' rbfox2_plot <- plot_logo(rbfox2_pwms[[1]], method = "probability")
#' print(rbfox2_plot)
#'
#' # Plot multiple logos
#' rbfox2_plot <- plot_logo(rbfox2_pwms)
#' print(rbfox2_plot)
plot_logo <- function(pwm, ncol = NULL, method = c("bits","probability"),
                      title = "Weighted Sequence Logo", name = NULL, subtitle = NULL) {
  if(is.null(subtitle)) {
    facet.text <- element_blank()
  } else {
    if(!length(pwm)==length(subtitle)) {
      stop("'subtitle' and 'pwm' must be the same length")
    }
    names(pwm) <- subtitle
    facet.text <- element_text(face="bold", angle=0, vjust=0.5)
  }

  if(is.null(title)) {
    title.info <- element_blank()
  } else {
    title.info <- element_text(hjust = 0.5)
  }

  if(!is.null(name)) {
    name.info <- element_text(face="bold", angle=0, vjust = 0.5)
  } else {
    name <- method
    name.info <- element_text()
  }

  sequence_plot <- ggseqlogo(pwm, ncol = ncol, method = method[1]) +
    ggtitle(title) +
    scale_y_continuous(expand = expansion(mult = c(0, .1))) +
    ylab(name) +
    theme(plot.title = title.info,
          axis.title.y = name.info,
          axis.line = element_line(color="black", linewidth=0.5, linetype="solid"),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          strip.background = element_blank(),
          strip.text = facet.text,
          strip.text.y.left = facet.text,
          strip.placement = "outside")

  return(sequence_plot)
}


#' Align k-mers in a dataframe
#'
#' @param df Dataframe containing k-mers
#' @param is_log2 Get 2^enrichments
#' @param limit Limit the number of k-mers to align
#'
#' @return Dataframe with aligned k-mers
#' @export
#' @importFrom stringr str_pad
#' @importFrom stats na.omit
#'
#' @examples
#' # Load sample data
#' data(rbfox2_enrichments)
#'
#' # Align all k-mers in data.frame
#' aligned <- align_kmers(rbfox2_enrichments)
#' head(aligned)
#'
#' # Align only the first five k-mers in data.frame
#' aligned <- align_kmers(rbfox2_enrichments, limit = 5)
#' head(aligned)
align_kmers <- function(df, is_log2 = FALSE, limit = Inf) {
  # Ensure the dataframe has the correct structure
  if (!all(c("kmer", "rval") %in% colnames(df))) {
    stop("Dataframe must contain 'kmer' and 'rval' columns")
  }

  df <- df[1:min(nrow(df),limit),]
  df <- na.omit(df)
  topk <- df[which.max(df[['rval']]),1]
  klen <- nchar(topk)

  align_kmer <- function(current_kmer) {
    if(hamming_distance(topk, current_kmer) <= 1) {
      return(current_kmer)
    }

    # Compare with 1 pad to the left
    curk <- str_pad(current_kmer, width = klen+1, side = "left", pad = "-")
    curk_str <- substr(curk, 2, klen)
    topk_str <- substr(topk, 2, klen)
    if(hamming_distance(topk_str, curk_str) <= 1) {
      return(curk)
    }

    # Compare with 1 pad to the right
    curk <- str_pad(current_kmer, width = klen+1, side = "right", pad = "-")
    curk_str <- substr(curk, 2, klen)
    topk_str <- substr(topk, 1, klen-1)
    if(hamming_distance(topk_str, curk_str) <= 1) {
      return(curk)
    }

    # Compare with 2 pads to the left
    curk <- str_pad(current_kmer, width = klen+2, side = "left", pad = "-")
    curk_str <- substr(curk, 3, klen)
    topk_str <- substr(topk, 3, klen)
    if(hamming_distance(topk_str, curk_str) == 0) {
      return(curk)
    }

    # Compare with 2 pads to the right
    curk <- str_pad(current_kmer, width = klen+2, side = "right", pad = "-")
    curk_str <- substr(curk, 3, klen)
    topk_str <- substr(topk, 1, klen-2)
    if(hamming_distance(topk_str, curk_str) == 0) {
      return(curk)
    }

    return(NA)
  }

  df$kmer <- sapply(df$kmer, align_kmer)
  df <- na.omit(df)

  # Calculate max_lpad and max_rpad
  calculate_padding <- function(kmer) {
    lpad <- 0
    rpad <- 0

    i <- 1
    while(substr(kmer,i,i) == "-") {
      i <- i+1
      lpad <- lpad+1
    }
    i <- nchar(kmer)
    while(substr(kmer,i,i) == "-") {
      i <- i-1
      rpad <- rpad+1
    }

    return(list(lpad = lpad, rpad = rpad))
  }

  paddings <- lapply(df$kmer, calculate_padding)
  max_lpad <- max(sapply(paddings, function(x) x$lpad))
  max_rpad <- max(sapply(paddings, function(x) x$rpad))
  df$kmer <- lapply(df$kmer, pad_kmer, max_lpad=max_lpad, max_rpad=max_rpad)

  df$rval <- as.numeric(df$rval)
  if(is_log2) {
    df["rval"] <- 2^df["rval"]
  }

  return(df)
}

################################################################################
#                               HELPER FUNCTIONS                               #
################################################################################


create_sequence_logo <- function(df, type = c("DNA", "RNA")) {
  # Ensure the dataframe has the correct structure
  if (!all(c("kmer", "rval") %in% colnames(df))) {
    stop("Dataframe must contain 'kmer' and 'rval' columns")
  }

  # Initialize a position weight matrix (PWM)
  kmer_length <- nchar(df$kmer[1])
  pwm <- matrix(0, nrow = 4, ncol = kmer_length)
  if(type[1] == "RNA") {
    rownames(pwm) <- c("A", "C", "G", "U")
  } else {
    rownames(pwm) <- c("A", "C", "G", "T")
  }

  # Substitute all 'T' with 'U'
  if(type[1] == "RNA") {
    df$kmer <- gsub('T', 'U', df$kmer)
  } else {
    df$kmer <- gsub('U', 'T', df$kmer)
  }

  # Fill the weighted frequency matrix
  for (i in 1:nrow(df)) {
    kmer <- unlist(strsplit(as.character(df$kmer[i]), split = ""))
    weight <- df$rval[i]
    for (j in 1:length(kmer)) {
      if(kmer[j] == "-") next
      pwm[kmer[j], j] <- pwm[kmer[j], j] + weight
    }
  }

  # Kill unaligned values
  r_threshold <- max(colSums(pwm))/3
  for(i in kmer_length:1) {
    if(substr(df$kmer[1],i,i) == "-" && sum(pwm[,i]) < r_threshold) {
      pwm <- pwm[,-i]
    }
  }

  # Normalize PWM
  for(i in 1:ncol(pwm)) {
    pwm[,i] <- pwm[,i]/sum(pwm[,i])
  }

  return(pwm)
}


#' Hamming distance of two strings
#'
#' @param str1 First string
#' @param str2 Second string
#'
#' @return Hamming distance between both strings
hamming_distance <- function(str1, str2) {
  return(sum(strsplit(str1,"")[[1]]!=strsplit(str2,"")[[1]]))
}


#' pad a kmer string
#'
#' @param kmer      String to pad
#' @param max_lpad  Amount to pad to the left
#' @param max_rpad  Amount to pad to the right
#'
#' @importFrom stringr str_pad
#'
#' @return Padded string `kmer`
pad_kmer <- function(kmer, max_lpad, max_rpad) {
  lpad <- 0
  rpad <- 0

  i <- 1
  while (substr(kmer,i,i) == "-") {
    i <- i+1
    lpad <- lpad+1
  }
  i <- nchar(kmer)
  while(substr(kmer,i,i) == "-") {
    i <- i-1
    rpad <- rpad+1
  }

  if(max_lpad-lpad)
    kmer <- str_pad(kmer, width = (nchar(kmer)+(max_lpad-lpad)), side = "right", pad = "-")
  kmer <- str_pad(kmer, width = (nchar(kmer)+(max_rpad-rpad)), side = "left", pad = "-")

  return(kmer)
}

################################################################################
#                              KATSS SAMPLE DATA                               #
################################################################################

#' Data frame containing the 5-mer enrichments of RNA-binding protein RBFOX2
#'
#' @name rbfox2_enrichments
#' @docType data
#' @keywords data
#' @usage data(rbfox2_enrichments)
"rbfox2_enrichments"

#' List of matrices containing the 5-mer positional weight matrix of RNA-binding
#' protein RBFOX2
#'
#' @name rbfox2_pwms
#' @docType data
#' @keywords data
#' @usage data(rbfox2_pwms)
"rbfox2_pwms"
