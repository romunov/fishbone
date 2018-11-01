#' Automatically call or flag alleles from NGS data
#'
#' Call or flag candidate alleles to construct a consensus genotype.
#'
#' Function works on an individual run (one PCR reaction). Apply this to \code{sample * locus * run} combination.
#' The data should come from an NGS run as processed by de Barba et al. (2016).
#'
#' For algorithm used to find stutters, see \code{\link{findStutter}}.
#'
#' De Barba, M., Miquel, C., Lobréaux, S., Quenette, P. Y., Swenson, J. E., & Taberlet, P. (2016).
#' High-throughput microsatellite genotyping in ecology: improved accuracy, efficiency, standardization
#' and success with low-quantity and degraded DNA. Molecular Ecology Resources, 1-16.
#' https://doi.org/10.1111/1755-0998.12594
#'
#' Key for abbreviations used in (pseudo)code:
#' \itemize{
#'   \item A = allele
#'   \item S = stutter
#'   \item R = relative low threshold
#'   \item L = low count threshold
#'   \item D = disbalance
#' }
#' The result is appended three columns; one for called A, one for flagged alleles and if read is a stutter.
#' Possible flags are:
#' \itemize{
#'   \item L = low amplification threshold (if for some reason, number of total reads is very low, alleles get a flag)
#'   \item N = no stutter (if there was enough reads but no stutter was present)
#'   \item D = disbalance - alleles not in balance (expecting 1:1 for heterozygotes, those out of balance flagged)
#'   \item M = multiple alleles (self explanatory)
#' }
#'
#' Algorithm is as follows:
#'
# for each A:
# compare everything according to the highest read count, calculate relative size to \code{maxA -> rs}
#'
#'\itemize{
#'   \item 0. find max allele height
#'   \item 1. if allele has number of reads < L, flag it as "L"
#'   \item 2. see if allele has stutter
#'   \item 2a. if yes, mark as called
#'   \item 2aa. if A in disbalance (A < D), flag as "D"
#'   \item 2ab. mark stutter as such $stutter = TRUE
#'   \item 2b. if no, check AlleleWithNoStutterHeight
#'   \item 2ba. if x > AlleleWithNoStutterHeight, add flag "N"
#'   \item 2bb. if x < AlleleWithNoStutterHeight, ignore allele
#'   \item 3. if number of unflagged alleles is more than 2 (those marked with D are not counted), add flag "M" to all
#' }
#'
#' @param fb A \code{fishbone} object.
#' @param tbase A data.frame with thresholds which are locus specific. Thresholds are:
#' \itemize{
#'   \item stutter (if lower than this, allele is ignored as stutter)
#'   \item disbalance (if heterozygous alleles are not in 1:1 ratio)
#'   \item low count (anything below this threshold gets flagged as light on read count)
#'   \item allele with no stutter height (if no stutter is found, how many reads do we allow for alleles to be called)
#' }
#' @param clean Logical. If \code{TRUE} (default), it will return only called alleles and their stutters.
#' @param verbose Logical. If \code{TRUE}, it will print which sample is being processed.
#'
#' @return Output should be all alleles and their stutters.
#' @export
#' @importFrom data.table ":="
#' @import data.table

callAllele <- function(fb, tbase = NULL, clean = TRUE, verbose = FALSE) {
  # In case the object is not a data.table, make it one.
  if (all(class(fb) != "data.table")) {
    fb <- as.data.table(fb)
  }

  if (is.null(tbase)) stop("Please provide `tbase` object.")
  # This is the function which implements core of the algorithm explained in the help file.
  # This function runs on sample * locus * run combination, which means only one marker per plate.
  locus <- unique(fb$Marker)
  sn <- unique(fb$Sample_Name)
  plate <- unique(fb$Plate)
  rn <- unique(fb$Run_Name)
  ps <- unique(fb$Position)

  proc.info <- sprintf("sample: %s; locus: %s; plate: %s; library: %s", sn, locus, plate, rn)
  if (verbose) {
    message(proc.info)
  }

  stopifnot(length(locus) == 1)
  stopifnot(length(plate) == 1)

  motif <- tbase[tbase$Marker %in% locus, "Repeat"]

  # prepare columns to be used for calling/flagging of allele(s)
  # define because R CMD check produces a NOTE if these variables are not defined somewhere
  called <- NULL
  flag <- NULL
  stutter <- NULL

  fb[, called := FALSE]
  fb[, flag := ""]
  fb[, stutter := FALSE]

  # Prepare an empty object in case there are no reads, one way or the other.
  out.blank <- fb[0, ]

  out.ord <- c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name", "length",
               "Position", "called", "flag", "stutter", "Sequence", "TagCombo")

  out.blank <- out.blank[, out.ord, with = FALSE]

  if (nrow(fb) == 0) {
    return(out.blank)
  }

  # We need this because data.table doesn't support row names, see https://stackoverflow.com/a/24246819/322912
  fb$fbid <- as.character(1:nrow(fb))

  # prepare statistics
  S <- fetchTH(tbase, stat = "Stutter", locus = locus)
  D <- fetchTH(tbase, stat = "Disbalance", locus = locus)
  L <- fetchTH(tbase, stat = "LowCount", locus = locus)
  N <- fetchTH(tbase, stat = "AlleleWithNoStutterHeight", locus = locus)

  # 0. find max allele height
  maxA <- max(fb$Read_Count)

  # 1. if allele has number of reads < L, flag it as "L"
  # (performed towards the end)

  for (i in 1:nrow(fb)) {
    rh <- (fb$Read_Count[i]/maxA)

    # 2. see if allele has stutter
    find.stt <- findStutter(fb[i, ], fb = fb, motif = motif)

    # 2a. if yes, mark as called
    if (is.character(find.stt) && length(find.stt) == 1 && rh > S) {
      fb[i, "called"] <- TRUE

      # 2aa. if A in disbalance (A < D), flag as "D"
      if (rh < D) {
        fb[i, "flag"] <- paste(fb[i, "flag"], "D", sep = "")
      }

      # 2ab. mark stutter as such
      fb[fb$fbid %in% find.stt, "stutter"] <- TRUE
    } else {
      # 2b. if no, check AlleleWithNoStutterHeight
      # 2ba. if x > AlleleWithNoStutterHeight, add flag "N"
      if (fb[i, "Read_Count"] >= N) {
        fb[i, "flag"] <- paste(fb[i, "flag"], "N", sep = "")
        # 2bb. if x < AlleleWithNoStutterHeight, ignore allele
      }
    }
  }

  # 3. if number of unflagged A > 2, add flag "M" to all alleles
  #    Do not include alleles with D as true called alleles, this
  #    should minimize number of called multiple alleles where stutters
  #    are actually higher than expected as specified in `pars.csv` table.
  calledA <- fb$called

  # if alleles have flag D, do not consider them as true alleles because
  # they may be just "too high" stutters
  calledA[grepl("D", fb$flag)] <- FALSE

  if (sum(calledA, na.rm = TRUE) > 2) {
    fb[i = calledA, flag := paste(fb[calledA, j = flag], "M", sep = "")]
  }

  fb <- fb[order(fb$length, decreasing = TRUE), ]

  # Sort columns in a more readable fashion (by e.g. keeping Sequence last).
  fb <- fb[, out.ord, with = FALSE]

  # 1. if allele has number of reads < L, flag it as "L"
  ss <- fb$Read_Count < as.numeric(L) & fb$called == TRUE
  fb[ss, flag := paste(fb$flag[ss], "L", sep = "")]

  # If clean == TRUE, return only sequences which were tagged as allele or stutter
  if (clean) {
    out <- fb[fb$called | fb$stutter, ]
  } else {
    out <- fb
  }

  if (nrow(out) == 0) {
    return(out.blank)
  }
  out
}
