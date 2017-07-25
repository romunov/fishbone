#' Automatically call or flag alleles from NGS data
#'
#' Call or flag candidate alleles to construct a consensus genotype.
#'
#' Function works on an individual run (one PCR reaction). Apply this to sample * locus * run combination.
#' The data should come from an NGS run as processed by de Barba et al. (2016).
#'
#' For algorithm used to find stutters, see `?findStutter`.
#'
#' De Barba, M., Miquel, C., Lobréaux, S., Quenette, P. Y., Swenson, J. E., & Taberlet, P. (2016).
#' High-throughput microsatellite genotyping in ecology: improved accuracy, efficiency, standardization
#' and success with low-quantity and degraded DNA. Molecular Ecology Resources, 1–16.
#' https://doi.org/10.1111/1755-0998.12594
#'
#' Key for abbreviations used in (pseudo)code:
#' A = allele
#' S = stutter
#' R = relative low threshold
#' L = low count threshold
#' D = disbalance
#'
#' The result is appended three columns; one for called A, one for flagged alleles and if read is a stutter.
#' Possible flags are:
#' * L = low amplification threshold (if for some reason, number of total reads is very low, alleles get a flag)
#' * N = no stutter (if there was enough reads but no stutter was present)
#' * D = disbalance - alleles not in balance (expecting 1:1 for heterozygotes, those out of balance flagged)
#' * M = multiple alleles (self explanatory)
#'
#' Algorithm is as follows:
#'
# for each A:
# compare everything according to the highest read count, calculate relative size to maxA -> rs
#'
#' 0. find max allele height
#' 1. check that at least max A is above L threshold
#' 1a. if not, flag all alleles with "L"
#' 2. see if allele has stutter
#' 2a. if yes, mark as called
#' 2aa. if A in disbalance (A < D), flag as "D"
#' 2ab. mark stutter as such $stutter = TRUE
#' 2b. if no, check AlleleWithNoStutterHeight
#' 2ba. if x > AlleleWithNoStutterHeight, add flag "N"
#' 2bb. if x < AlleleWithNoStutterHeight, ignore allele
#' 3. if number of unflagged alleles is more than 2, add flag "M" to all
#'
#' Output should be all alleles and their stutters.
#'
#' @param fb A `fishbone` object.
#' @param tbase A data.frame with thresholds which are locus specific. Thresholds are:
#' - stutter (if lower than this, allele is ignored as stutter)
#' - disbalance (if heterozygous alleles are not in 1:1 ratio)
#' - low count (anything below this threshold gets flagged as light on read count)
#' - allele with no stutter height (if no stutter is found, how many reads do we allow for alleles
#' to be called)
#' @param clean Logical. If TRUE (default), it will return only called alleles and their stutters.
#' @export
#' @importFrom data.table ":="
#' @import data.table

callAllele <- function(fb, tbase = NULL, clean = TRUE) {
  if (is.null(tbase)) stop("Please provide `tbase` object.")
  # This is the function which implements core of the algorithm explained in the help file.
  # This function runs on sample * locus * run combination, which means only one marker per plate.
  locus <- unique(fb$Marker)
  stopifnot(length(locus) == 1)
  stopifnot(length(unique(fb$Plate)) == 1)

  motif <- tbase[tbase$Marker %in% locus, "Repeat"]

  # prepare columns to be used for calling/flagging of allele(s)
  fb[, called := FALSE]
  fb[, flag := ""]
  fb[, stutter := FALSE]

  # browser()
  # We need this because data.table doesn't support row names, see https://stackoverflow.com/a/24246819/322912
  fb$fbid <- as.character(1:nrow(fb))

  # prepare statistics
  S <- fetchTH(tbase, stat = "Stutter", locus = locus)
  D <- fetchTH(tbase, stat = "Disbalance", locus = locus)
  L <- fetchTH(tbase, stat = "LowCount", locus = locus)
  N <- fetchTH(tbase, stat = "AlleleWithNoStutterHeight", locus = locus)

  # 0. find max allele height
  maxA <- max(fb$Read_Count)

  # 1. check that at least max A is above L threshold
  if (maxA <= L) {
    fb$flag <- paste(fb$flag, "L", sep = "")
  }

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
  calledA <- fb$called
  if (sum(calledA, na.rm = TRUE) > 2) {
    fb[i = calledA, flag := paste(fb[calledA, j = flag], "M", sep = "")]
  }

  fb <- fb[order(fb$length, decreasing = TRUE), ]
  out.ord <- c("Sample_Name", "Plate", "Read_Count", "Marker", "Run_Name", "length",
               "Position", "called", "flag", "stutter", "Sequence")
  fb <- fb[, ..out.ord]

  if (clean) {
    fb[fb$called | fb$stutter, ]
  } else {
    fb
  }
}
