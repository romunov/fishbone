#' Automatically call or flag alleles
#'
#' Call or flag candidate alleles to construct a consensus genotype.
#'
#' Function works on an individual run (one PCR reaction). Use R functionality to
#' apply this to sample * locus * run combination.
#'
#' Key for abbreviations:
#' A = allele
#' S = stutter
#' IT = ignore threshold
#' LRT = low run
#' B = balance
#'
#' The result is appended two columns, one for called A and another for flagged.
#' Possible flags are:
#' * L = low run threshold
#' * N = no stutter
#' * B = alleles not in balance
#' * M = multiple alleles
#'
#' Algorithm is as follows:
#' - descendingly sort A according to number of reads
#' - find number of reads for highest A -> maxA
#'
#' for each A:
#' 1. compare everything according to the highest read count, calculate relative size to maxA -> rs
#' 2. if signal is below ignore threshold (IT), remove this allele
#' 3. if below low run threshold (LRT), add flag "L"
#' 4. if stutter for A found, call, otherwise flag as "N"
#' 5. if rs <= B, add flag "B"
#' 6. if number of unflagged A > 2, add flag "M" to all alleles
#'
#' @param fb A `fishbone` object.
#' column names.
#' @param tbase A data.frame with thresholds which are locus specific. Thresholds are balance of alleles (B),
#' ignore threshold (IT), low run threshold (LRT) and relative stutter height (S)

callGenotype <- function(fb, tbase) {
  # This is the function which implements core of the algorithm explained in the help file (or see
  # roxygen2 comments above).

  # To process data row-wise, we can spit it first for convenience (because apply(x, MARGIN = 1, ...)
  # coerces input to vector).
  x.split <- split(fb, f = 1:nrow(fb), drop = TRUE)

  cg <- function(x, maxA, fb) {
    # prepare thresholds
    locus <- unique(fb$Marker)
    stopifnot(length(locus) == 1)
    IT <- fetchTH(tbase, stat = "IT", locus = locus)
    B <- fetchTH(tbase, stat = "B", locus = locus)
    S <- fetchTH(tbase, stat = "S", locus = locus)
    LRT <- fetchTH(tbase, stat = "LRT", locus = locus)

    # 1. compare everything according to the highest read count, calculate relative size to maxA -> rs
    rl <- x/maxA

    if (rl <= IT) {
      x
    }
  }
  sapply(x.split, FUN = cg, maxA = max(x$Read_Count))
}

fetchTH <- function(x, stat, locus) {
  out <- x[x$stat == stat & x$L == locus, "value"]
  stopifnot(nrow(out) == 1)
  out
}