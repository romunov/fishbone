#' Automatically call or flag alleles from NGS data
#'
#' Call or flag candidate alleles to construct a consensus genotype.
#'
#' Function works on an individual run (one PCR reaction). Use R functionality to
#' apply this to sample * locus * run combination. The data should come from an NGS run as
#' processed by de Barba et al. (2016).
#'
#' A note o finding stutters: For each candidate allele the sequence is shortened for the first
#' motif found given the locus.
#'
#' De Barba, M., Miquel, C., Lobréaux, S., Quenette, P. Y., Swenson, J. E., & Taberlet, P. (2016).
#' High-throughput microsatellite genotyping in ecology: improved accuracy, efficiency, standardization
#' and success with low-quantity and degraded DNA. Molecular Ecology Resources, 1–16.
#' https://doi.org/10.1111/1755-0998.12594
#'
#' Key for abbreviations:
#' A = allele
#' S = stutter
#' R = relative low threshold
#' L = low count threshold
#' D = disbalance
#'
#' The result is appended two columns, one for called A and another for flagged.
#' Possible flags are:
#' * L = low amplification threshold
#' * N = no stutter
#' * D = alleles not in balance
#' * M = multiple alleles
#'
#' Algorithm is as follows:
#' - descendingly sort A according to number of reads
#' - find number of reads for highest A -> maxA
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
#' column names.
#' @param tbase A data.frame with thresholds which are locus specific. Thresholds are balance of alleles (B),
#' ignore threshold (IT), low run threshold (LRT) and relative stutter height (S).
#' @param motif A data.frame of loci motifs. Expected two columns with loci names (`locus`) and actual motifs (`motif`).

callGenotype <- function(fb, tbase = NULL, motif = NULL) {
  if (is.null(tbase)) stop("Please provide `tbase` object.")
  if (is.null(motif)) stop("Please provide `motif` object.")
  # This is the function which implements core of the algorithm explained in the help file.

  # This function runs on sample * locus * run combination, which means only one marker per plate.
  locus <- unique(fb$Marker)
  motif <- motif[motif$locus == locus, "motif"]

  stopifnot(length(locus) == 1)
  stopifnot(length(unique(fb$Plate)) == 1)

  # prepare columns to be used for calling/flagging of allele(s)
  fb$call <- ""
  fb$flag <- ""
  fb$stutter <- ""

  # prepare statistics
  S <- fetchTB(tbase, stat = "Stutter", locus = locus)
  D <- fetchTH(tbase, stat = "Disbalance", locus = locus)
  L <- fetchTH(tbase, stat = "LowCount", locus = locus)
  N <- fetchTB(tbase, stat = "AlleleWithNoStutterHeight", locus = locus)

  # 0. find max allele height
  maxA <- max(fb$Read_Count)

  # 1. check that at least max A is above L threshold
  if (maxA <= L) {
    fb$flag <- paste(fb$flag, "L", sep = "")
  }

  for (i in 1:nrow(fb)) {
    # 2. see if allele has stutter
    find.stt <- findStutter(fb[i, ], fb = fb, motif = motif)

    # 2a. if yes, mark as called
    if (is.character(find.stt) && length(find.stt) == 1) {
      fb[i, "call"] <- TRUE

      # 2aa. if A in disbalance (A < D), flag as "D"
      if ((fb[i, "Read_Count"]/maxA) < D) {
        fb[i, "flag"] <- paste(fb[i, "flag"], "D", sep = "")
      }

      # 2ab. mark stutter as such $stutter = TRUE
      fb[rownames(fb) %in % find.stt, "stutter"] <- TRUE
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
  if (sum(fb[, "called"]) > 2) {
    calledA <- fb[, "called"]
    fb[calledA, "flag"] <- paste(fb[calledA, "flag"], "M", sep = "")
  }
  fb
}

#' Find locus specific stutter.
#'
#' From a `fishbone` object, try to find a stutter.
#' @param x A `fishbone` object.
#' @param fb A `fishbone` object with cancidate stutter alleles.
#' @return A character string of the row name of the stutter(s). If stutter is not found, it
#' returns NULL.
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)

findStutter <- function(x, fb, motif) {
  # 1. shorten allele for a given motif
  new.stt <- sub(motif, replacement = "", x = x$Sequence)

  # return FALSE if allele was not shortened into a stutter - indicative that something,
  # somewhere, somehow went wrong
  if (nchar(new.stt) >= x$Sequence) {
    warning(sprintf("findStutter: unable to shorten allele into a stutter (sequence %s of motif %s)",
                    unique(x$Sequence), motif)
    )
    return(FALSE)
  }

  # 2. Compare to candidate stutters and see if it matches structurally
  # candidate stutters of motif length shorter. Since we are using tetra- and penta-
  # repeats, stutters of 2*motif_lengths are not very common.
  motif.length <- nchar(motif)
  cand.stt <- fb[fb$lengths == (x$lengths - motif.length), ]

  if (nrow(cand.stt) == 0) {
    return(NULL)
  }

  found.stt <- cand.stt[new.stt == cand.stt$Sequence, ]

  # 3. return candidate stutter allele rowname
  if (nrow(found.stt) == 0) {
    message(sprintf("No stutter from allele %s", x$Allele))
    return(NULL)
  }
  if (nrow(found.stt) == 1) {
    message(sprintf("From allele %s, stutter %s", x$Allele, rownames(found.stt)))
    return(rownames(found.stt))
  }
  if (nrow(found.stt) > 1) {
    warning(sprintf("findStutter: Found multiple stutters for sample %s", x$Sample_Name))
    return(rownames(found.stt))
  }
}

#' Fetch threshold values from some database.
#' @param x A slice of `fishbone` object.
#' @param stat Which statistic to fetch, e.g. "RelativeLow", "Disbalance", "Stutter", "LowCount".
#' @param locus Character. For which locus to fetch statistics?
#' @author Roman Lustrik (roman.lustri@@biolitika.si)

fetchTH <- function(x, stat, locus) {
  marker <- gsub("^.*(\\d+)$")
  out <- x[x$Marker %in% locus, ..stat]
  if (length(out) != 1) stop(sprintf("Unable to fetch %s for locus %s.", stat, locus))
  out
}
