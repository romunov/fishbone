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
#' LRT = low run threshold
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
#' compare everything according to the highest read count, calculate relative size to maxA -> rs
#' 1. if signal is below ignore threshold (IT), remove this allele
#' 2. if below low run threshold (LRT), add flag "L"
#' 3. if stutter for A found, call, otherwise flag as "N"
#' 4. if rs <= B, add flag "B"
#' 5. if number of unflagged A > 2, add flag "M" to all alleles
#'
#' @param fb A `fishbone` object.
#' column names.
#' @param tbase A data.frame with thresholds which are locus specific. Thresholds are balance of alleles (B),
#' ignore threshold (IT), low run threshold (LRT) and relative stutter height (S).
#' @param motif A data.frame of loci motifs. Expected two columns with loci names (`locus`) and actual motifs (`motif`).

callGenotype <- function(fb, tbase, motif) {
  # This is the function which implements core of the algorithm explained in the help file (or see
  # roxygen2 comments above).

  # To process data row-wise, we can spit it first for convenience (because apply(x, MARGIN = 1, ...)
  # coerces input to vector).
  x.split <- split(fb, f = 1:nrow(fb), drop = TRUE)

  cg <- function(x, maxA, fb, tbase) {
    id <- rownames(x)
    fb <- fb[!(rownames(fb) %in% id), ]

    # prepare thresholds
    locus <- unique(fb$Marker)
    stopifnot(length(locus) == 1)
    IT <- fetchTH(tbase, stat = "IT", locus = locus)
    B <- fetchTH(tbase, stat = "B", locus = locus)
    S <- fetchTH(tbase, stat = "S", locus = locus)
    LRT <- fetchTH(tbase, stat = "LRT", locus = locus)

    # prepare columns to be used for calling/flagging of allele(s)
    x$call <- NA
    x$flag <- ""

    # compare everything according to the highest read count, calculate relative size to maxA -> rs
    rs <- x$lengths/maxA

    # 1. if signal is below ignore threshold (IT), remove this allele
    if (rs <= IT) {
      return(NA)
    }

    # 2. if below low run threshold (LRT), add flag "L"
    if (rl <= LRT) {
      x$flag <- paste(x$flag, "L", sep = "")
    }

    # 3. if stutter for A found, call, otherwise flag as "N"
    find.stt <- findStutter(x, locus = locus, fb = fb, motif = motif[motif$locus == locus])

    if (find.stt == TRUE) {
      x$call <- "called"
    } else {
      x$flag <- paste(x$flag, "N", sep = "")
    }

    # 4. if rs <= B, add flag "B"
    # 5. if number of unflagged A > 2, add flag "M" to all alleles
  }
  sapply(x.split, FUN = cg, maxA = max(x$Read_Count), fb = fb, tbase = tbase)
}

#' Find locus specific stutter.
#'
#' From a `fishbone` object, try to find a stutter.
#' @param x A `fishbone` object.
#' @param locus Character. Since threshold values can be locus specific, which locus?
#' @param fb A `fishbone` object with cancidate stutter alleles.
findStutter <- function(x, locus, fb, motif) {
  # motif = "attt"
  # fb = data.frame(Sample_Name, Allele, lengths...)
  # locus = "03"

  # 1. fetch motif for given locus
  mt <- motif[motif$locus = locus, "motif"]
  # 2. shorten allele for a given motif
  new.stt <- sub(mt, replacement = "", x = x$Sequence)

  # issue a warning if allele was not shortened into a stutter - indicative that something,
  # somewhere, somehow went wrong
  if (nchar(new.stt) >= x$Sequence) {
    warning(sprintf("findStutter: unable to shorten allele into a stutter (sequence %s of locus %s and motif %s)",
                    x$Sequence, locus, motif))
  }

  # 3. compare to candidate stutters and see if it matches structurally
  # candidate stutters can be n*length of motif shorter
  motif.length <- nchar(motif)
  cand.stt <- x[x$lengths == (x$lengths - motif.length), ]

  if (nrow(cand.stt) < 1) {
    return("No")
  }

  found.stt <- cand.stt[new.stt == cand.stt$Sequence, ]

  # 4. return candidate stutter allele
  if (nrow(found.stt) == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}
#' Fetch threshold values from some database.
#' @param x A slice of `fishbone` object.
#' @param stat Which statistic to fetch, e.g. "IT", "B", "S", "LRT".
#' @param locus Character. For which locus to fetch statistics?
#' @author Roman Lustrik (roman.lustri@@biolitika.si)
fetchTH <- function(x, stat, locus) {
  out <- x[x$stat == stat & x$L == locus, "value"]
  stopifnot(nrow(out) == 1)
  out
}
