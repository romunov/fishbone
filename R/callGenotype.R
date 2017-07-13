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
#' 3. if stutter for A found (relative height correct), call, otherwise flag as "N"
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

  # prepare columns to be used for calling/flagging of allele(s)
  fb$call <- NA
  fb$flag <- ""

  # Check for LRT from the get go. A calling proceeds like for any other run.
  # 2. if below low run threshold (LRT), add flag "L"
  locus <- unique(fb$Marker)
  stopifnot(length(locus) == 1)
  LRT <- fetchTH(tbase, stat = "LRT", locus = locus)

  if (max(fb$Read_count) <= LRT) {
    fb$flag <- paste(fb$flag, "L", sep = "")
  }

  # To process data row-wise, we can spit it first for convenience (because apply(x, MARGIN = 1, ...)
  # coerces input to vector which can hold only one type).
  x.split <- split(fb, f = 1:nrow(fb), drop = TRUE)

  cg <- function(x, maxA, fb, tbase) {
    id <- rownames(x)
    fb <- fb[!(rownames(fb) %in% id), ]

    # prepare thresholds
    locus <- unique(fb$Marker)
    stopifnot(length(locus) == 1)
    IT <- fetchTH(tbase, stat = "IT", locus = locus)
    B <- fetchTH(tbase, stat = "B", locus = locus)
    LRT <- fetchTH(tbase, stat = "LRT", locus = locus)

    # calculate relative size to maxA -> rs
    rs <- x$lengths/maxA

    # 1. If signal is below the ignore threshold (IT), remove/ignore A.
    if (rs <= IT) {
      return(NA)
    }

    # 2. It makes sense to do this step before iterating over all alleles. See code in the
    # immediate definition of `callGenotype`.

    # 3. Find structural stutter for A. If found, call it, otherwise flag as "N".
    find.stt <- findStutter(x, locus = locus, fb = fb, motif = motif[motif$locus == locus])

    if (is.character(find.stt)) {
      x$call <- "called"
    } else {
      x$flag <- paste(x$flag, "N", sep = "")
    }

    # 4. Check if alleles are in high disequilibrium. We would expect for heterozygots
    # to have equal number of reads per allele, but do not due to whatever reasons.
    if (rs <= B) {
      x$flag <- paste(x$flag, "B", sep = "")
    }
    x
  }

  run.A <- sapply(x.split, FUN = cg, maxA = max(x$Read_Count), fb = fb, tbase = tbase, simplify = FALSE)

  # Remove all A which are flagged as non A
  run.A <- run.A[!sapply(run.A, FUN = is.na)]

  run.A <- do.call(rbind, run.A)

  # 5. if number of unflagged A > 2, add flag "M" to all alleles
  if (length(run.A) > 2) {
    run.A$flag <- paste(run.A$flag, "M", sep = "")
  }
  run.A
}

#' Find locus specific stutter.
#'
#' From a `fishbone` object, try to find a stutter.
#' @param x A `fishbone` object.
#' @param locus Character. Since threshold values can be locus specific, which locus?
#' @param fb A `fishbone` object with cancidate stutter alleles.
#' @return A character string of the row name of the stutter(s).
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)

findStutter <- function(x, locus, fb, motif) {
  # motif = "attt"
  # fb = data.frame(Sample_Name, Allele, lengths...)
  # locus = "03"

  # 1. fetch motif for given locus
  mt <- motif[motif$locus == locus, "motif"]
  # 2. shorten allele for a given motif
  new.stt <- sub(mt, replacement = "", x = x$Sequence)

  # return FALSE if allele was not shortened into a stutter - indicative that something,
  # somewhere, somehow went wrong
  if (nchar(new.stt) >= x$Sequence) {
    warning(sprintf("findStutter: unable to shorten allele into a stutter (sequence %s of locus %s and motif %s)",
                    unique(x$Sequence), locus, motif)
    )
    return(FALSE)
  }

  # 3. Compare to candidate stutters and see if it matches structurally
  # candidate stutters of motif length shorter. Since we are using tetra- and penta-
  # repeats, stutters of 2*motif_lengths are not very common.
  motif.length <- nchar(motif)
  cand.stt <- x[x$lengths == (x$lengths - motif.length), ]

  if (nrow(cand.stt) < 1) {
    return(FALSE)
  }

  found.stt <- cand.stt[new.stt == cand.stt$Sequence, ]

  # 4. return candidate stutter allele rowname
  if (nrow(found.stt) == 1) {
    return(rownames(found.stt))
  } else {
    warning(sprintf("findStutter: Found multiple stutters for sample %s and locus %s",
                    x$Sample_Name, locus))
    return(rownames(found.stt))
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
