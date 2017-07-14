#' Automatically call or flag alleles from NGS data
#'
#' Call or flag candidate alleles to construct a consensus genotype.
#'
#' Function works on an individual run (one PCR reaction). Use R functionality to
#' apply this to sample * locus * run combination. The data should come from an NGS run as
#' processed by de Barba et al. (2016).
#'
#' De Barba, M., Miquel, C., Lobréaux, S., Quenette, P. Y., Swenson, J. E., & Taberlet, P. (2016).
#' High-throughput microsatellite genotyping in ecology: improved accuracy, efficiency, standardization
#' and success with low-quantity and degraded DNA. Molecular Ecology Resources, 1–16.
#' https://doi.org/10.1111/1755-0998.12594
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

callGenotype <- function(fb, tbase = NULL, motif = NULL) {
  if (is.null(tbase)) stop("Please provide `tbase` object.")
  if (is.null(motif)) stop("Please provide `motif` object.")

  # This is the function which implements core of the algorithm explained in the help file (or see
  # roxygen2 comments above).

  # This function runs on sample * locus * run combination, which means only one marker per plate.
  locus <- unique(fb$Marker)

  stopifnot(length(locus) == 1)
  stopifnot(length(unique(fb$Plate)) == 1)


  # prepare columns to be used for calling/flagging of allele(s)
  fb$call <- NA
  fb$flag <- ""

  # 2. if below low run threshold (LRT), add flag "L"
  # Check for LRT from the get go. A calling proceeds like for any other run.
  LRT <- fetchTH(tbase, stat = "LRT", locus = locus)

  if (max(fb$Read_Count) <= LRT) {
    fb$flag <- paste(fb$flag, "L", sep = "")
  }

  # To process data row-wise, we can spit it first for convenience (because apply(x, MARGIN = 1, ...)
  # coerces input to vector which can hold only one type).
  x.split <- split(fb, f = 1:nrow(fb), drop = TRUE)

  cg <- function(x, maxA, fb, tbase, motif) {
    # exclude x from `fb` dataset so as not to compare against itself
    id <- rownames(x)
    fb <- fb[!(rownames(fb) %in% id), ]

    # prepare thresholds and find locus used down the line
    locus <- unique(fb$Marker)
    stopifnot(length(locus) == 1)
    IT <- fetchTH(tbase, stat = "IT", locus = locus)
    B <- fetchTH(tbase, stat = "B", locus = locus)

    # calculate relative size to maxA
    rs <- x$Read_Count/maxA

    # 1. If signal is below the ignore threshold (IT), remove/ignore A.
    if (rs <= IT) {
      return(NULL)
    }

    # 2. It makes sense to do this step before iterating over all alleles. See code in the
    # immediate definition of `callGenotype`.

    # 3. Find structural stutter for A. If one stutter found, call it, otherwise flag as "N".
    find.stt <- findStutter(x, fb = fb, motif = motif[motif$locus == locus, "motif"])

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

  run.A <- sapply(x.split, FUN = cg, maxA = max(x$Read_Count), fb = fb, tbase = tbase,
                  motif = motif, simplify = FALSE)

  # Remove all A which are flagged as non A
  run.A <- run.A[!sapply(run.A, FUN = is.null)]

  run.A <- do.call(rbind, run.A)

  # 5. if number of unflagged A > 2, add flag "M" to all alleles
  if (nrow(run.A) > 2) {
    run.A$flag <- paste(run.A$flag, "M", sep = "")
  }
  run.A
}

#' Find locus specific stutter.
#'
#' From a `fishbone` object, try to find a stutter.
#' @param x A `fishbone` object.
#' @param fb A `fishbone` object with cancidate stutter alleles.
#' @return A character string of the row name of the stutter(s).
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

  if (nrow(cand.stt) < 1) {
    return(FALSE)
  }

  found.stt <- cand.stt[new.stt == cand.stt$Sequence, ]

  # 3. return candidate stutter allele rowname
  if (nrow(found.stt) == 1) {
    message(sprintf("Stutter from allele %s", rownames(found.stt)))
    return(rownames(found.stt))
  } else {
    warning(sprintf("findStutter: Found multiple stutters for sample %s",
                    x$Sample_Name))
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
  if (length(out) != 1) stop(sprintf("Unable to fetch %s for locus %s.", stat, locus))
  out
}
