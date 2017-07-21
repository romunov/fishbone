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
#' 1. pogledamo, če je max A višji od L
#' 1a. če ni, flagaj vse "L"
#' 2. najdi najvišji alel
#' 3. preveri, če ima stutter
#' 3a. če ima stutter, je to kandidatni alel
#' 3b. če ima stutter, pa je prenizek (D), mu daj flag "D"
#' 3c. stutterju dodaj flag $stutter = TRUE (rezultat sta alel in stutter)
#' 3d. če nima stutterja, preveri absolutno višino
#' 3da. če je ( > AlleleWithNoStutterHeight) mu daš flag "N"
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

  stopifnot(length(locus) == 1)
  stopifnot(length(unique(fb$Plate)) == 1)

  # prepare columns to be used for calling/flagging of allele(s)
  fb$call <- ""
  fb$flag <- ""

  # prepare statistics
  L <- fetchTH(tbase, stat = "LowCount", locus = locus)
  R <- fetchTB(tbase, stat = "RelativeLow", locus = locus)

  # 1. remove all alleles which are very low (below [locus specific]% of max read)
  rel.2.max <- fb$Read_Count/max(fb$Read_Count)
  fb <- fb[rel.2.max > R, ]

  # 2. if below low run threshold (L), add flag "L"
  if (max(fb$Read_Count) <= L) {
    fb$flag <- paste(fb$flag, "L", sep = "")
  }

  # To process data row-wise, we can spit it first for convenience (because apply(x, MARGIN = 1, ...)
  # coerces input to vector which can hold only one type).
  x.split <- split(fb, f = 1:nrow(fb), drop = TRUE)

  # Call this function on every allele.
  cg <- function(x, maxA, fb, tbase, motif) {
    # prepare thresholds and find locus
    locus <- unique(fb$Marker)
    stopifnot(length(locus) == 1)
    IT <- fetchTH(tbase, stat = "RelativeLow", locus = locus)
    B <- fetchTH(tbase, stat = "Disbalance", locus = locus)

    # exclude x from `fb` dataset so as not to compare against itself
    id <- rownames(x)
    remove.self <- !(rownames(fb) %in% id)
    fb.noself <- fb[remove.self, ]

    # If there is only one candidate allele, it can't have a stutter, so add a flag.
    if (length(fb.noself) == 0) {
      fb$flag <- paste(fb$flag, "N")
      return(NULL)
    } else {
      fb <- fb.noself
    }

    # calculate relative size to maxA
    rs <- x$Read_Count/maxA

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
  run.A <- sapply(x.split, FUN = cg, maxA = max(fb$Read_Count), fb = fb, tbase = tbase,
                  motif = motif, simplify = FALSE)

  # Remove all A which are flagged as non A
  run.A <- run.A[!sapply(run.A, FUN = is.null)]

  run.A <- do.call(rbind, run.A)

  if (length(run.A) == 0) {
    message(sprintf("No genotypes called for %s (%s)", unique(fb$Sample_Name), unique(fb$Marker)))
    return(NA)
  }

  # 5. if number of unflagged A > 2, add flag "M" to all alleles

  if (sum(sapply(run.A$flag, nchar) == "") > 2) {
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
  if (nrow(found.stt) == 0) {
    message(sprintf("No stutter from allele %s", x$Allele))
    return(NA)
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
