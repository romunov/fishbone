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