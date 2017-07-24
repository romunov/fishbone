#' Find locus specific stutter given motif.
#'
#' From a `fishbone` object, try to find a stutter. It will try its best to find a stutter
#' no matter where in the sequence the motif may be missing.
#'
#' Algorithm:
#' 1. count number of motif repeats in a sequence
#' 2. count number of motif repeats in reverse of sequence
#' 3. if no motif repeats, signal that motif was not found
#' 4. find locations of all motif occurrences
#' 5. for each direction of sequence and for every location...
#' 5a. shorted candidate sequence
#' 5b. compare shortened sequence to possible stutter candidates
#' 5c. if found, signal stutter found and exit
#' 5d. if not found, remove second occurrence and compare to stutter candidates
#' 5e. repeat 3-3d until all positions tested
#' 6a. if match not found, declare no stutter found
#'
#' @param x A `fishbone` object.
#' @param fb A `fishbone` object with cancidate stutter alleles.
#' @return A character string of stutter(s) ID. If stutter is not found, it
#' returns NULL.
#' @importFrom stringi stri_locate_all_fixed stri_count stri_reverse stri_locate_all
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)

findStutter <- function(x, fb, motif) {
  # Create list of candidate stutters which are one repeat shorther
  # than the sequence in question.
  xl <- nchar(x$Sequence) # x length
  cand.st <- fb[fb$length == (xl - nchar(motif)), ] # list of candidate stutters

  if (nrow(cand.st) == 0) {
    # message(sprintf("Found no stutters for %s (ID: %s)", x$length, x$fbid))
    return(NULL)
  }

  # 1. count number of motif repeats in a sequence
  n.fwd <- stri_count(str = x$Sequence, fixed = motif)
  # 2. count number of motif repeats in reverse of sequence
  n.rev <- stri_count(str = stri_reverse(x$Sequence), fixed = motif)

  # 3. if no motif repeats, signal that motif was not found
  if (n.fwd == 0 & n.rev == 0) {
    # message(sprintf("Found no stutters for %s (ID: %s)", x$length, x$fbid))
    return(NULL)
  }

  # 5. for each direction of sequence and for every location...
  for (sq in list(x$Sequence, stri_reverse(x$Sequence))) {
    # 4. find locations of all motif occurrences
    sq.loc <- stri_locate_all(sq, fixed = motif)[[1]]

    for (i in 1:nrow(sq.loc)) {
      # browser()
      # 5a. shorted candidate sequence (paste together everything before/after found motif)
      newst <- paste(substr(sq, start = 1, stop = sq.loc[i, 1] - 1),
                     substr(sq, start = sq.loc[i, 2] + 1, xl),
                     sep = "")

      # 5b. compare shortened sequence to possible stutter candidates
      fnds <- cand.st[cand.st$Sequence == newst, ]

      # 5c. if found, signal stutter found and exit
      if (nrow(fnds) >= 1) {
        return(fnds$fbid)
      }
      # 5d. if not found, remove second occurrence and compare to stutter candidates
      if (nrow(fnds) == 0) {
        next
      }
      # 5e. repeat 3-3d until all positions tested
      return(NULL)
    }
  }
}
