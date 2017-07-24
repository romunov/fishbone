#' Fetch threshold values from some database.
#' @param x A slice of `fishbone` object.
#' @param stat Which statistic to fetch, e.g. "RelativeLow", "Disbalance", "Stutter", "LowCount".
#' @param locus Character. For which locus to fetch statistics?
#' @author Roman Lustrik (roman.lustri@@biolitika.si)

fetchTH <- function(x, stat, locus) {
  out <- x[x$Marker %in% locus, ..stat]
  if (length(out) != 1) stop(sprintf("Unable to fetch %s for locus %s.", stat, locus))
  out
}
