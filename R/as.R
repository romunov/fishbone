#' Coerces a data.frame (like object) to a `fishbone` object.
#'
#' Expected columns are Sample_Name, Plate, Marker, Read_Count, Allele, Sequence.
#'
#' @param x Object of class \code{data.frame} or similar.
#' @param motif Integer or data.frame or character. Length of the repeating motif, locus specific.
#' Motif (raw, e.g. "ctat") can also be provided in the form of a character. If a data.frame provided, it
#' should have two columns:
#' locus (exact locus name, e.g. "03")
#' motif (motif, eg.. "ctat")
#'
#' @return It adds an extra column with shortened sequence to enable pretty printing. It also
#' adds two more columns, namely `length` which reads sequence length and `poly` which is the
#' variation of the allele.
#'
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)
#' @export

as.fishbone <- function(x, motif) {
  if (!all(c("Sample_Name", "Plate", "Marker", "Read_Count", "Allele", "Sequence") %in% names(x))) {
    stop("Please see the documentation for mandatory column names.")
  }

  # In case `motif.length` is a data.frame, fetch the motif for the appropriate locus
  # and find the length of the repeat.
  if (is.data.frame(motif)) {
    motif <- nchar(motif[motif$locus == unique(x$Marker), "motif"])
  }

  # If raw motif is provided, count motif length.
  if (is.character(motif)) {
    motif <- nchar(motif)
  }

  x$lengths <- as.numeric(gsub("(^\\d+)_(\\d+)$", "\\1", x$Allele))
  x <- x[order(x$lengths, decreasing = TRUE), ]

  x <- as.data.table(x)
  attr(x, "motif.length") <- motif
  class(x) <- c("fishbone", class(x))
  x
}
