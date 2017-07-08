#' Coerces a data.frame (like object) to a `fishbone` object.
#'
#' Expected columns are Sample_Name, Plate, Marker, Read_Count, Allele, Sequence.
#'
#' @param x Object of class \code{data.frame} or similar.
#' @param motif.length Integer. Length of the repeating motif, locus specific.
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)
#' @export
as.fishbone <- function(x, motif.length) {
  if (!all(c("Sample_Name", "Plate", "Marker", "Read_Count", "Allele", "Sequence") %in% names(x))) {
    stop("Please see the documentation for mandatory column names.")
  }
  find.sequence <- which.max(sapply(x[1,], nchar))
  x <- as.data.frame(x)
  x$seq.preview <- sapply(x[find.sequence], substr, start = 1, stop = 5)
  x$seq.preview <- paste(x$seq.preview, "...", sep = "")

  x$lengths <- as.numeric(gsub("(^\\d+)_(\\d+)$", "\\1", x$Allele))
  x$poly <- as.numeric(gsub("(^\\d+)_(\\d+)$", "\\2", x$Allele))
  # remove.ln.poly <- which(names(x) %in% c("lengths", "poly"))
  # x <- x[order(x$lengths, rev(x$poly), decreasing = TRUE), -remove.ln.poly]
  x <- x[order(x$lengths, rev(x$poly), decreasing = TRUE), ]

  attr(x, "sequence") <- find.sequence
  attr(x, "motif.length") <- motif.length
  class(x) <- c("fishbone", class(x))
  x
}
