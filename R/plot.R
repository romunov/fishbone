#' Plot a list of \code{fishbone} objects.
#'
#' @param x A list of \code{fishbone} objects.
#'
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)
#' @export
#'
#' @importFrom ggplot2 ggplot
#' @importFrom gridExtra grid.arrange

plotRuns <- function(x) {
  xy.plots <- vector("list", length(x))

  for (i in 1:length(xy.plots)) {
    xy.plots[[i]] <- plot.fishbone(x[[i]])
  }

  grid.arrange(grobs = xy.plots, ncol = 4)
}

#' Plot fishbone object.
#'
#' Plots candidate alleles and colors those that are not expected given the
#' nucleotide repeat length. For instance, if a difference in length between two
#' alleles is not say 4 for tetranucleotides, it could be due to tag jump and the
#' allele is probably not what you're after.
#'
#' @param x \code{fishbone} object.
#' @param ... Currently not used.
#'
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)
#'
#' @export
#' @importFrom ggplot2 ggplot
#' @importFrom ggplot2 aes
#' @importFrom ggplot2 theme_bw
#' @importFrom ggplot2 geom_col
#' @importFrom ggplot2 scale_fill_brewer
#'
plot.fishbone <- function(x, ...) {
  stopifnot(any(class(x) %in% "fishbone"))

  # Find which allele is not a slip of length of the nucleotide repeat sequence.
  # For example, if you have a tetranucleotide, values should be a multiple of 4.
  rl <- attr(x, "motif.length")
  canal <- x$lengths[which.max(x$Read_Count)] # candidate allele with highest count
  x$group <- "ok"
  # find if lenght of others corresponds to the multiple of nucleotide length and
  # mark those as not OK
  x$group[((x$lengths - canal) %% rl) != 0] <- "not OK"

  # Coerce to factor to preserve correct ordering
  x$Allele <- factor(x$Allele, levels = rev(x$Allele))

  out <- ggplot(x, aes_string(x = "Allele", y = "Read_Count", fill = "group")) +
    theme_bw() +
    ylab("") +
    xlab(sprintf("Alleles for locus %s", unique(x$Marker))) +
    # assign blue to ok and red to not OK
    scale_fill_manual(values = c("not OK" = "#E41A1C", "ok" = "#377EB8"), guide = FALSE) +
    geom_col()
  out
}
