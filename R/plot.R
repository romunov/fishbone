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
  if (all(is.na(x))) {
    message("Nothing to draw.")
    return(NA)
  }

  xy.plots <- vector("list", length(x))

  for (i in 1:length(xy.plots)) {
    xy.plots[[i]] <- plot.fishbone(x[[i]])
  }

  out <- do.call(arrangeGrob, c(xy.plots, list(ncol = 4)))
  out
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

  # Coerce to factor to preserve correct ordering
  x$Allele <- factor(x$Allele, levels = rev(x$Allele))
  x$lengths <- factor(x$lengths, levels = rev(unique(x$lengths)))

  out <- ggplot(x, aes_string(x = "lengths", y = "Read_Count", group = "Allele")) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.25),
          axis.title.x = element_text(size = 10)) +
    ylab("") +
    geom_col(position = "dodge", width = 0.5) +
    geom_text(aes_string(x = "lengths", y = "Read_Count", label = "Allele"), position = position_dodge(0.9),
              vjust = -0.5) +
    xlab(sprintf("Marker: %s; Sample: %s; Run: %s", unique(x$Marker), unique(x$Sample_Name), unique(x$Plate)))
  out
}
