#' Improves printing of \code{fishbone} objects.
#'
#' To prevent cluttering of the output, this is a convencience functionality which
#' truncates Sequence to a few characters. Coerce to \code{data.frame} for the original
#' values.
#'
#' @param x Object of class \code{fishbone}.
#' @param ... Currently not used.
#'
#' @importFrom data.table as.data.table
#'
#' @author Roman Lustrik (roman.lustrik@@biolitika.si)
#' @export
# print.fishbone <- function(x, ...) {
#   print(head(as.data.table(x)))
# }

head.fishbone <- function(x, ...) {
  print(head(as.data.table(x), ...))
}