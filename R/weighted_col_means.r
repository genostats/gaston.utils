weighted.col.means <- function(x, alpha) {
  .Call('gg_weighted_col_means', PACKAGE = "gaston.utils", x@bed, alpha)
}
