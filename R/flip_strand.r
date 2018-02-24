flip.strand <- function(allele) {
  .Call('gg_flip_strand', PACKAGE = "gaston.utils", allele)
}
