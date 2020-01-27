
snp2na <- function(x, snp, w) {
  if(length(w) != nrow(x))
    stop("dim mismatch")
  .Call('snp2na', PACKAGE = "gaston.utils", x@bed, snp, w)
}


