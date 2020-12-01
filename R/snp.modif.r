
snp.modif <- function(x, snp, new.vals) {
  if(length(new.vals) != nrow(x))
    stop("dim mismatch")
  .Call('snp_modif_in_place', PACKAGE = "gaston.utils", x@bed, snp, new.vals)
}


