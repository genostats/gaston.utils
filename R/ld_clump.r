
LD.clump <- function(x, p, r2.threshold, p.threshold, max.dist = 500e3) {
  if(!is(x, "bed.matrix")) 
    stop("x is not a bed matrix")
  if(is.data.frame(p)) {
    p <- p$p
  } 
  if(length(p) != ncol(x))
    stop("Dimensions mismatch")

  or <- order(p) - 1 # -1 pour des indices qui démarrent à 0...

  I <- .Call("ld_clump", PACKAGE = "gaston.utils", x@bed, x@mu, x@sigma, r2.threshold, x@snps$pos, x@snps$chr, max.dist, or);
  data.frame( chr = x@snps$chr, id = x@snps$id, pos = x@snps$pos, p = p, index = I+1)
}


