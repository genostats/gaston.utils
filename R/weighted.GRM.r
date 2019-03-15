weighted.GRM <- function(x, weights = rep(1,ncol(x)), which.snps, autosome.only = TRUE, chunk = 1L) {
  if(missing(which.snps)) which.snps <- rep(TRUE, ncol(x))
  if(autosome.only) 
    which.snps <- which.snps & is.autosome(x@snps$chr)

  which.snps <- which.snps & (x@p > 0) & (x@p < 1) & (weights != 0)
  if(!is.logical(which.snps) | length(which.snps) != ncol(x))
    stop("which.snps must be a Logical vector of length ncol(x)")
  
  weights <- weights/sum(weights)
  w <- weights/(2*x@p*(1-x@p)) 
  K <- .Call('gg_weighted_Kinship_w', PACKAGE = "gaston.utils", x@bed, 2*x@p[which.snps], w[which.snps], which.snps, chunk) 

  if(!is.null(x@ped$id)) {
    if(anyDuplicated(x@ped$id) == 0) 
      rownames(K) <- colnames(K) <- x@ped$id
  }

  K
}
