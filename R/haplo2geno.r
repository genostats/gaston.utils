
haplo2geno <- function(x, H1, H2) {
  if(any(H1 > ncol(x)) | any(H1 > ncol(x)) | any(H1 <= 0) | any(H2 <= 0)) 
    stop("Wrong haplotype numbers")

  bed <- .Call('haplo2geno', PACKAGE = "gaston.utils", x@bed, H1, H2)

  ped <- data.frame(famid = 1:length(H1), id = 1:length(H2), father = 0, mother = 0, sex = 0, pheno = NA)

  x <- new("bed.matrix", bed = bed, snps = x@snps[,1:6], ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x)
  x
}
