
read.ms.file <- function(filename, rep = 1L) {
  L <- .Call('read_ms_file', PACKAGE = "gaston.utils", filename, rep)

  snp <- data.frame(chr = 1, id = paste("id", 1:L$snps), dist = L$pos, pos = NA_integer_, A1 = '0', A2 = '1')
  ped <- data.frame(famid = rep(NA_character_, L$inds), id = NA_character_, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)

  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x

}
