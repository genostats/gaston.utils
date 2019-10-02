
read.ms.file <- function(filename) {
  L <- .Call('read_ms_file', PACKAGE = "gaston.utils", filename)

  nrep <- length(L$reps);

  snp <- data.frame(chr = 1, id = sprintf("snp%0*d", nchar(L$snps), 1:L$snps), dist = NA_real_, pos = NA_integer_, 
                    A1 = '0', A2 = '1', rep = sprintf("rep%0*d", nchar(nrep), unlist(mapply(rep, 1:nrep, each = L$reps))),
                    stringsAsFactors = FALSE)
  ped <- data.frame(famid = 1:L$haps, id = 1:L$haps, father = 0, mother = 0, sex = 0, pheno = NA)

  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x)
  x
}

