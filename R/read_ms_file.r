
read.ms.file <- function(filename, rep = 1L) {
  L <- .Call('read_ms_file', PACKAGE = "gaston.utils", filename, rep)

  snp <- data.frame(chr = 1, id = sprintf("Ga%0*d", nchar(L$snps), 1:L$snps), dist = NA_real_, pos = as.integer(1e5*L$pos), 
                    A1 = '0', A2 = '1', stringsAsFactors = FALSE)
  ped <- data.frame(famid = 1:L$inds, id = 1:L$inds, father = 0, mother = 0, sex = 0, pheno = NA)

  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x)
  x
}

read.ms.file.haplotypes <- function(filename, rep = 1L) {
  L <- .Call('read_ms_file_haplotype', PACKAGE = "gaston.utils", filename, rep)

  snp <- data.frame(chr = 1, id = sprintf("Ga%0*d", nchar(L$snps), 1:L$snps), dist = NA_real_, pos = as.integer(1e5*L$pos), 
                    A1 = '0', A2 = '1', stringsAsFactors = FALSE)
  ped <- data.frame(famid = 1:L$haps, id = 1:L$haps, father = 0, mother = 0, sex = 0, pheno = NA)

  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x)
  x
}
