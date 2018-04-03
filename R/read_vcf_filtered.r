read.vcf.filtered <- function(file, x, by = "chr:pos:alleles", chr.ids, samples, get.info = FALSE, verbose = getOption("gaston.verbose",TRUE)) {
  filename <- path.expand(file)

  if(missing(chr)) {
    chr <- -1L 
    low <- high <- -1L;
  } else {
    if(missing(range)) {
      low <- high <- -1L;
    } else {
      low  <- range[1]
      high <- range[2]
      if(low > high) 
        stop("Bad range")
    }
  }

  if(missing(chr.ids))
    set.chr.ids() 
  else
    set.chr.ids(chr.ids)

  if(missing(samples)) 
    samples <- character(0)

  L <- .Call("gg_read_vcf_chr_range", PACKAGE = "gaston.utils", filename, get.info, chr, low, high, samples)

  snp <- data.frame(chr = L$chr, id = L$id, dist = rep(0, length(L$chr)), pos = L$pos , A1 = L$A1, A2 = L$A2,
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)

  if(get.info) {
    w <- substr(names(L),1,4) == "info"
    snp[ ,names(L)[w] ] <- L[w]
  }

  ped <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}
