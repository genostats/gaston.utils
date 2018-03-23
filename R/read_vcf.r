read.vcf.1 <- function(file, chr, range, get.info = FALSE, convert.chr = TRUE, verbose = getOption("gaston.verbose",TRUE)) {
  filename <- path.expand(file)

  if(missing(chr)) {
    chr <- ""
    low <- high <- 0;
  } else {
#    chr <- as.character(chr)
    if(missing(range))
      low <- high <- -1
    else {
      low <- range[1]
      high <- range[2]
      if(low > high) 
        stop("Bad range")
    }
  }

  set.chr.ids()
  L <- .Call("gg_read_vcf_chr_range", PACKAGE = "gaston.utils", filename, get.info, chr, low, high)

  snp <- data.frame(chr = L$chr, id = L$id, dist = rep(0, length(L$chr)), pos = L$pos , A1 = L$A1, A2 = L$A2,
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)

  if(get.info) {
    w <- substr(names(L),1,4) == "info"
    snp[ ,names(L)[w] ] <- L[w]
  }

  if(convert.chr) {
    chr <- as.integer(L$chr)
    chr[L$chr == "X"  | L$chr == "x"]  <- getOption("gaston.chr.x")[1]
    chr[L$chr == "Y"  | L$chr == "y"]  <- getOption("gaston.chr.y")[1]
    chr[L$chr == "MT" | L$chr == "mt"] <- getOption("gaston.chr.mt")[1]
    if(any(is.na(chr)))
      warning("Some unknown chromosomes id's (try to set convert.chr = FALSE)")
    snp$chr <- chr
  }

  ped <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}
