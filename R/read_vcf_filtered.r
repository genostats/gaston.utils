read.vcf.filtered <- function(file, x, by = "chr:pos:alleles", chr.ids, samples, get.info = FALSE, haplotypes = FALSE, verbose = getOption("gaston.verbose",TRUE)) {
  filename <- path.expand(file)

  if(class(x) == "bed.matrix")
    x <- x@snps
  if(class(x) != "data.frame")
    stop("x should be a bed.matrix or a data frames")

  if(missing(chr.ids))
    set.chr.ids() 
  else
    set.chr.ids(chr.ids)

  if(missing(samples)) 
    samples <- character(0)

  if(haplotypes) {
    if(by == "chr:pos") 
      L <- .Call("gg_read_vcf_chr_pos_haplo", PACKAGE = "gaston.utils", filename, get.info, x$chr, x$pos, samples)
    else if(by == "chr:pos:alleles")
      L <- .Call("gg_read_vcf_chr_pos_al_haplo", PACKAGE = "gaston.utils", filename, get.info, x$chr, x$pos, x$A1, x$A2, samples)
    else
      stop("Parameter 'by' should be 'chr:pos' or 'chr:pos:alleles'")
    famid <- rep(L$sample, each = 2)
    id <- paste0( rep(L$sample, each = 2), c(".1", ".2"))
  } else {
    if(by == "chr:pos") 
      L <- .Call("gg_read_vcf_chr_pos", PACKAGE = "gaston.utils", filename, get.info, x$chr, x$pos, samples)
    else if(by == "chr:pos:alleles")
      L <- .Call("gg_read_vcf_chr_pos_al", PACKAGE = "gaston.utils", filename, get.info, x$chr, x$pos, x$A1, x$A2, samples)
    else
      stop("Parameter 'by' should be 'chr:pos' or 'chr:pos:alleles'")
    id <- famid <- L$samples
  }

  snp <- data.frame(chr = L$chr, id = L$id, dist = rep(0, length(L$chr)), pos = L$pos , A1 = L$A1, A2 = L$A2,
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)

  if(get.info) {
    w <- substr(names(L),1,4) == "info"
    snp[ ,names(L)[w] ] <- L[w]
  }

  ped <- data.frame(famid = famid, id = id, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)
  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}
