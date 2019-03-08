read.vcf.1000genomes <- function(dir, x, by = "chr:pos:alleles", samples, populations, super.populations, phase = 3, get.info = FALSE, verbose = getOption("gaston.verbose",TRUE)) {

  # pour l'instant on ne lit que les autosomes !! (c'est le bordel avec le chr Y)
  if(phase == 3) {
    filename <- sprintf("%s/ALL.chr%d.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz", path.expand(dir), 1:22)
  } else if(phase == 1) {
    filename <- sprintf("%s/ALL.chr%d.SHAPEIT2_integrated_phase1_v3.20101123.snps_indels_svs.genotypes.all.vcf.gz", path.expand(dir), 1:22)
  } else {
    stop("Parameter 'phase' should be 1 or 3 (default)")
  }

  if(!missing(x)) {
    filter.snps <- TRUE
    if(class(x) == "bed.matrix")
      x <- x@snps
    if(class(x) != "data.frame")
      stop("x should be a bed.matrix or a data frames")
  } else {
    filter.snps <- FALSE
  }

  # pour 1000 genomes le défaut doit être ok
  set.chr.ids() 

  if(missing(samples)) {
    if(!missing(populations))
      samples <- KG.samples$sample[ KG.samples$population %in% populations ]
    else if(!missing(super.populations))
      samples <- KG.samples$sample[ KG.samples$super.population %in% super.populations ]
    else
      samples <- character(0)
  }

  if(filter.snps) {
    if(by == "chr:pos") 
      L <- .Call("gg_read_vcf_chr_pos", PACKAGE = "gaston.utils", filename, get.info, x$chr, x$pos, samples)
    else if(by == "chr:pos:alleles")
      L <- .Call("gg_read_vcf_chr_pos_al", PACKAGE = "gaston.utils", filename, get.info, x$chr, x$pos, x$A1, x$A2, samples)
    else
      stop("Parameter 'by' should be 'chr:pos' or 'chr:pos:alleles'")
   } else {
   L <- .Call("gg_read_vcf_chr_range", PACKAGE = "gaston.utils", filename, get.info, -1L, 0L, 0L, samples)
  }

  snp <- data.frame(chr = L$chr, id = L$id, dist = rep(0, length(L$chr)), pos = L$pos , A1 = L$A1, A2 = L$A2,
                    quality = L$quality, filter = factor(L$filter), stringsAsFactors = FALSE)

  if(get.info) {
    w <- substr(names(L),1,4) == "info"
    snp[ ,names(L)[w] ] <- L[w]
  }

  ped <- data.frame(famid = L$samples, id = L$samples, father = 0, mother = 0, sex = 0, pheno = NA, stringsAsFactors = FALSE)

  # ajoute population / super.population
  m <- match(ped$id, KG.samples$sample)
  ped$population       <- droplevels( KG.samples$population[m] )
  ped$super.population <- droplevels( KG.samples$super.population[m] )

  x <- new("bed.matrix", bed = L$bed, snps = snp, ped = ped,
           p = NULL, mu = NULL, sigma = NULL, standardize_p = FALSE,
           standardize_mu_sigma = FALSE )

  if(getOption("gaston.auto.set.stats", TRUE)) x <- set.stats(x, verbose = verbose)
  x
}
