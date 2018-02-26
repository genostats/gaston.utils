
snplist_file <- "/WORKING_DIRECTORY/hperdry/LDSCORE/eur_w_ld_chr/w_hm3.snplist"
ldscores <- "/WORKING_DIRECTORY/hperdry/LDSCORE/eur_w_ld_chr/"

add.ld.scores <- function(snp.table, snp.list.file, ld.score.dir, ld.score.files = sprintf("%s/%d.l2.ldscore.gz", ld.score.dir, 1:22)) {
  T <- snp.table

  noms <- names(T)
  id  <- if("id" %in% noms) "id"  else "SNP"
  #pos <- if("id" %in% noms) "pos" else "BP"

  if(!all(c(id, "A1", "A2") %in% noms))
    stop("Columns id or SNP, A1, and A2 are mandatory")

  # +/- le travail de munge.py
  # on lit dans le fichier snp list ceux qui ont le bon id... 
  cat(nrow(T), "SNP given\n")
  cat("Reading SNP list reference file", snp.list.file, "\n")
  SnpList <- read.snp.list(snp.list.file, T[,id])
  w <- match(T[,id],  SnpList$SNP)
  T$in_ref <- !is.na(w)

  cat(sum(T$in_ref), "SNPs of the table are in reference\n")
  SnpList <- SnpList[w, ]

  # vérif alleles
  L <- check.alleles(T$A1, T$A2, SnpList$A1, SnpList$A2)
  cat(sum(L$OK), "SNPs alleles can be matched with alleles from reference \n(")
  cat(sum(L$swap), "swaps, ", sum(L$flip), "strand flips)\n")

  T$allele_ok <- L$OK
  T$swap <- L$swap
  T$flip <- L$flip

  # T$Z <- T$beta/T$sd
  # T$Z <- ifelse( L$swap, -T$Z, T$Z)

  # le début du travail de ldsc
  cat("Readind LD scores\n")
  S <- read.ld.score(ld.score.files, T[T$in_ref, id] )    
  cat(nrow(S$ld.scores), "of SNPs in reference have LD scores (total available LD scores = ", S$total.nb.snps, ")\n")
  Scores <- S$ld.scores

  w <- match(T[,id], Scores$SNP)
  T$ld.score <- Scores$L2[w]

  T
}

