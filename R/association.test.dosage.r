association.test.dosage <- function(filename, Y, X, method = c("lm", "lmm"), response = c("quantitative", "binary"), 
                                    test = c("score", "wald", "lrt"), K, eigenK, beg, end, p = 0, 
                                    tol = .Machine$double.eps^0.25, ...) {
  filename <- path.expand(filename)
  dims <- dim.dosage.file(filename)
  nb.inds <- dims[1]
  nb.snps <- dims[2]

  if(missing(beg)) beg <- 1
  if(missing(end)) end <- nb.snps
  
  if(beg < 1 || end > nb.snps) stop("range too wide")
  if(length(Y) != nb.inds) stop("Dimension of Y and #individuals in ", filename, " mismatch")
 
  if(missing(X)) X <- rep(1, nb.inds); # default = intercept 
  X <- as.matrix(X)

  if(nrow(X) != nb.inds) stop("Dimensions of Y and #individuals in ", filename, " mismatch")

  # check dimensions before anything
  if(!missing(K)) {
    if(nb.inds != nrow(K) | nb.inds != ncol(K)) stop("K dimensions and #individuals in ", filename, " mismatch")
  }
  if(!missing(eigenK)) {
    if(nb.inds != nrow(eigenK$vectors) | nb.inds != ncol(eigenK$vectors) | nb.inds != length(eigenK$values)) 
      stop("eigenK dimensions and #individuals in ", filename, " mismatch")
  }


  response <- match.arg(response)
  test <- match.arg(test)
  method <- match.arg(method)

  # preparation de X 
  if(p > 0) {
    if((method == "lmm" & response == "quantitative" & test == "score") |
       (method == "lmm" & response == "binary") |
       (method == "lm")) { # il faut ajouter les PCs à X
      X <- cbind(X, eigenK$vectors[,seq_len(p)])
      X <- gaston:::trans.X(X, mean.y = mean(Y))
    } else {
      X <- gaston:::trans.X(X, eigenK$vectors[,seq_len(p)], mean(Y))
    }
  } else {
    X <- gaston:::trans.X(X, mean.y = mean(Y))
  }
  
  # random effect
  if(method == "lmm") { 

    # if(response == "binary" & test != "score") {
    #  warning('Binary phenotype and method = "lmm" force test = "score"')
    #  test <- "score"
    # }

    if(test == "score" | response == "binary") {
      if(missing(K)) stop("For a score test and for binary traits, argument K is mandatory")
      # avec le score test on peut gérer les données manquantes dans Y
      if( any(is.na(Y)) ) {
        w <- !is.na(Y)
        X <- as.matrix(X[w,])
        Y <- Y[w]
        if (is.list(K)) K <- lapply(K, function(x) x[w,w]) else K <- K[w,w]
        warning(sum(!w), 'individuals with missing phenotype are ignored.\n')
      } 
    } else {
      if(missing(eigenK)) 
        stop("For quantitative Wald and LRT tests, argument eigenK is mandatory")
      if( any(is.na(Y)) ) 
        stop("Can't handle missing data in Y, please recompute eigenK for the individuals with non-missing phenotype")
    }

    if(response == "quantitative") { # score (argument K), wald ou lrt (eigen K) possibles
      if(test == "score") {
        model <- lmm.aireml(Y, X = X, K, get.P = TRUE, ... )
        t <- .Call("gg_GWAS_lmm_score_dosages", PACKAGE = "gaston.utils", filename, model$Py, model$P, beg, end)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE)
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        t <- .Call("GWAS_dosage_lmm_wald", PACKAGE = "gaston.utils", filename, Y, X, p, eigenK$values, eigenK$vectors, beg, end, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
      } else { # test == "lrt"
        X <- cbind(X, 0) # space for the SNP
        t <- .Call("GWAS_dosage_lmm_lrt", PACKAGE = "gaston.utils", filename, Y, X, p, eigenK$values, eigenK$vectors, beg, end, tol)
        t$p <- pchisq( t$LRT, df = 1, lower.tail=FALSE)
      }
    } else { # response == "binary", seulement le score test, avec argument K
      if(test == "score") {
        model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
        omega <- model$BLUP_omega
        if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta
        pi <- 1/(1+exp(-omega))
        t <- .Call("gg_GWAS_lmm_score_dosages", PACKAGE = "gaston.utils", filename, Y-pi, model$P, beg, end)
        t$p <- pchisq( t$score, df = 1, lower.tail=FALSE) 
      } else if(test == "wald") {
        X <- cbind(X, 0) # space for the SNP
        t <- .Call("GWAS_dosage_logitmm_wald_f", PACKAGE = "gaston.utils", filename, Y, X, K, beg, end, tol)
        t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
      } else stop("LRT test for binary trait not available")
    }
  }

  # only fixed effects
  if(method == "lm") {
    if(test != "wald") warning('Method = "lm" force test = "wald"')
    if( any(is.na(Y)) ) 
      stop("Can't handle missing data in Y")
    if(response == "quantitative") {
      t <- .Call("GWAS_dosage_lm_quanti", PACKAGE = "gaston.utils", filename, Y, Q, beg, end);
      t$p <- pt( abs(t$beta/t$sd), df = nb.inds - ncol(Q) - 1, lower.tail=FALSE)*2
    }
    if(response == "binary") {
      X <- cbind(X,0)
      t <- .Call("GWAS_dosage_logit_wald_f", PACKAGE = "gaston.utils", filename, Y, X, beg, end, tol);
      t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail=FALSE)
    }
  }

  data.frame( t )
}


