#' @description Copied from skedastic::blus (GPL-3)
#' @source https://rdrr.io/github/tjfarrar/skedastic/src/R/blus.R


blus <- function(mainlm, omit = c("first", "last", "random"), keepNA = TRUE,
                 exhaust = NA, seed = 1234) {
  
  processmainlm(m = mainlm, needy = FALSE)
  
  n <- nrow(X)
  if (!is.na(seed)) set.seed(seed)
  omitfunc <- do_omit(omit, n, p, seed)
  Xmats <- do_Xmats(X, n, p, omitfunc$omit_ind)
  singular_matrix <- FALSE
  
  if (is.singular.mat(Xmats$X_ord_sq) ||
      is.singular.mat(Xmats$X0)) {
    
    singular_matrix <- TRUE
    message("Passed `omit` argument resulted in singular matrix; BLUS residuals
          could not be computed. Randomly chosen combinations of indices to
          omit will be attempted according to `exhaust` argument passed.")
    ncombn <- choose(n, p)
    if ((is.na(exhaust) || is.null(exhaust))) {
      dosample <- (ncombn > 1e4)
      numsample <- 1e4
    } else if (exhaust <= 0) {
      stop("`exhaust` is not positive; no attempts will be made to find subset to omit")
    } else {
      dosample <- (ncombn > exhaust)
      numsample <- min(ncombn, exhaust)
    }
    
    if (dosample) {
      subsetstotry <- unique(t(replicate(numsample, sort(sample(x = n, size = p)))))
      rowstodo <- 1:nrow(subsetstotry)
    } else {
      subsetstotry <- t(utils::combn(n, p))
      maxrow <- min(exhaust, nrow(subsetstotry), na.rm = TRUE)
      rowstodo <- sample(1:nrow(subsetstotry), maxrow, replace = FALSE)
    }
    
    for (r in rowstodo) {
      omitfunc <- do_omit(subsetstotry[rowstodo[r], , drop = FALSE], n, p)
      Xmats <- do_Xmats(X, n, p, omitfunc$omit_ind)
      if (!is.singular.mat(Xmats$X_ord_sq) &&
          !is.singular.mat(Xmats$X0)) {
        singular_matrix <- FALSE
        message(paste0("Success! Subset of indices found that does not yield singular
                   matrix: ", paste(omitfunc$omit_ind, collapse = ",")))
        break
      }
    }
  }
  if (singular_matrix) stop("No subset of indices to omit was found that
                            avoided a singular matrix.")
  
  keep_ind <- setdiff(1:n, omitfunc$omit_ind)
  G <- Xmats$X0 %*% solve(Xmats$X_ord_sq) %*% t(Xmats$X0)
  Geig <- eigen(G, symmetric = TRUE)
  lambda <- sqrt(Geig$values)
  q <- as.data.frame(Geig$vectors)
  Z <- Reduce(`+`, mapply(`*`, lambda / (1 + lambda),
                          lapply(q, function(x) tcrossprod(x)),
                          SIMPLIFY = FALSE))
  e0 <- e[omitfunc$omit_ind]
  e1 <- e[keep_ind]
  e_tilde <- c(e1 - Xmats$X1 %*% solve(Xmats$X0) %*% Z %*% e0)
  
  if (keepNA) {
    rval <- rep(NA_real_, n)
    rval[keep_ind] <- e_tilde
  } else {
    rval <- e_tilde
  }
  rval
}