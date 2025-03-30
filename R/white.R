#' @description Copied from skedastic::white (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/white.R


white <- function(mainlm, interactions = FALSE, statonly = FALSE) {
  
  processmainlm(m = mainlm, needy = FALSE)
  
  hasintercept <- columnof1s(X)
  if (!hasintercept[[1]]) {
    message("Intercept included in auxiliary design matrix")
  } else {
    X <- X[, -hasintercept[[2]], drop = FALSE]
  }
  
  n <- nrow(X)
  
  Z <- cbind(1, X, X ^ 2)
  if (interactions) {
    Z <- cbind(Z, generate_interactions(X))
  }
  
  esq <- e ^ 2
  auxlm <- stats::lm.fit(Z, esq)
  iota <- rep(1, n)
  N <- diag(n) - 1 / n * (tcrossprod(iota))
  e_aux <- auxlm$residuals
  teststat <- as.double(n * (1 - crossprod(e_aux) / (t(esq) %*% N %*% esq)))
  if (statonly) return(teststat)
  
  df <- ncol(Z) - 1
  pval <- stats::pchisq(teststat, df = df, lower.tail = FALSE)
  
  rval <- structure(list(statistic = teststat, parameter = df, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater",
                         method = "White's Test"), class = "htest")
  broom::tidy(rval)
}

white_lm <- white