#' @description Copied from skedastic::anscombe (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/anscombe.R


anscombe <- function(mainlm, studentise = TRUE, statonly = FALSE) {
  
  processmainlm(m = mainlm, needyhat = TRUE)
  
  n <- nrow(X)
  M <- fastM(X, n)
  
  s_sq <- sum(e ^ 2) / (n - p)
  tbar <- sum(diag(M) * yhat) / (n - p)
  if (studentise) {
    method <- "Anscombe-Bickel (studentised)"
    sigma_tilde_sq_xsq <- sum((yhat - tbar) ^ 2) * sum((e ^ 2 -
                                                          mean(e ^ 2)) ^ 2) / (n - p)
    teststat <- sum(e ^ 2 * (yhat - tbar)) / sqrt(sigma_tilde_sq_xsq)
  } else {
    method <- "Anscombe (non-studentised)"
    sigma_tilde_sq <- 2 * (n - p) / (n - p + 2) * s_sq ^ 2 *
      t(yhat - tbar) %*% (M ^ 2) %*% (yhat - tbar)
    teststat <- sum(e ^ 2 * (yhat - tbar)) / sqrt(sigma_tilde_sq)
  }
  
  if (statonly) return(teststat)
  
  pval <- 2 * stats::pnorm(abs(teststat), lower.tail = FALSE)
  
  rval <- structure(list(statistic = teststat, parameter = NULL,
                         p.value = pval,
                         null.value = 0,
                         alternative = "two.sided", method = method), class = "htest")
  broom::tidy(rval)
}