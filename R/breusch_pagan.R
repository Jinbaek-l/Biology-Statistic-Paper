#' @description Copied from skedastic::breusch_pagan (GPL-3)
#' @source https://rdrr.io/cran/skedastic/src/R/breusch_pagan.R


breusch_pagan <- function(mainlm, auxdesign = NA, koenker = TRUE,
                          statonly = FALSE) {
  
  auxfitvals <- ifelse(all(is.na(auxdesign)) | is.null(auxdesign), FALSE,
                       auxdesign == "fitted.values")
  processmainlm(m = mainlm, needy = auxfitvals, needyhat = auxfitvals,
                needp = FALSE)
  
  if (all(is.na(auxdesign)) || is.null(auxdesign)) {
    Z <- X
  } else if (is.character(auxdesign)) {
    if (auxdesign == "fitted.values") {
      Z <- t(t(yhat))
    } else stop("Invalid character value for `auxdesign`")
  } else {
    Z <- auxdesign
    if (nrow(auxdesign) != nrow(X)) stop("No. of observations in `auxdesign`
                                         must match\nno. of observations in
                                         original model.")
  }
  
  hasintercept <- columnof1s(Z)
  if (!hasintercept[[1]]) {
    Z <- cbind(1, Z)
    message("Column of 1's added to `auxdesign`")
  }
  
  q <- ncol(Z) - 1
  n <- nrow(Z)
  sigma_barsq <- sum(e ^ 2) / n
  w_hat <- e ^ 2 - sigma_barsq
  
  if (koenker) {
    method <- "Koenker (studentised)"
    teststat <- n * sum(stats::lm.fit(Z, w_hat)$fitted.values ^ 2) / sum(w_hat ^ 2)
  } else {
    method <- "Breusch-Pagan (non-studentised)"
    teststat <- sum(stats::lm.fit(Z, w_hat)$fitted.values ^ 2) / (2 * sigma_barsq ^ 2)
  }
  if (statonly) return(teststat)
  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = method), class = "htest")
  broom::tidy(rval)
}