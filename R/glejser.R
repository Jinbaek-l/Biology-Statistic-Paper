#' @description Copied from skedastic::glejser (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/glejser.R

glejser <- function(mainlm, auxdesign = NA,
                    sigmaest = c("main", "auxiliary"), statonly = FALSE) {
  
  sigmaest <- match.arg(sigmaest, c("main", "auxiliary"))
  
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
  auxresponse <- abs(e)
  auxres <- stats::lm.fit(Z, auxresponse)$residuals
  if (sigmaest == "main") {
    sigma_hatsq <- sum(e ^ 2) / n
  } else if (sigmaest == "auxiliary") {
    sigma_hatsq <- sum(auxres ^ 2) / n
  }
  
  teststat <- (sum(auxresponse ^ 2) - n * mean(auxresponse) ^ 2 -
                 sum(auxres ^ 2)) / (sigma_hatsq * (1 - 2 / pi))
  if (statonly) return(teststat)
  
  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater"), class = "htest")
  broom::tidy(rval)
}