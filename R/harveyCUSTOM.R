#' @description Adapted from skedastic::harvey (GPL-3)
#' @source https://rdrr.io/github/tjfarrar/skedastic/src/R/harvey.R
#' @include myutils.R
#' @include psi.R
#' @include broom package

#' @details In the auxiliary regression, the response variable was log-transformed.  Since taking the log of zero or a negative value results in -Inf, which causes errors, any -Inf values were replaced with the most negative double-precision number representable in R

harveyCUSTOM <- function(mainlm, auxdesign = NA, statonly = FALSE) {
  
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
  
  p <- ncol(Z) - 1
  n <- nrow(Z)
  auxresponse <- log(e ^ 2)
  
  auxresponse[auxresponse == -Inf] <- -.Machine$double.xmax  ## Code modified!
  
  auxres <- stats::lm.fit(Z, auxresponse)$residuals
  
  teststat <- (sum(auxresponse ^ 2) - n * mean(auxresponse) ^ 2
               - sum(auxres ^ 2)) / pracma::psi(1, 1 / 2)
  if (statonly) return(teststat)
  
  pval <- stats::pchisq(teststat, df = p, lower.tail = FALSE)
  
  rval <- structure(list(statistic = teststat, parameter = p, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater"), class = "htest")
  broom::tidy(rval)
}