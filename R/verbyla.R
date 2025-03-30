#' @description Copied from skedastic::verbyla (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/verbyla.R

verbyla <- function(mainlm, auxdesign = NA, statonly = FALSE) {
  
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
  if (hasintercept[[1]]) {
    Z <- Z[, -hasintercept[[2]], drop = FALSE]
  }
  
  p <- ncol(X)
  q <- ncol(Z)
  n <- nrow(Z)
  
  Z <- cbind(1, Z)
  
  M <- fastM(X, n)
  sigma_hatbar <- sum(e ^ 2) / (n - p)
  term1 <- t(t(e ^ 2 / sigma_hatbar - diag(M)))
  
  mat_to_invert <- t(Z) %*% (M ^ 2) %*% Z
  if (det(mat_to_invert) == 0) {
    Z <- Z[, -1, drop = FALSE]
    mat_to_invert <- t(Z) %*% (M ^ 2) %*% Z
    if (det(mat_to_invert) == 0) stop("Linear dependency detected in auxiliary matrix")
    message("Intercept not included in auxiliary design in order to avoid linear dependency")
  }
  
  teststat <- as.double(1 / 2 * t(term1) %*% Z %*%
                          solve(mat_to_invert)
                        %*% t(Z) %*% term1)
  if (statonly) return(teststat)
  
  pval <- stats::pchisq(teststat, df = q, lower.tail = FALSE)
  rval <- structure(list(statistic = teststat, parameter = q, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = "Verbyla"),
                    class = "htest")
  broom::tidy(rval)
}