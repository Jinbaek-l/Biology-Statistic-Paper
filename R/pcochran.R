#' @description Copied from PMCMRplus::pcochran (GPL-3)
#' @source https://github.com/cran/PMCMRplus/blob/master/R/Cochran.R


pcochran <-function (q, k, n, lower.tail = TRUE, log.p = FALSE)
{
  qF <- (1 / q - 1) / (k - 1)
  pval <- 1 - k * pf(q = qF,
                     df1 = (k - 1) * (n - 1),
                     df2 = n - 1,
                     lower.tail=TRUE)
  
  if (!lower.tail) {
    pval <- 1 - pval
  }
  if (log.p){
    pval <- log(pval)
  }
  return(pval)
}