#' @description Copied from skedastic::yuce (GPL-3)
#' @source https://github.com/tjfarrar/skedastic/blob/master/R/yuce.R


yuce <- function(mainlm, method = c("A", "B"), statonly = FALSE) {
  
  processmainlm(m = mainlm, needy = FALSE)
  
  method <- match.arg(toupper(method), c("A", "B"))
  n <- length(e)
  
  if (method == "A") {
    alt <- "greater"
    teststat <- sum(e ^ 2) / ((sum(abs(e)) ^ 2 * pi) /
                                (pi - 2 + 2 * n * (n - p)))
    if (statonly) return(teststat)
    fullmethod <- "Chi-Squared Test"
    pval <- stats::pchisq(teststat, df = n - p, lower.tail = FALSE)
  } else if (method == "B") {
    alt <- "two.sided"
    abspart <- sum(abs(e)) ^ 2 * pi / (pi - 2 + 2 * n * (n - p))
    teststat <- (sum(e ^ 2) / (n - p) - abspart) /
      sqrt(2 / (n - p) * abspart ^ 2)
    if (statonly) return(teststat)
    fullmethod <- "t Test"
    pval <- 2 * stats::pt(abs(teststat), df = n - p, lower.tail = FALSE)
  } else stop("Invalid `method` argument")
  
  rval <- structure(list(statistic = teststat, p.value = pval,
                         null.value = "Homoskedasticity", alternative = alt,
                         method = fullmethod), class = "htest")
  broom::tidy(rval)
}