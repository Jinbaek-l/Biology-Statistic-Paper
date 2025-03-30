#' @description Copied from skedastic::bamset (GPL-3)
#' @source https://github.com/cran/skedastic/blob/master/R/bamset.R

bamset <- function(mainlm, k = 3, deflator = NA, correct = TRUE,
                   omitatmargins = TRUE, omit = NA,
                   categorical = FALSE, statonly = FALSE) {
  
  processmainlm(m = mainlm, needy = FALSE)
  
  hasintercept <- columnof1s(X)
  if (inherits(mainlm, "list")) {
    if (hasintercept[[1]]) {
      if (hasintercept[[2]] != 1) stop("Column of 1's must be first column
                                         of design matrix")
      colnames(X) <- c("(Intercept)", paste0("X", 1:(p - 1)))
    } else {
      colnames(X) <- paste0("X", 1:p)
    }
  }
  
  checkdeflator(deflator, X, p, hasintercept[[1]])
  
  if (!is.na(deflator) && !is.null(deflator)) {
    if (!is.na(suppressWarnings(as.integer(deflator)))) {
      deflator <- as.integer(deflator)
    }
    e <- e[order(X[, deflator])]
    X <- X[order(X[, deflator]), , drop = FALSE]
  }
  
  n <- nrow(X)
  nprime <- n - p
  
  if (categorical) {
    if (is.factor(X[, deflator])) {
      lev <- levels(X[, deflator])
    } else {
      lev <- unique(sort(X[, deflator]))
    }
    k <- length(lev)
    min_subset_size <- as.integer(nprime / k) # called r_1 in Ramsey (1969)
    nprime_modk <- nprime %% k
    v <- table(X[, deflator])
    
  } else {
    min_subset_size <- as.integer(nprime / k) # called r_1 in Ramsey (1969)
    nprime_modk <- nprime %% k
    
    # slightly different definition of v_i compared with Ramsey (1969)
    # this one makes the subsets more equitable in size (not differing by more than 1)
    v <- rep(min_subset_size, k) + c(rep(1, nprime_modk), rep(0, k - nprime_modk))
  }
  sub_ind <- vector("list", k)
  sub_ind[[1]] <- seq.int(from = 1, by = 1, length.out = v[1])
  for (j in 2:k) {
    sub_ind[[j]] <- seq.int(from = max(sub_ind[[j - 1]]) + 1, by = 1,
                            length.out = v[j])
  }
  
  if (omitatmargins)
    omit <- margin_indices(v, p, sub_ind,
                           categ = categorical)
  
  res <- blus(mainlm = list("X" = X, "e" = e), omit)
  
  s_sq <- unlist(lapply(sub_ind, function(i, e)
    sum(e[i] ^ 2, na.rm = TRUE), e = res)) / v
  s_sq_tot <- sum(res ^ 2, na.rm = TRUE) / nprime
  
  teststat <- nprime * log(s_sq_tot) - sum(v * log(s_sq))
  if (correct) {
    scaling_constant <- 1 + (sum(1 / v) - 1 / nprime) / (3 * (k - 1))
    teststat <- teststat / scaling_constant
  }
  
  if (statonly) return(teststat)
  
  df <- k - 1
  pval <- stats::pchisq(teststat, df = df, lower.tail = FALSE)
  
  rval <- structure(list(statistic = teststat, parameter = df, p.value = pval,
                         null.value = "Homoskedasticity",
                         alternative = "greater", method = "BAMSET"),
                    class = "htest")
  broom::tidy(rval)
}