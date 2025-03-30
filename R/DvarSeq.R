#' DvarSeq: Differential variability analysis for sequencing data 
#' 
#' @description
#' This function integrates 24 different methods that have been developed for variability difference tests from 1935 to 2024
#' Some tests are customized, and most cite the original code as it is
#' Each test supports parallel processing for multiple genes and minimizes package dependence
#' 
#' @param y	A matrix consisting of genes in the row and samples in the column, which is bulk RNA-seq data. Pseudo single cell RNA-seq data is also allowed
#' @param group	Corresponding group information for each sample; only two groups are allowed
#' @param method	Character. Specifies the variability difference test to use. See the documentation for the list of supported method names.
#' @param parallel	Logical. Whether to enable parallel processing. If TRUE, computations across genes will be distributed across multiple threads. Requires the `threads` parameter to be set.
#' @param threads	Integer. The number of threads to use when parallel is TRUE. Ignored when parallel is FALSE.
#' 
#' @return Returns the test statistic, p-value and false discovery rate-adjusted p-value for whole genes

DvarSeq <- function(y,group,method,parallel=FALSE,threads=1){

  ## make geneName vector
  geneName <- rownames(y)
  
  ## check group type is factor
  if (!is.factor(group)) {
    warning(deparse(substitute(group)), " coerced to factor.")
    group <- as.factor(group)
  }
  
  ## check NA or Inf in data frame
  ind <- which(as.vector(apply(y, 1, function(x) any(any(is.na(x) | is.infinite(x))))))
  if(!identical(ind,integer(0))){
    stop(paste0("NA or Inf must be controlled\n",toString(geneName[ind])))
  }
  
  ## check groups with no more than two samples
  tmp <- names(which(table(group) <= 2))
  if(!identical(tmp,character(0))){
    stop(paste0("Variance cannot be obtained if the number of samples is 2 or less\n",toString(tmp)))
  }
  
  ## check gene with zero variance in whole group (ex.000/000,111/333)
  geneError <- geneName[which(apply(y,1,function(x) var(x)==0))]
  
  ## provide gene list to be removed
  if(length(geneError)!=0){
    stop("All groups of following gene have zero variance. It must be removed\n",paste0(geneError,", "),call. = FALSE)
  }
  
  ## Defines the corresponding function for each method
  tool <- FALSE
  switch(method,
         "Ftest" = { # F-test
           test_func <- function(x) { var.test(unlist(x[which(group==1)], use.names = FALSE), unlist(x[which(group==2)], use.names = FALSE)) }
         },
         "Hartley" = { # Hartley's F max test customized
           test_func <- function(x) { hartleyCUSTOM(x, group) }
         },
         "Bartlett" = { # Bartlett's test
           test_func <- function(x) { bartlett.test(x ~ group) }
         },
         "CochranG" = { # Cochran's G test
           test_func <- function(x) { cochranTest(x ~ group) }
         },
         "Levene" = { # Levene's test
           test_func <- function(x) { leveneTest(x ~ group, center = "mean") }
         },
         "TrimmedLevene" = { # Trimmed(0.25) Levene's test
           test_func <- function(x) { leveneTest(x ~ group, center = "mean", trim = 0.25) }
         },
         "BrownForsythe" = { # Brown-Forsythe test
           test_func <- function(x) { leveneTest(x ~ group, center = "median") }
         },
         "OBrienHines" = { # O'Brien, Hines test
           test_func <- function(x) { levene.test(x, group, correction.method = "zero.correction") }
         },
         "FlignerKilleen" = { # Fligner-Killeen's test
           test_func <- function(x) { fligner.test(x ~ group) }
         },
         "Koenker" = { # Studentized Breusch-Pagan test
           test_func <- function(x) { breusch_pagan(lm(x ~ group), koenker = TRUE) }
         },
         "BreuschPagan" = { # Breusch-Pagan test
           test_func <- function(x) { breusch_pagan(lm(x ~ group), koenker = FALSE) }
         },
         "GoldfeldQuandt" = { # Goldfeld-Quandt F test
           test_func <- function(x) { goldfeld_quandt(lm(x ~ group), alternative = "two.sided") }
         },
         "Glejser" = { # Glejser test
           test_func <- function(x) { glejser(lm(x ~ group)) }
         },
         "White" = { # White's test
           test_func <- function(x) { white(lm(x ~ group)) }
         },
         "Anscombe" = { # Anscombe's test
           test_func <-  function(x) { anscombe(lm(x ~ group), studentise = FALSE) }
         },
         "RamseyBAMSET" = { # Ramsey's BAMSET test
           test_func <- function(x) { bamset(lm(x ~ group)) }
         },
         "CookWeisberg" = { # Cook-Weisberg score test
           test_func <- function(x) { cook_weisberg(lm(x ~ group)) }
         },
         "Harvey" = { # Harvey test Customized
           test_func <- function(x) { harveyCUSTOM(lm(x ~ group)) }
         },
         "Horn" = { # Horn's test
           test_func <- function(x) { horn(lm(x ~ group), alternative = "two.sided") }
         },
         "Verbyla" = { # Verbyla's test
           test_func <- function(x) { verbyla(lm(x ~ group)) }
         },
         "Yuce" = { # Yuce's test
           test_func <- function(x) { yuce(lm(x ~ group, method = "B")) }
         },
         "Zhou" = { # Zhou, Song, and Thompson’s Test Customized
           test_func <- function(x) { zhouCUSTOM(lm(x ~ group)) }
         },
         
         "DiffVar" = { # DiffVar test from limma
           tool <- TRUE
           design <- model.matrix(~group)
           fit <- varFit(y, design,coef=c(1,2))
           result <- topVar(fit, coef=2, n=nrow(y), sort = FALSE)
           table <- data.frame(STAT = result$t, P = result$P.Value)
         },
         
         "MDSeq" = { # MDSeq
           tool <- TRUE
           contrasts <- get.model.matrix(group)
           if(parallel){
             fit <- MDSeq(y, contrast=contrasts, mc.cores = threads)
           }else {
             fit <- MDSeq(y, contrast=contrasts, mc.cores = 1)
           }
           result <- extract.ZIMD(fit, compare=list(A="1",B="2"))
           table <- data.frame(STAT = result$Statistics.dispersion, P = result$Pvalue.dispersion)
         },
         
         warning("No corresponding method exist!!!")
  )
  
  if(!tool){
    if(parallel){
      cl <- makeCluster(threads)
      suppressWarnings(
        clusterEvalQ(cl, {
          source("DvarSeq.R")
          source("fileLoad.R")
        })
      )
      
      suppressWarnings(
        clusterExport(cl, c("group"))
      )
      suppressWarnings(
        result <- parApply(cl, y, 1, function(x) test_func(x))
      )
      stopCluster(cl)
    }else{
      result <- apply(y, 1, function(x) test_func(x))
    }
    
    if(method %in% c("Levene","TrimmedLevene","BrownForsythe")){
      table <- data.frame(STAT = sapply(result, function(x) x$`F value`[1]), P = sapply(result, function(x) x$`Pr(>F)`[1]))
    }else{
      table <- data.frame(STAT = sapply(result, function(x) x$statistic), P = sapply(result, function(x) x$p.value))
    }
  }
  rownames(table) <- geneName
  
  # Calculate FDR adjusted p-value
  adjPvlTable <- p.adjust(table$P, method = "BH") # FDR correction
  
  table <- table %>% dplyr::mutate(AdjP=adjPvlTable)

  return(table)
}