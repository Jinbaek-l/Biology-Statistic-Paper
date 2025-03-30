#' simDvarSeq: Differential variability analysis for simulated data 
#' 
#' @description
#' This function loads all simulation datasets from the specified directory and applies the DvarSeq function to each dataset. Each simulation dataset is processed in parallel.
#' After performing differential variability analysis, the p-values are adjusted using empirical FDR correction. The performance of each method is integrated, summarized into one CSV file for each simulation dataset, and saved in a specified directory.
#'
#' @param args	Arguments for the simulation experiment design
#' @param method	Character. Specifies the variability difference test to use. See the documentation for the list of supported method names.
#' @param dir	Directory containing each simulation data set and saving analysis results

simDvarSeq <- function(args,method,dir)
{
  # Make directory
  if(!dir.exists(paste0(dir,args$case,"/simResult/",method))){
    dir.create(paste0(dir,args$case,"/simResult/",method),recursive = TRUE)
  }
  
  # Read simulated data
  counts.files <- list.files(paste0(dir,args$case,"/simCount/"), pattern = paste(paste0("_SS",args$nsm,"_"),collapse = "|"))
  groups.files <- list.files(paste0(dir,args$case,"/simGroup/"), pattern = paste(paste0("_SS",args$nsm,"_"),collapse = "|"))
  ind.files <- list.files(paste0(dir,args$case,"/simDDGind/"), pattern = paste(paste0("_SS",args$nsm,"_"),collapse = "|"))
  
  # Setting for simulation parallel
  cl <- makeCluster(args$thread)
  registerDoParallel(cl)
  
  packages <- c("doParallel","data.table","stringr","tidyr","utils","parallel","dplyr","stats","SuppDists","MDSeq")
  
  print("================================")
  print(paste0(method," test parallel start"))
  
  # Parallel for nBCV x nSample x nsim
  bind.result <- foreach(p=1:length(counts.files), .packages = packages, .combine = "rbind") %dopar% {
    source("DvarSeq.R")
    source("fileLoad.R")
    load(paste0(dir,args$case,"/simCount/",counts.files[p]))
    load(paste0(dir,args$case,"/simGroup/",groups.files[p]))
    load(paste0(dir,args$case,"/simDDGind/",ind.files[p]))
    
    bcv <- gsub(".*BCV(\\d+(?:\\.\\d+)?).*", "\\1", counts.files[p])
    ss <- gsub(".*SS(\\d+).*", "\\1", counts.files[p])
    sim <- gsub(".*sim(\\d+).*", "\\1", counts.files[p])
    
    # Handling NaN, Inf
    rawcount[is.na(rawcount) | rawcount=="Inf"] = NA
    
    if(method=="MDSeq"){
      # Shuffling the gene order
      rawcount <- rawcount[rev(1:nrow(rawcount)), ]
      DDGind <- match(DDGind, rev(1:10000))
    }
    
    # Perform DvarSeq analysis
    start_time <- Sys.time()
    result <- DvarSeq(y = rawcount, group = group, method = method, threads = args$thread)
    end_time <- Sys.time()

    write.csv(result,file = paste0(dir, args$case,"/simResult/",method,"/",args$case,"_BCV",bcv,"_SS",ss,"_sim",sim,"_pvlTable.csv"), row.names = FALSE)
    
    # Perform empirical FDR correction
    ngenes <- nrow(rawcount)
    ranked <- data.frame(matrix(0,nrow=ngenes,ncol=1))
    colnames(ranked) <- method
    fd <- ranked
    nd <- c(rep(0,1))
    names(nd) <- method
    nd.fd <- nd
    
    status <- rep(0,ngenes) # Denote DDG index
    status[DDGind] <- 1
    
    nd <- nd + sum(result$AdjP<0.05,na.rm = TRUE) # Find empirical prior with ranked p value
    od <- order(result$P)
    ranked <- status[od]
    
    fd <- fd + cumsum(1-ranked) # Calculate empirical FDR
    nd.fd <- fd[ceiling(nd+0.00001),1]
    
    tp <- nd-nd.fd # Calculate tp by total positive - false positive
    tp[which(tp==-1)] <- 0
    
    
    # fd: false discovery / nd.fd: false positive / nd: total positive
    df <- data.frame(BCV = as.numeric(bcv), n = as.numeric(ss), sim = as.numeric(sim),
                     design = args$case, method = method,
                     TP = tp, FP = nd.fd, FN = length(DDGind)-tp, TN = ngenes - length(DDGind) - nd.fd, ET = round(difftime(end_time, start_time, units = "mins"),5))
    rownames(df) <- NULL
    
    return(df)
  }

  stopCluster(cl)
  write.csv(bind.result,file = paste0(dir,args$case,"/simResult/",method,"/",method,"_simResult.csv"), row.names = FALSE)
}
