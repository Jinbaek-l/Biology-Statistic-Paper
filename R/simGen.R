#' simGen: Generate simulation count data  
#' 
#' @description
#' This function generates simulated RNA-seq data controlling for nine conditions (arguments) using a data generation method built by the WEHI Institute. 
#'
#' @param args	Arguments set to design a simulation experiment
#' @param dir	Directory to store each simulation data

simGen <- function(args,dir){
  
  # Make directory
  if(!dir.exists(paste0(dir,args$case))){
    dir.create(paste0(dir,args$case))
    dir.create(paste0(dir,args$case,"/simCount"))
    dir.create(paste0(dir,args$case,"/simDDGind"))
    dir.create(paste0(dir,args$case,"/simGroup"))
  }
  
  # Set basic parameter for parallel
  base <- list()
  for(x in args$BCVw){
    for(y in args$nsm){
      for(z in seq_len(args$nsim))
      base[[paste0("BCV",x,"N",y,"NSIM",z)]] <- c(x,y,z)
    }
  }
  
  # Setting for parallel
  cl <- makeCluster(args$thread)
  registerDoParallel(cl)

  packages <- c("data.table","dplyr","tidyr","stringr","limma","edgeR","DESeq2","doParallel","foreach")
  functions <- c("qAbundanceDist")
  
  print("Simulated data generation start")
  
  foreach(p=1:length(base), .packages = packages, .export = functions) %dopar% {
    
    # Set basic paramter
    case <- args$case
    BCVw <- base[[p]][1]
    dnsample <- base[[p]][2]
    nsim <- base[[p]][3]
    
    # Set design parameter
    ngroup <- args$ngroup
    ecInd <- args$ecInd
    unblratio <- args$unblratio
    uneqLSw <- args$uneqLSw
    TCVd <- args$TCVd
    RNsd <- args$RNsd
    ZIratio <- args$ZIratio
    Norm <- args$Norm
    
    # BEGIN SIM
    ngenes <- 10000
    nDDG <- ngenes*0.02 
    
    # Random DDG index
    ind <- sample(ecInd[1]:ecInd[2],nDDG)
    oind <- c(1:ngenes)[-ind]
    
    # Baseline proportion
    baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
    baselineprop <- baselineprop/sum(baselineprop)
    
    # Sample size
    if(!is.null(unblratio)){
      nsample <- c(round(dnsample*(unblratio[1]/10),0),round(dnsample*(unblratio[2]/10),0))
    }else{
      nsample <- rep(round(dnsample/2,0),ngroup)
    }
    
    # Group design
    group <- factor(rep(seq_len(ngroup),nsample))
    
    # Library size
    if(!is.null(uneqLSw)){
      expected.lib.size <- rep(c(22e6,22e6*uneqLSw),sum(nsample)/ngroup)
    }else{
      expected.lib.size <- rep(11e6,sum(nsample))
    }
    
    # Calculate Expected counts 
    baselinepropGroup <- list()
    for(i in seq_len(ngroup)){
      baselinepropGroup[[i]] <- baselineprop # Create baseline proportion by group
    }
    groupExpcount <- list()
    groupExpcount[[1]] <- matrix(baselinepropGroup[[1]],ngenes,1) %*% matrix(expected.lib.size[which(group==1)],1,nsample[1])
    groupExpcount[[2]] <- matrix(baselinepropGroup[[2]],ngenes,1) %*% matrix(expected.lib.size[which(group==2)],1,nsample[2])
    expectedCount <- Reduce(cbind,groupExpcount)
    
    # NB dispersion
    BCV0 <- 0.2+1/sqrt(expectedCount)
    BCV <- matrix(0,ngenes,sum(nsample))
    
    ## Biological variation 
    BCV[oind,] <- BCV0[oind,]*exp( rnorm(length(oind),mean=0,sd=0.25)/2 )
    
    # Control Biological variation in each group to generate DDG 
    BCV[ind,which(group==1)] <- BCV0[ind,which(group==1)]*exp( rnorm(length(ind),mean=0,sd=BCVw)/2 )
    BCV[ind,which(group==2)] <- BCV0[ind,which(group==2)]*exp( rnorm(length(ind),mean=0,sd=BCVw)/2 )
    shape <- 1/BCV^2
    scale <- expectedCount/shape
    mu <- matrix(rgamma(ngenes*sum(nsample),shape=shape,scale=scale),ngenes,sum(nsample))
    
    # Technical variation
    if(!is.null(TCVd)){
      # [Mu] is the number of failures until [size] times of successes, probability can be calculated with two parameters
      counts <- matrix(rnbinom(ngenes * sum(nsample), size = mu / (TCVd - 1), mu = mu), ngenes, sum(nsample))
      # NA will produced if size or mu is to small to be valid, Change NA in to 0 
      counts[is.na(counts)] = 0
    }else{
      counts <- matrix(rpois(ngenes * sum(nsample), lambda = mu), ngenes, sum(nsample))
    }
    
    # Zero-inflated pattern applied or not
    if(!is.null(ZIratio)){
      for(gene in seq_len(ngenes)){
        ZeroSample <- which(counts[gene,]==0)
        TrueZeroLen <- length(ZeroSample) # True number of zero samples already exist in data
        DemZeroLen <- sum(nsample)*ZIratio # Demanded number of zero samples
        diff <- TrueZeroLen-DemZeroLen
        if(diff<0){ 
          # No enough zero samples in simulation data: Take "diff" out of the nonzero samples and change it to zero
          allInd <- c(1:sum(nsample))
          if(length(ZeroSample)!=0){
            nonZero <- allInd[-ZeroSample]
          }else{
            nonZero <- allInd
          }
          IndToChange <- sample(nonZero,abs(diff))
          counts[gene,IndToChange] <- 0
        }else if(diff>0){
          # To many zero samples in simulation data: Taking "diff" out of zero samples and converting them into non-zero values
          IndToChange <- sample(ZeroSample,diff)
          counts[gene,IndToChange] <- counts[gene,IndToChange] + 1
        }
      }
    }
    
    # Random noise or not
    if(!is.null(RNsd)){
      counts <- apply(counts,1,function(x) x <- x * exp(rnorm(sum(nsample), mean = 0, sd = RNsd)))
      counts <- as.data.frame(t(counts))
    }else{
      counts <- as.data.frame(counts)
    }
    
    # Filtering genes: Variance zero in both groups
    rmInd <- which(apply(counts,1,function(x) (var(x[which(group==1)]) == 0) && var(x[which(group==2)]) == 0)==TRUE)
    for(k in rmInd){
      counts[k,which.max(mu[k,which(group==1)])] <- counts[k,which.max(mu[k,which(group==1)])] + 1
      counts[k,nsample[1]+which.max(mu[k,which(group==2)])] <- counts[k,nsample[1]+which.max(mu[k,which(group==2)])] + 1
    }
    
    # Apply normalization method or not
    if(!is.null(Norm)){
      if(Norm=="TMM"){
        ## TMM(Trimmed mean of m-values) with "edgeR"
        obj <- edgeR::DGEList(counts=counts,group=group)
        obj <- edgeR::calcNormFactors(obj,method = "TMM")
        counts <- as.data.frame(edgeR::cpm(obj,log = TRUE))
        
      }else if(Norm=="RLE"){
        ## RLE(Relative log expression) with "edgeR"
        obj <- edgeR::DGEList(counts=counts,group=group)
        obj <- edgeR::calcNormFactors(obj,method = "RLE")
        counts <- as.data.frame(edgeR::cpm(obj,log = TRUE))
        
      }else if(Norm=="UQ"){
        ## UQ(Upper quartile) with "edgeR"
        obj <- edgeR::DGEList(counts=counts,group=group)
        obj <- edgeR::calcNormFactors(obj,method = "upperquartile")
        counts <- as.data.frame(edgeR::cpm(obj,log = TRUE))
        
      }else if(Norm=="CSS"){
        ## CSS(Cumulative sum scaling) with "metagenomeSeq"
        obj <- metagenomeSeq::newMRexperiment(counts)
        obj <- metagenomeSeq::cumNorm(obj)
        counts <- as.data.frame(MRcounts(obj, norm=TRUE,log=TRUE))
        
      }else if(Norm=="VSD"){
        ## VSD (Variance stabilized DESeq) normalization
        colData <- data.frame(group = group)
        dds <- DESeq2::DESeqDataSetFromMatrix(counts, colData = colData, design = ~ group)
        dds <- BiocGenerics::estimateSizeFactors(dds)
        vsd <- DESeq2::vst(dds)
        counts <- assay(vsd)
        
      }else if(Norm=="CLR"){
        ## CLR (Centered log ratio) transformation
        counts <- t(compositions::clr(t(counts)))
        
      }else{
        print("Wrong normalization method")
      }
    }
    # END SIM
    DDGind <- ind 
    rawcount <- counts
    
    save(rawcount, file = paste0(dir,case,"/simCount/",case,"_BCV",BCVw,"_SS",dnsample,"_sim",nsim,"_rawcount.RData"))
    save(DDGind, file = paste0(dir,case,"/simDDGind/",case,"_BCV",BCVw,"_SS",dnsample,"_sim",nsim,"_DDGind.RData"))
    save(group, file = paste0(dir,case,"/simGroup/",case,"_BCV",BCVw,"_SS",dnsample,"_sim",nsim,"_group.RData"))
  }
  stopCluster(cl)
  
  print("Simulating data done!")
}