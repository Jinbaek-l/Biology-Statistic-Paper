#' simControl: Main controller for setting the nine conditions / generating simulation data (simGen) / Applying DVG analysis (simDvarSeq)

## Set repositories
setRepositories(ind=1:8)

## Set directory 
dir="my directory"
setwd(dir)
getwd()

## Load library
source("simDvarSeq.R")
source("simGen.R")
source("DvarSeq.R")

library(doParallel)
library(foreach) # parallel for simulation DvarSeq
library(dplyr)
library(parallel) # parallel for common DvarSeq
library(edgeR)

## Generate baseline proportions for desired number of genes
load(url("http://bioinf.wehi.edu.au/voom/qAbundanceDist.RData"))

## Default setting
args <- list()
args[["case"]] <- "Set case"
args[["thread"]] <- 50
args[["nsim"]] <- 30
args[["nsm"]] <- c(10,20,30,50,100,200,300,500,1000,3000)
args[["BCVw"]] <- c(0.1,0.25,0.5,1.0,1.5,2.0)
args[["ngroup"]] <- 2
args[["ecInd"]] <- c(1,10000)
args[["unblratio"]] <- NULL
args[["uneqLSw"]] <- NULL
args[["TCVd"]] <- NULL
args[["RNsd"]] <- NULL
args[["ZIratio"]] <- NULL
args[["Norm"]] <- NULL


#' @description [Pre-screening filter 1]
#' @details BP server + 3 time simulations + sample size 500 + BCV 1.0 + 32 threads parallel
methods <- c("Harvey","BreuschPagan","CookWeisberg","Ftest","Hartley","Bartlett","CochranG","RamseyBAMSET","Levene","Glejser","DiffVar","OBrienHines","BrownForsythe","TrimmedLevene","GoldfeldQuandt","Koenker","FlignerKilleen","Horn","Yuce","White","Verbyla","Anscombe","Zhou","LiYao","Bickel","RackauskasZuokas","WilcoxKeselman","Szroeter","EvansKing","DiblasiBowman","CarapetoHolt","HarrisonMcCabe","MDSeq","DiffDist","SimonoffTsai")

length(methods)

args[["case"]] <- "PreScreening_1"
args[["nsim"]] <- 3
args[["nsm"]] <- c(500)
args[["BCVw"]] <- c(1.0)
args[["thread"]] <- 32

simGen(args, dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}


#' @description [Pre-screening filter 2]
#' @details IU server + 3 time simulations + sample size 3000 + BCV 1.0 + 32 threads parallel
methods <- c("Harvey","BreuschPagan","CookWeisberg","Ftest","Hartley","Bartlett","CochranG","RamseyBAMSET","Levene","Glejser","DiffVar","OBrienHines","BrownForsythe","TrimmedLevene","GoldfeldQuandt","Koenker","FlignerKilleen","Horn","Yuce","White","Verbyla","Anscombe","Zhou","LiYao","MDSeq","WilcoxKeselman")

length(methods)

args[["case"]] <- "PreScreening_2"
args[["nsim"]] <- 3
args[["nsm"]] <- c(3000)
args[["BCVw"]] <- c(1.0)
args[["thread"]] <- 32

simGen(args, dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}


#' @description [MDSeq test]
#' @details IU server + 1 time simulations + various sample size with BCV 1.0 + 32 threads parallel

args[["case"]] <- "MDSeq_test"
args[["nsim"]] <- 1
args[["nsm"]] <- c(10,20,30,50,100,200,300,500,550,600)
args[["BCVw"]] <- c(1.0)
args[["thread"]] <- 32

simGen(args, dir = "../")
simDvarSeq(args,method = "MDSeq", dir = "../")



#' @description [Horn test]
#' @details RV server + 1 time simulations + various sample size with BCV 1.0 + 32 threads parallel

args[["case"]] <- "Horn_test"
args[["nsim"]] <- 1
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
args[["BCVw"]] <- c(1.0)
args[["thread"]] <- 32

simGen(args, dir = "../")
simDvarSeq(args,method = "Horn", dir = "../")


#######################################################################################

methods <- c("Harvey","BreuschPagan","CookWeisberg","Ftest","Hartley","Bartlett","CochranG","RamseyBAMSET","Levene","Glejser","DiffVar","OBrienHines","BrownForsythe","TrimmedLevene","GoldfeldQuandt","Koenker","FlignerKilleen","Yuce","Anscombe","Verbyla","White","Zhou")
length(methods)

#' @description [Baseline]
args[["case"]] <- "Baseline"
simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")


#' @description [TCVd_X] X = {1,2,3,5,10}
args[["case"]] <- "TCVd_10"
args[["TCVd"]] <- 10
simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")


#' @description [RN_sdX] X = {0.1,0.25,0.5,0.75,1}
args[["case"]] <- "RN_sd1"
args[["RNsd"]] <- 1
simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = "Zhou", dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")


#' @description [LSw_X] X = {0.1,0.25,0.5,0.75}
args[["case"]] <- "LSw_0.75"
args[["uneqLSw"]] <- 0.75
simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = "Zhou", dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")


#' @description [Normalization_X] X = {TMM,RLE,UQ,VSD,CLR}
args[["case"]] <- "Normalization_CLR"
args[["Norm"]] <- "CLR"
simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")


#' @description [SSratioX] X = {64,73,82,91}
args[["case"]] <- "SSratio91"
args[["unblratio"]] <- c(9,1)
args[["nsm"]] <- c(30,50,100,200,300,500,1000,3000)

simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}
simDvarSeq(args,method = "Horn", dir = "../")
args[["nsm"]] <- c(30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")


#' @description [RE_logCP_X] X = {underNEG1,NEG1to0,0to1,1to2,2to3,3to4}
ngenes <- 10000
baselineprop <- qAbundanceDist( (1:ngenes)/(ngenes+1) )
cpm <- log(cpm(baselineprop))
param <- list()
param[[1]] <- c(min(which(cpm < -1)),max(which(cpm < -1)))
for(i in c(2:6)){
  param[[i]] <- c( min(which( (i-3)<=cpm & cpm<(i-2))), max(which( (i-3)<=cpm & cpm<(i-2))))
}
title <- c("RE_logCPM_underNEG1","RE_logCPM_NEG1to0","RE_logCPM_0to1","RE_logCPM_1to2","RE_logCPM_2to3")
names(param) <- title

args[["case"]] <- "RE_logCPM_2to3"
args[["ecInd"]] <- param$RE_logCPM_2to3
simDvarSeq(args,method = "Zhou", dir = "../")


simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")


#' @description [ZI_X] X = {0.1,0.3,0.5,0.7,0.9}
args[["case"]] <- "ZI_0.9"
args[["ZIratio"]] <- 0.9
simGen(args,dir = "../")
for(meth in methods){
  simDvarSeq(args,method = meth, dir = "../")
}
args[["nsm"]] <- c(10,20,30,50,100,200,300,500)
simDvarSeq(args,method = "MDSeq", dir = "../")
args[["nsm"]] <- c(20,30,50,100,200,300,500,1000,3000)
simDvarSeq(args,method = "Horn", dir = "../")

