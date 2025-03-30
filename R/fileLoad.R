#' Packages and functions that must be imported

#' @description Package for DvarSeq (required by each test are specified)
library(stats) # Ftest, Bartlett, FlignerKilleen
source("hartleyCUSTOM.R") # Hartley
library(SuppDists) # Hartley
source("cochranTest.R") # CochranG
source("pcochran.R") # CochranG
source("leveneTest.R") # Levene, TrimmedLevene, BrownForsythe
source("levene.test.R") # OBrienHines
source("breusch_pagan.R") # Koenker, BreuschPagan
source("myutils.R") # Koenker, BreuschPagan, GoldfeldQuandt, Glejser, White, Anscombe, RamseyBAMSET, CookWeisberg, Horn, Verbyla, Yuce, Zhou
source("goldfeld_quandt.R") # GoldfeldQuandt
source("twosidedpval.R") # GoldfeldQuandt
source("glejser.R") # Glejser
source("white.R") # White
source("anscombe.R") # Anscombe
source("bamset.R") # RamseyBAMSET
source("blus.R") # RamseyBAMSET
source("cook_weisberg.R") # CookWeisberg
source("harveyCUSTOM.R") # Harvey
source("horn.R") # Horn
source("trend.R") # Horn
source("verbyla.R") # Verbyla
source("yuce.R") # Yuce
source("zhouCUSTOM.R") # Zhou
source("DiffVar.R") # DiffVar
library(MDSeq) # MDSeq