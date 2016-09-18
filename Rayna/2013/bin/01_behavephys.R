# Part 1: Reading and analyzing beahvior and physiology data

## wrangle the data ----
setwd("~/Github/qPCR-mouse/Rayna/2013/data")
wt <- read.csv("01_WTbehavephys.csv", header=TRUE, stringsAsFactors = FALSE)
fmr1 <- read.csv("01_FMR1behavephys.csv", header=TRUE, stringsAsFactors = FALSE)

## fixing columns both fmr1 and wt are as equivilent as possible
wt$filename <- as.factor(paste(wt$filename, wt$X, sep="_")) ## mgerge 2 columns 
wt$X <- NULL ## delete now obsolete column
wt$X.1 <- NULL
fmr1$X <- NULL
wt$next. <- NULL
fmr1$next. <- NULL
fmr1$PopSpike.max........mV. <- NA  ## add column
fmr1$PospSpike.Probability <- NA  ## add column
fmr1 <- fmr1[-c(75:80), ] ## removing rows
names(wt)
names(fmr1)

## give all non-summary columns with all a better name!!! ----
names(wt)[3] <- "TotalTime"
names(wt)[4] <- "TotalPath"
names(wt)[5] <- "Entrances"
names(wt)[6] <- "TimeToFirstEntrance"
names(wt)[7] <- "PathToFirstEntrance"
names(wt)[8] <- "SpeedToFirstEntrance"
names(wt)[9] <- "EntrancePerDistance"
names(wt)[10] <- "TotalShocks"
names(wt)[11] <- "TimeToFirstShock"
names(wt)[12] <- "PathToFirstShock"
names(wt)[13] <- "Speed"
names(wt)[14] <- "SDofSpeed"
names(wt)[15] <- "LineArity"
names(wt)[16] <- "MaxAvoidTime"
names(wt)[17] <- "MaxAvoidPath"
names(wt)[18] <- "TimeToSecondEntrance"
names(wt)[19] <- "PathToSecondEntrance"
names(wt)[20] <- "SpeedToSecondEntranc"
names(wt)[21] <- "TimeTARG"
names(wt)[22] <- "pTimeTARG"
names(wt)[23] <- "pTimeCCW"
names(wt)[24] <- "pTimeOPP"
names(wt)[25] <- "pTimeCW"
names(wt)[26] <- "RayleigLength"
names(wt)[27] <- "RayleigAngle"
names(wt)[28] <- "PolarAveVal"
names(wt)[29] <- "PolarSdVal"
names(wt)[30] <- "PolarMinVal"
names(wt)[31] <- "PolarMinBin"
names(wt)[32] <- "MinLoBin"
names(wt)[33] <- "MinHiBin"
names(wt)[34] <- "PolarMaxVal"
names(wt)[35] <- "PolarMaxBin"
names(wt)[36] <- "MaxLoBin"
names(wt)[37] <- "MaxHiBin"
names(wt)[38] <- "AnnularMinVal"
names(wt)[39] <- "AnnularMinBin"
names(wt)[40] <- "AnnularMaxVal"
names(wt)[41] <- "AnnularMaxBin"
names(wt)[42] <- "AnnularAvg"
names(wt)[43] <- "AnnularSD"
names(wt)[44] <- "AnnularSkewnes"
names(wt)[45] <- "AnnularKurtosis"
names(wt)[46] <- "ind"
names(wt)[47] <- "genotype"
names(wt)[48] <- "APA"
names(wt)[49] <- "TotalLocomotorActivity"
names(wt)[50] <- "LastLocomotorActivity"
names(wt)[51] <- "TotalPunishment"
names(wt)[53] <- "LastPunishment"
names(wt)[54] <- "LastTotalEntrances"
names(wt)[55] <- "T1Retention"
names(wt)[56] <- "T2Retention"
names(wt)[57] <- "IO_Max"
names(wt)[58] <- "LTP_Baseline"
names(wt)[59] <- "LTP_Baseline_SD"

## repeat for fmr1 columsn
names(fmr1)[3] <- "TotalTime"
names(fmr1)[4] <- "TotalPath"
names(fmr1)[5] <- "Entrances"
names(fmr1)[6] <- "TimeToFirstEntrance"
names(fmr1)[7] <- "PathToFirstEntrance"
names(fmr1)[8] <- "SpeedToFirstEntrance"
names(fmr1)[9] <- "EntrancePerDistance"
names(fmr1)[10] <- "TotalShocks"
names(fmr1)[11] <- "TimeToFirstShock"
names(fmr1)[12] <- "PathToFirstShock"
names(fmr1)[13] <- "Speed"
names(fmr1)[14] <- "SDofSpeed"
names(fmr1)[15] <- "LineArity"
names(fmr1)[16] <- "MaxAvoidTime"
names(fmr1)[17] <- "MaxAvoidPath"
names(fmr1)[18] <- "TimeToSecondEntrance"
names(fmr1)[19] <- "PathToSecondEntrance"
names(fmr1)[20] <- "SpeedToSecondEntranc"
names(fmr1)[21] <- "TimeTARG"
names(fmr1)[22] <- "pTimeTARG"
names(fmr1)[23] <- "pTimeCCW"
names(fmr1)[24] <- "pTimeOPP"
names(fmr1)[25] <- "pTimeCW"
names(fmr1)[26] <- "RayleigLength"
names(fmr1)[27] <- "RayleigAngle"
names(fmr1)[28] <- "PolarAveVal"
names(fmr1)[29] <- "PolarSdVal"
names(fmr1)[30] <- "PolarMinVal"
names(fmr1)[31] <- "PolarMinBin"
names(fmr1)[32] <- "MinLoBin"
names(fmr1)[33] <- "MinHiBin"
names(fmr1)[34] <- "PolarMaxVal"
names(fmr1)[35] <- "PolarMaxBin"
names(fmr1)[36] <- "MaxLoBin"
names(fmr1)[37] <- "MaxHiBin"
names(fmr1)[38] <- "AnnularMinVal"
names(fmr1)[39] <- "AnnularMinBin"
names(fmr1)[40] <- "AnnularMaxVal"
names(fmr1)[41] <- "AnnularMaxBin"
names(fmr1)[42] <- "AnnularAvg"
names(fmr1)[43] <- "AnnularSD"
names(fmr1)[44] <- "AnnularSkewnes"
names(fmr1)[45] <- "AnnularKurtosis"
names(fmr1)[46] <- "ind"
names(fmr1)[47] <- "genotype"
names(fmr1)[48] <- "APA"
names(fmr1)[49] <- "TotalLocomotorActivity"
names(fmr1)[50] <- "LastLocomotorActivity"
names(fmr1)[51] <- "TotalPunishment"
names(fmr1)[53] <- "LastPunishment"
names(fmr1)[54] <- "LastTotalEntrances"
names(fmr1)[55] <- "T1Retention"
names(fmr1)[56] <- "T2Retention"
names(fmr1)[57] <- "IO_Max"
names(fmr1)[58] <- "LTP_Baseline"
names(fmr1)[59] <- "LTP_Baseline_SD"

names(wt)
names(fmr1)

## remove non-data rows
wt <- filter(wt, grepl("oom", filename))
fmr1 <- filter(fmr1, grepl("oom", filename))

## merge the two ----
behavephy <- rbind(wt, fmr1)

## subset just summary data and just the every trial data
summary <- filter(behavephy, grepl("Retest|retention", filename))
summary <- summary[c(1,46:59)]

notsummarydata <- behavephy[c(1:48)]
  