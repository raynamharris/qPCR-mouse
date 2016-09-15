##=============================================================================================================#
## Script created by Rayna Harris, rayna.harris@utexas.edu
## Script created in version R 3.3.1 
## This script is for analyzing data related to the mouse qpcr projects
##=============================================================================================================#

#install.packages("dplyr")
#install.packages("plyr")
#install.packages("MCMC.qpcr")
#install.packages("reshape2")
#install.packages("car")
library(dplyr)
library(plyr)
library(MCMC.qpcr)
library(reshape2)
library(car)


## set path to data dir
setwd("/Volumes/nsb/NSB_2016/4_MouseMolecular/Final Project/standardcurvedata")
#setwd("Z:/NSB_2016/4_MouseMolecular/Final\ Project/standardcurvedata")

## 1. Read standard curve files

OtxrChH <- read.csv("2016-07-27_StandardCurve-OtxR-CrH.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
V1arNr3c1 <- read.csv("2016-07-27_StandardCurve-V1ar-Nr3c1.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
V1brmtor <- read.csv("2016-07-28_StandardCurve-V1b-mtor.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
prkci <- read.csv("2016-07-21_213928_Rachel_Rpl19_Fmr1_Prkci.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
BdnfCamk2a <- read.csv("2016-07-17_184859_ArctoFmr1.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
Rpl19CREB1 <- read.csv("2016-07-17_184859_FostoRpl19.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
Hdac2 <- read.csv("2016-07-23_154302_Hdac2_Bdnf_eff_NLB.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]
prkcz <- read.csv("2016-07-20_115248_ShaynNicole_CREB_PkcZ_Htt.csv", header=TRUE, stringsAsFactors = FALSE)[ ,1:24]


## 2. Use rbind to bind data into a fill "raw data" dataframe
rawdata <- rbind(OtxrChH, V1arNr3c1, V1brmtor, prkci, BdnfCamk2a, Rpl19CREB1, Hdac2, prkcz)
names(rawdata)
str(rawdata)

## 3. Clean data

cleandata <- rawdata
cleandata <- rename(cleandata, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 
names(cleandata)
str(cleandata)
cleandata$dna <- as.numeric(cleandata$dna, na.rm = TRUE)
cleandata$cq <- as.numeric(cleandata$cq, na.rm = TRUE)
cleandata$gene <- as.factor(cleandata$gene)
str(cleandata)

## 4. Create dilutions dataframe with quantiry, target name, and ct. 

dilutions <- cleandata %>%
  filter(Task == "STANDARD") %>%
  select(dna, cq, gene) 
str(dilutions)

dilutions$gene <- recode(dilutions$gene, "'Nr3c1' = 'GR'")

## 5. Calculate primer efficiences with MCMC.qpcr PrimEFF function
PrimEff(dilutions) # makes a plot with the primer efficiencies
amp.eff <- PrimEff(dilutions) #creates a table with the primer efficiencies

