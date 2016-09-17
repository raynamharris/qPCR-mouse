## Want one script to wrangle and anlayze all the wt frm1 data for the qpcr project

# Parts
## Part 1 : behavior and physiology data
## Part 2 : qpcr data
## Part 3 : integraive analysis

# Part 2: Reading and analyzing qPCR data ----

## wrangle the gene expression qpcr data ----
setwd("~/Github/qPCR-mouse/Rayna/2013/data")
qpcr <- read.csv("02_qpcrdata.csv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
str(qpcr)
summary(qpcr)

## setting the factors
qpcr$sample <- as.factor(qpcr$sample)
qpcr$ind <- as.factor(qpcr$ind)
qpcr$time <- as.factor(qpcr$time)
qpcr$region <- as.factor(qpcr$region)
qpcr$APA <- as.factor(qpcr$APA)
qpcr$strained <- as.factor(qpcr$strained)
str(qpcr)
summary(qpcr)
head(qpcr)

## rename headers with rename() and factors with revalue() ----
# install.packages("dplyr")
library(dplyr)
qpcr <- rename(qpcr, genotype = strained) #note, this function masked when plyr installed. have to restart to rerun command

library(plyr)
qpcr$genotype <- revalue(qpcr$genotype, c("fmr1" = "FMR1-KO")) 
qpcr$genotype <- revalue(qpcr$genotype, c("wt" = "WT"))
head(qpcr)

## subset the data and droplevels! to have FMR1 and WT in separate files ----

## just FRM1 trained and untrained
FMR1KO <- filter(qpcr, genotype == "FMR1-KO")
FMR1KO <- droplevels(FMR1KO)
str(FMR1KO)

## just WT CA3
WT <- filter(qpcr, genotype == "WT", region == "CA3")
WT <- droplevels(WT)
str(WT)

## just WT and FRM1 CA1
CA1_3genes <- qpcr[c(1:10)]

## gene expression analysis with mcmc.qpcr ----
#install.packages("MCMC.qpcr")
library(MCMC.qpcr)

## read in dilution series, drop & renmae things and calculate gene effeciencies
dilutions <- read.csv("02_dilutions_CA1CA3.csv", header = TRUE)
head(dilutions)

dilutions <- filter(dilutions, gene != "dlg4", gene != "pkmz.conc", gene != "rpl19.conc", gene != "grim")
dilutions <- droplevels(dilutions)
dilutions$gene <- revalue(dilutions$gene, c("dlg4.conc" = "dlg4")) 
dilutions$gene <- revalue(dilutions$gene, c("grim.conc" = "grim"))
dilutions$gene <- revalue(dilutions$gene, c("fmr1.conc" = "fmr1"))
dilutions$gene <- revalue(dilutions$gene, c("rRNA18s" = "rRNA18S"))

PrimEff(dilutions) -> eff

## Analyze FMR1 data with cq2counts function and naive model ----
dd <- cq2counts(data=FMR1KO, genecols=c(8:19), condcols=c(1:7), effic=eff)
head(dd)

## use soft norm model to fit
soft <- mcmc.qpcr(
  data=dd,
  fixed="APA",
  controls=c("rRNA18S"),
  normalize = TRUE,
  pr=T,pl=T)

diagnostic.mcmc(model=soft, col="grey50", cex=0.8)
HPDsummary(soft, dd, relative=TRUE)
HPDsummary(soft, dd) 
dumm

## Analyze WT data with cq2counts function and naive model ----
ddwt <- cq2counts(data=WT, genecols=c(8:19), condcols=c(1:7), effic=eff)
head(ddwt)

## use soft norm model to fit
soft_wt <- mcmc.qpcr(
  data=ddwt,
  fixed="APA",random="sample",
  controls=c("rRNA18S"),
  normalize = TRUE,
  pr=T,pl=T)

diagnostic.mcmc(model=soft, col="grey50", cex=0.8)
HPDsummary(soft_wt, ddwt, relative=TRUE)
HPDsummary(soft_wt, ddwt)


## Analyze 3 gene data with cq2counts function and naive model ----
dd3genes <- cq2counts(data=CA1_3genes, genecols=c(8:10), condcols=c(1:7), effic=eff)
head(dd3genes)
