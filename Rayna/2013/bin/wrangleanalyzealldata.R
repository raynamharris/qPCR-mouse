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
FMR1KO <- filter(qpcr, genotype == "FMR1-KO", APA != "homecage")
FMR1KO <- droplevels(FMR1KO)
str(FMR1KO)

WT <- filter(qpcr, genotype == "WT", region == "CA3")
WT <- droplevels(WT)
str(WT)


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
dd <- cq2counts(data=FMR1KO, genecols=c(7:18), condcols=c(1:6), effic=eff)
head(dd)

## use naive model to fit
naive <- mcmc.qpcr(
  data=dd,
  fixed="APA",random="ind",
  pr=T, pl=T)

diagnostic.mcmc(model=naive, col="grey50", cex=0.8)
HPDsummary(naive, dd, relative=TRUE)

## use soft norm model to fit
soft <- mcmc.qpcr(
  data=dd,
  fixed="APA",random="ind",
  controls=c("rRNA18S"),
  normalize = TRUE,
  pr=T,pl=T)

diagnostic.mcmc(model=soft, col="grey50", cex=0.8)
HPDsummary(soft, dd, relative=TRUE)
HPDsummary(soft, dd)


## Analyze WT data with cq2counts function and naive model ----
ddwt <- cq2counts(data=WT, genecols=c(7:18), condcols=c(1:6), effic=eff)
head(ddwt)

## use naive model to fit
naive_wt <- mcmc.qpcr(
  data=ddwt,
  fixed="APA",random="sample",
  pr=T, pl=T)

diagnostic.mcmc(model=naive_wt, col="grey50", cex=0.8)
HPDsummary(naive_wt, ddwt, relative=TRUE)
HPDsummary(naive_wt, ddwt)

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


## Analyze All data (minus that one CA1 WT sample) with cq2counts function and naive model ----
## combine the WT and FRM1 dataset because some very tiny groups have been removed from qpcr dataframe
FMR1KO_WT <- rbind(FMR1KO, WT)
FMR1KO_WT$genotype <- factor(FMR1KO_WT$genotype, levels = c("WT", "FMR1-KO"))

## cq to counds and 
fmr1ko_wt <- cq2counts(data=FMR1KO_WT, genecols=c(7:18), condcols=c(1:6), effic=eff)
head(fmr1ko_wt)
fmr1ko_wt$genotype <- factor(fmr1ko_wt$genotype, levels = c("WT", "FMR1-KO"))
levels(fmr1ko_wt$genotype)

## use naive model to fit
naive_fmr1ko_wt <- mcmc.qpcr(
  data=fmr1ko_wt,
  fixed="APA+genotype+APA:genotype",
  pr=T, pl=T)

diagnostic.mcmc(model=naive_fmr1ko_wt, col="grey50", cex=0.8)
HPDsummary(naive_fmr1ko_wt, fmr1ko_wt, relative=TRUE)
HPDsummary(naive_fmr1ko_wt, fmr1ko_wt)


## use soft norm model to fit
soft_fmr1ko_wt <- mcmc.qpcr(
  data=fmr1ko_wt,
  fixed="APA+genotype+APA:genotype",
  controls = c("rRNA18S"),
  pr=T, pl=T)

diagnostic.mcmc(model=soft_fmr1ko_wt, col="grey50", cex=0.8)
HPDsummary(soft_fmr1ko_wt, fmr1ko_wt, relative=TRUE)
HPDsummary(soft_fmr1ko_wt, fmr1ko_wt)

## use soft norm model to fit
soft_fmr1ko_wt <- mcmc.qpcr(
  data=fmr1ko_wt,
  fixed="region+APA+region:APA",
  controls = c("rRNA18S"),
  pr=T, pl=T)

diagnostic.mcmc(model=soft_fmr1ko_wt, col="grey50", cex=0.8)
HPDsummary(soft_fmr1ko_wt, fmr1ko_wt, relative=TRUE)
sHPDsummary(soft_fmr1ko_wt, fmr1ko_wt)
