## Want one script to wrangle and anlayze all the wt frm1 data for the qpcr project

# Parts
## Part 1 : behavior and physiology data
## Part 2 : qpcr data
## Part 3 : integraive analysis



# Reading and analyzing qPCR data


## wrangle the gene expression qpcr data ----
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
qpcr <- rename(qpcr, genotype = strained)

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

## read in dilution series and calculate gene effeciencies
dilutions <- read.csv("02_dilutions_CA1CA3.csv", header = TRUE)
str(dilutions)
head(dilutions)
PrimEff(dilutions)

## dropping the poorer std curves and renameing genes. redo primeff()
dilutions <- filter(dilutions, gene != "dlg4", gene != "pkmz.conc", gene != "rpl19.conc", gene != "grim")
dilutions <- droplevels(dilutions)
dilutions$gene <- revalue(dilutions$gene, c("dlg4.conc" = "dlg4")) 
dilutions$gene <- revalue(dilutions$gene, c("grim.conc" = "grim"))
dilutions$gene <- revalue(dilutions$gene, c("fmr1.conc" = "fmr1"))
dilutions$gene <- revalue(dilutions$gene, c("rRNA18s" = "rRNA18S"))
PrimEff(dilutions) -> eff


## Analyze FMR1 data with cq2counts function
dd <- cq2counts(data=FMR1KO, genecols=c(7:18), condcols=c(1:6), effic=eff)
head(dd)

## use naive model to fit
naive <- mcmc.qpcr(
  fixed="APA",
  random="ind",
  data=dd,
  pr=T,
  pl=T
)

diagnostic.mcmc(
  model=naive,
  col="grey50",
  cex=0.8
)

HPDsummary(naive, dd, relative=TRUE)

## use soft norm model to fit
soft <- mcmc.qpcr(
  fixed="APA",
  random="ind",
  data=dd,
  controls=c("rRNA18S","rpl19"),
  normalize = TRUE,
  pr=T,
  pl=T
)

diagnostic.mcmc(
  model=soft,
  col="grey50",
  cex=0.8
)

HPDsummary(soft, dd, relative=TRUE)
