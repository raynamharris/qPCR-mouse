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

## rename some columsn and factors

# install.packages("dplyr") #note, this function masked when plyr installed. have to restart to rerun command
library(dplyr)
qpcr <- rename(qpcr, genotype = strained) 

library(plyr)
qpcr$genotype <- revalue(qpcr$genotype, c("fmr1" = "FMR1-KO")) 
qpcr$genotype <- revalue(qpcr$genotype, c("wt" = "WT"))
head(qpcr)

qpcr$region.genotype <- as.factor(paste(qpcr$region, qpcr$genotype, sep="_"))
qpcr <- qpcr[c(1:7,20,8:19)]
head(qpcr)

## calculating gene efficiencies & rename genes ----
#install.packages("MCMC.qpcr")
library(MCMC.qpcr)

dilutions <- read.csv("02_dilutions_CA1CA3.csv", header = TRUE)
head(dilutions)

dilutions <- filter(dilutions, gene != "dlg4", gene != "pkmz.conc", gene != "rpl19.conc", gene != "grim")
dilutions <- droplevels(dilutions)
dilutions$gene <- revalue(dilutions$gene, c("dlg4.conc" = "dlg4")) 
dilutions$gene <- revalue(dilutions$gene, c("grim.conc" = "grim"))
dilutions$gene <- revalue(dilutions$gene, c("fmr1.conc" = "fmr1"))
dilutions$gene <- revalue(dilutions$gene, c("rRNA18s" = "rRNA18S"))

PrimEff(dilutions) -> eff

## Create "all but homecage" dataframe, anlayze with cq2counts function and naive model ----
nohome <- filter(qpcr, APA != "homecage")
nohome <- droplevels(nohome)

ddnohome <- cq2counts(data=nohome, genecols=c(9:20), condcols=c(1:8), effic=eff)
head(ddnohome)

naive_nohome <- mcmc.qpcr(
  data=ddnohome,
  fixed="APA+region.genotype+APA:region.genotype",
  pr=T,pl=T, singular.ok=TRUE,  geneSpecRes=FALSE)
diagnostic.mcmc(model=naive_nohome, col="grey50", cex=0.8)
HPDsummary(naive_nohome, nohome) -> summarynohome
trellisByGene(summarynohome,xFactor="region.genotype",groupFactor="APA", nrow=4)+xlab("Group") 
#printed and saved as 12genes-3x4.png

## some extra analyzes
nd <- getNormalizedData(naive_nohome,data=ddnohome) #export results
nd_normData <- as.data.frame(nd$normData)
nd_conditions <- as.data.frame(nd$conditions)

## Create "all CA1 but no homecage" dataframe, anlayze with cq2counts function and naive model ----
nohomeCA1 <- filter(qpcr, APA != "homecage", region != "CA3")
nohomeCA1 <- droplevels(nohomeCA1)

ddnohomeCA1 <- cq2counts(data=nohomeCA1, genecols=c(9:20), condcols=c(1:8), effic=eff)
head(ddnohome)

naive_nohomeCA1 <- mcmc.qpcr(
  data=ddnohomeCA1,
  fixed="APA+region.genotype+APA:region.genotype",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_nohomeCA1, col="grey50", cex=0.8)
HPDsummary(naive_nohomeCA1, ddnohomeCA1) -> summarynohomeCA1
trellisByGene(summarynohomeCA1,xFactor="region.genotype",groupFactor="APA")+xlab("group")

## Subset FMR1 data then anlyze with cq2counts function and naive model ----
FMR1KO <- filter(qpcr, genotype == "FMR1-KO")
FMR1KO <- droplevels(FMR1KO)
str(FMR1KO)

dd_FMR1KO <- cq2counts(data=FMR1KO, genecols=c(9:20), condcols=c(1:8), effic=eff)
head(dd_FMR1KO)

naive_FMR1 <- mcmc.qpcr(
  data=dd,
  fixed="APA",
  pr=T,pl=T, geneSpecRes=FALSE
  )
diagnostic.mcmc(model=naive_FMR1, col="grey50", cex=0.8)
HPDsummary(naive_FMR1, dd_FMR1KO) -> summaryFMR1  

nd_FMR1 <- getNormalizedData(naive_FMR1,data=dd_FMR1KO) #export results
head(nd_FMR1$normData)
tail(nd_FMR1$conditions)

par(mfrow=c(3,4))
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$cam2kd)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$creb)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$dlg4)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$fmr1)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$fos)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$gria)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$grim)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$grin)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$nsf)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$pkmz)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$rpl19)
plot(nd_FMR1$conditions$APA, nd_FMR1$normData$rRNA18S)

## Subset WT-CA3 data with cq2counts function and naive model ----
WT <- filter(qpcr, genotype == "WT", region == "CA3")
WT <- droplevels(WT)

ddwt <- cq2counts(data=WT, genecols=c(9:20), condcols=c(1:8), effic=eff)

soft_wt <- mcmc.qpcr(
  data=ddwt,
  fixed="APA",random="sample",
  controls=c("rRNA18S"),
  normalize = TRUE,
  pr=T,pl=T)
diagnostic.mcmc(model=soft, col="grey50", cex=0.8)
HPDsummary(soft_wt, ddwt)

## subset WT-CA1-only then analyze data with cq2counts function and naive model ----
CA1_year <- qpcr[c(1:11)] %>% 
  filter(APA != "homecage") %>% 
  filter(genotype != "FMR1-KO") %>% 
  filter(region == "CA1") 
CA1_year <- droplevels(CA1_year)
str(CA1_year)

ddCA1year <- cq2counts(data=CA1_year, genecols=c(9:11), condcols=c(1:8), effic=eff)

naive_ddCA1year <- mcmc.qpcr(
  data=ddCA1year,
  fixed="year+APA+APA:year",random="sample",
  pr=T,pl=T)
diagnostic.mcmc(model=naive_ddCA1year, col="grey50", cex=0.8)
HPDsummary(naive_ddCA1year, ddCA1year)

## Suset 3 gene data then anlyzewith cq2counts function and naive model ----
CA1_3genes <- qpcr[c(1:11)] %>% filter( APA != "homecage")
CA1_3genes <- droplevels(CA1_3genes)
str(CA1_3genes)

dd3genes <- cq2counts(data=CA1_3genes, genecols=c(9:11), condcols=c(1:8), effic=eff)

naive_3genes <- mcmc.qpcr(
  data=dd3genes,
  fixed="region.genotype+APA+APA:region.genotype",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_3genes, col="grey50", cex=0.8)
HPDsummary(naive_3genes, dd3genes) -> sumary3genes
trellisByGene(sumary3genes,xFactor="region.genotype",groupFactor="APA")+xlab("group")

