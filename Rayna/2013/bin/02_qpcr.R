# Part 2: Reading and analyzing qPCR data

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
dilutions$gene <- revalue(dilutions$gene, c("dlg4.conc" = "dlg4")) 
dilutions$gene <- revalue(dilutions$gene, c("grim.conc" = "grim"))
dilutions$gene <- revalue(dilutions$gene, c("fmr1.conc" = "fmr1"))
dilutions$gene <- revalue(dilutions$gene, c("rRNA18s" = "rRNA18S"))
dilutions <- droplevels(dilutions)

PrimEff(dilutions) -> eff

## Create "all but homecage" dataframe, anlayze with cq2counts function and naive model ----
nohomecage <- filter(qpcr, APA != "homecage")
nohomecage <- droplevels(nohomecage)

dd_nohomecage <- cq2counts(data=nohomecage, genecols=c(9:20), condcols=c(1:8), effic=eff)
head(dd_nohomecage)

naive_dd_nohomecage <- mcmc.qpcr(
  data=dd_nohomecage,
  fixed="APA+region.genotype+APA:region.genotype",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_dd_nohomecage, col="grey50", cex=0.8)
HPDsummary(naive_dd_nohomecage, dd_nohomecage) -> summarynohomecage
trellisByGene(summarynohomecage,xFactor="region.genotype",groupFactor="APA", nrow=4)+xlab(NULL) 
#saved as 12genes-3x4.png

dd_nohomecage %>% 
  ggplot(aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene, scales = "free_y")

dd_nohomecage %>% filter(gene %in% c("rpl19", "grim", "pkmz")) %>%
  ggplot(aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene, scales = "free_y")

dd_nohomecage %>% filter(gene %in% c("rpl19", "rRNA18S", "cam2kd")) %>%
  ggplot(aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene)

dd_nohomecage %>% filter(gene %in% c("fmr1", "fos", "nsf", "pkmz", "creb", "dlg4")) %>%
  ggplot(aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene)

dd_nohomecage %>% filter(gene %in% c("gria", "grim", "grin")) %>%
  ggplot(aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene)

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
  data=dd_FMR1KO,
  fixed="APA",
  pr=T,pl=T, geneSpecRes=FALSE
  )
diagnostic.mcmc(model=naive_FMR1, col="grey50", cex=0.8)
HPDsummary(naive_FMR1, dd_FMR1KO) -> summaryFMR1  

#get the normalize dataframes (using getNormalizedData funct.) and combine them into one 
nd_FMR1 <- getNormalizedData(naive_FMR1,data=dd_FMR1KO) #export results
cbind(nd_FMR1$conditions, nd_FMR1$normData) -> nd_FMR1

# melt the data and rename column to 
library(reshape2)
nd_FMR1 <- melt(nd_FMR1, id.vars = c("ind", "time", "region", "APA", "genotype", "year", "region.genotype"))
names(nd_FMR1)[8] <- "gene"

nd_FMR1$gene <- factor(nd_FMR1$gene, levels = c("cam2kd", "creb", "dlg4", "fmr1", "fos", "gria", "grim", "grin", "nsf", "pkmz", "rpl19", "rRNA18S"))
head(nd_FMR1)

ggplot(dd_FMR1KO, aes(x=genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) + 
  scale_y_log10() +
  facet_wrap(~gene)

## Subset WT-CA3-only data with cq2counts function and naive model ----
WT <- filter(qpcr, genotype == "WT", region == "CA3")
WT <- droplevels(WT)

ddwt <- cq2counts(data=WT, genecols=c(9:20), condcols=c(1:8), effic=eff)

naive_wt <- mcmc.qpcr(
  data=ddwt,
  fixed="APA",random="sample",
  pr=T,pl=T, geneSpecRes=FALSE)
diagnostic.mcmc(model=soft, col="grey50", cex=0.8)
HPDsummary(naive_wt, ddwt, relative=TRUE)

#get the normalize dataframes (using getNormalizedData funct.) and combine them into one 
nd_naive_wt <- getNormalizedData(naive_wt,data=ddwt) #export results
cbind(nd_naive_wt$conditions, nd_naive_wt$normData) -> nd_naive_wt

# melt the data and rename column to 
library(reshape2)
nd_naive_wt <- melt(nd_naive_wt, 
                    id.vars = c("ind", "time", "region", "APA", "genotype", "year", "region.genotype")
)
names(nd_naive_wt)[8] <- "gene"

nd_naive_wt$gene <- factor(nd_naive_wt$gene, 
                           levels = c("cam2kd", "creb", "dlg4", "fmr1", 
                                      "fos", "gria", "grim", "grin",
                                       "nsf", "pkmz", "rpl19", "rRNA18S")
)
head(nd_naive_wt)

ggplot(nd_naive_wt, aes(x=region.genotype, y=value)) + 
  geom_boxplot(aes(fill=APA)) + 
  facet_wrap(~gene)

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
  pr=T,pl=T, singular.ok=TRUE,  geneSpecRes=FALSE)
diagnostic.mcmc(model=naive_3genes, col="grey50", cex=0.8)
HPDsummary(naive_3genes, dd3genes) -> sumary3genes
trellisByGene(sumary3genes,xFactor="region.genotype",groupFactor="APA", nrow=3)+xlab(NULL)
#saved as 3genes-1x3.png

#get the normalize dataframes (using getNormalizedData funct.) and combine them into one 
nd_naive_3genes <- getNormalizedData(naive_3genes,data=dd3genes) #export results
cbind(nd_naive_3genes$conditions, nd_naive_3genes$normData) -> nd_naive_3genes

# melt the data and rename column to 
library(reshape2)
nd_naive_3genes <- melt(nd_naive_3genes, 
                        id.vars = c("ind", "time", "region", "APA", "genotype", "year", "region.genotype")
                        )
names(nd_naive_3genes)[8] <- "gene"
head(nd_naive_3genes)
str(nd_naive_3genes)

library(ggplot2)
ggplot(dd3genes, aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene, scales = "free_y")

dd3genes %>% filter(gene != "rpl19") %>% 
  ggplot(aes(x=region.genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10() +
  facet_wrap(~gene)