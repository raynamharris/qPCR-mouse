# Part 2: Reading and analyzing qPCR data

## load libraries ----
library(dplyr) # for renaming columns
library(plyr) # for renaming factors
library(MCMC.qpcr) # for qpcr analysis
library(reshape2) # for melting data
library(reshape) #for making data wide
library(cowplot) #for multiple plots
library(gplots) # for heatmap.2
library(corrplot) # simple corrplot with size-circles
library(PerformanceAnalytics) # fancy ass correlation plot
library(Hmisc) # for correlation stats
library(cowplot) # for multip pannel plots

## wrangle the gene expression qpcr data ----
setwd("~/Github/qPCR-mouse/Rayna/2013/data")
qpcr <- read.csv("02_qpcrdata.csv", header = TRUE, na.strings = "NA", stringsAsFactors = FALSE)
str(qpcr)
summary(qpcr)
head(qpcr)

## setting the factors
qpcr$sample <- as.factor(qpcr$sample)
qpcr$ind <- as.factor(qpcr$ind)
qpcr$time <- as.factor(qpcr$time)
qpcr$region <- as.factor(qpcr$region)
qpcr$APA <- as.factor(qpcr$APA)
qpcr$strain <- as.factor(qpcr$strain)
qpcr$X <- NULL
str(qpcr)
summary(qpcr)
head(qpcr)

## rename some columsn and factors
qpcr <- dplyr::rename(qpcr, genotype = strain) 
qpcr$genotype <- revalue(qpcr$genotype, c("fmr1" = "FMR1-KO")) 
qpcr$genotype <- revalue(qpcr$genotype, c("wt" = "WT"))
names(qpcr)[7] <- "prkcz"


## rename "WT AB" to "FMR1-KO AB"
head(qpcr, 12)
qpcr$ind <- revalue(qpcr$ind, c("BL AB" = "FMR1 AB"))
qpcr[10:12, 6] = "FMR1-KO"

## create new columns to join factors
qpcr$region.genotype <- as.factor(paste(qpcr$region, qpcr$genotype, sep="_"))
qpcr$genoAPA <- as.factor(paste(qpcr$genotype,qpcr$APA, sep="_"))
qpcr$genoAPAregion <- as.factor(paste(qpcr$genoAPA,qpcr$region, sep="_"))
names(qpcr)

## reorder the dataframe
qpcr <- qpcr[c(1:6,19:21,7:18)]
names(qpcr)

## calculating gene efficiencies & rename genes ----
dilutions <- read.csv("02_dilutions_CA1CA3.csv", header = TRUE)
head(dilutions)

dilutions <- filter(dilutions, gene != "dlg4", gene != "pkmz.conc", gene != "rpl19.conc", gene != "grim")
dilutions$gene <- revalue(dilutions$gene, c("dlg4.conc" = "dlg4")) 
dilutions$gene <- revalue(dilutions$gene, c("grim.conc" = "grim"))
dilutions$gene <- revalue(dilutions$gene, c("fmr1.conc" = "fmr1"))
dilutions$gene <- revalue(dilutions$gene, c("rRNA18s" = "rRNA18S"))
dilutions$gene <- revalue(dilutions$gene, c("pkmz" = "prkcz"))

dilutions <- droplevels(dilutions)

PrimEff(dilutions) -> eff

## Analyze just FMR1 data ----
head(qpcr)
FMR1 <- dplyr::filter(qpcr, APA != "home", genotype != "WT")
FMR1 <- droplevels(FMR1)
str(FMR1)

dd_FMR1 <- cq2counts(data=FMR1, genecols=c(10:21), condcols=c(1:9), effic=eff)
head(dd_FMR1)

naive_dd_FMR1 <- mcmc.qpcr(
  data=dd_FMR1,
  fixed="APA",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_dd_FMR1, col="grey50", cex=0.8)
HPDsummary(naive_dd_FMR1, dd_FMR1) -> summaryFMR1
summary(naive_dd_FMR1)

# saved as 2-FMR1CA1qpcrdata.png
FMFR1Palette <- c('grey50','darkorange')
dd_FMR1 %>%  
  ggplot(aes(x=APA, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10(name="Gene Expression (Log10 Counts)") +
  facet_wrap(~gene, scales = "free_y") +
  scale_x_discrete(name="") +
  theme(legend.position="none", 
        strip.text = element_text(face = "italic")) + 
  scale_fill_manual(values = FMFR1Palette)

## Correlations

# first widen the data
head(dd_FMR1) 
dd_FMR1_wide <- dcast(dd_FMR1, ind + APA ~ gene, value.var= "count", fun.aggregate=mean)

## then create a matrix of just gene expression data with NAs ommited
dd_FMR1_wide_matrix <- dd_FMR1_wide # prepare for matrix
dd_FMR1_wide_matrix$APAind <- as.factor(paste(dd_FMR1_wide_matrix$APA,dd_FMR1_wide_matrix$ind, sep="_"))
rownames(dd_FMR1_wide_matrix) <- dd_FMR1_wide_matrix$APAind  # set $genoAPAsession as rownames
names(dd_FMR1_wide_matrix)
dd_FMR1_wide_matrix <- dd_FMR1_wide_matrix[-c(1:2,15)]  ## remove non-numeric columns
dd_FMR1_wide_matrix <- as.matrix(dd_FMR1_wide_matrix)
dd_FMR1_wide_matrix <- na.omit(dd_FMR1_wide_matrix)
head(dd_FMR1_wide_matrix)

# then create a correlation matrix
dd_FMR1_wide_matrix_cor <- round(cor(dd_FMR1_wide_matrix),2) 
head(dd_FMR1_wide_matrix_cor,12)

# plot the correlation matrix
## saved as 2-FMR1qpcrcorrelationmatrix.png
dev.off()
corrplot(dd_FMR1_wide_matrix_cor, type="lower", order="hclust", tl.col="black", tl.srt=45)

# ploting one correlation at a time by group
a1 <-  ggplot(dd_FMR1_wide, aes(x = dlg4, y = fos, colour = APA)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) + 
  scale_colour_manual(values = FMFR1Palette) +
  theme(legend.position="none")
b1 <- ggplot(dd_FMR1_wide, aes(x = prkcz, y = gria, colour = APA)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) + 
  scale_colour_manual(values = FMFR1Palette) +
  theme(legend.position="none")
c1 <- ggplot(dd_FMR1_wide, aes(x = gria, y = grin, colour = APA)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) + 
  scale_colour_manual(values = FMFR1Palette) +
  theme(legend.position="none")
d1 <- ggplot(dd_FMR1_wide, aes(x = cam2kd, y = nsf, colour = APA)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) + 
  scale_colour_manual(values = FMFR1Palette) +
  theme(legend.position="none")

a2 <- ggplot(dd_FMR1_wide, aes(x = dlg4, y = fos)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) 
b2 <- ggplot(dd_FMR1_wide, aes(x = gria, y = prkcz)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) + 
  scale_colour_manual(values = FMFR1Palette)
c2 <- ggplot(dd_FMR1_wide, aes(x = gria, y = grin)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) 
d2 <- ggplot(dd_FMR1_wide, aes(x = cam2kd, y = nsf)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm) 

plot_grid(a1, b1, c1, d1, a2, b2, c2, d2, nrow = 2)


## Analyzing just the WT CA3 data ----
head(qpcr)
CA3only <- dplyr::filter(qpcr, region == "CA3")
CA3only <- droplevels(CA3only)
names(CA3only)

dd_CA3only <- cq2counts(data=CA3only, genecols=c(10:21), condcols=c(1:9), effic=eff)
head(dd_CA3only)

naive_dd_CA3only <- mcmc.qpcr(
  data=dd_CA3only,
  fixed="APA",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_dd_CA3only, col="grey50", cex=0.8)
HPDsummary(naive_dd_CA3only, dd_CA3only) -> summaryCA3
summary(naive_dd_CA3only)

#saved as 2-WTCA3geneexpression.png
WTPalette <- c('black','red')
dd_CA3only %>%  
  ggplot(aes(x=APA, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10(name="Gene Expression (Log10 Counts)") +
  facet_wrap(~gene, scales = "free_y") +
  scale_x_discrete(name="") +
  theme(legend.position="none", 
        strip.text = element_text(face = "italic")) + 
  scale_fill_manual(values = WTPalette)



### merging the FMR1 CA1 and WT CA3 datasets
names(dd_CA3only)
names(dd_FMR1)
dds <-  dplyr::bind_rows(dd_FMR1, dd_CA3only)
str(dds)
str(dd_FMR1)

dds$sample <- as.factor(dds$sample)
dds$ind <- as.factor(dds$ind)
dds$time <- as.factor(dds$time)
dds$region <- as.factor(dds$region)
dds$APA <- as.factor(dds$APA)
dds$genotype <- as.factor(dds$genotype)
dds$region.genotype <- as.factor(dds$region.genotype)
dds$genoAPA <- as.factor(dds$genoAPA)
dds$genoAPAregion <- as.factor(dds$genoAPAregion)
dds$genoAPA <- factor(dds$genoAPA, levels = c("WT_control", "WT_trained", "FMR1-KO_control", "FMR1-KO_trained"))
str(dds)

## saved as 2-WTCA3_FMR1CA1data
FentonPalette2 <- c('black','red','grey50','darkorange')
dds %>%  
  ggplot(aes(x=genoAPA, y=count, fill= genoAPA)) + 
  geom_boxplot() +  scale_y_log10(name="Gene Expression (Log10 Counts)") +
  facet_wrap(~gene, scales = "free_y") +
  scale_x_discrete(name="") +
  theme(legend.position="bottom", 
        strip.text = element_text(face = "italic"),
        axis.ticks = element_blank(), axis.text.x = element_blank()) + 
  scale_fill_manual(values = FentonPalette2,
                    breaks=c("WT_control", "WT_trained", "FMR1-KO_control", "FMR1-KO_trained"),
                    labels=c("WT CA3 control", "WT CA3 trained", "FMR1-KO CA1 control", "FMR1-KO CA1 trained")) + 
  guides(fill=guide_legend(title=NULL)) #removed legend title


## Suset 3 gene data then anlyzewith cq2counts function and naive model ----
head(qpcr)
3genes <- qpcr[c(1:12)] %>% 
  filter(APA != "home")
CA1_3genes <- droplevels(CA1_3genes)
str(CA1_3genes)

dd3genes <- cq2counts(data=CA1_3genes, genecols=c(10:12), condcols=c(1:9), effic=eff)

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

## correlations 3 genes
dd3genes_wide <- dcast(dd3genes, ind + region.genotype + genoAPA + region ~ gene, value.var= "count", fun.aggregate=mean)
head(dd3genes_wide)

## ploting one correlation at a time
ggplot(dd3genes_wide, aes(x = pkmz, y = grim, colour = region.genotype)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm)

## See http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software for code
## making a matrix to plot many correlations
dd3genes_wide_matrix <- dd3genes_wide
dd3genes_wide_matrix$genoAPAregionInd <- as.factor(paste(dd3genes_wide_matrix$genoAPA,dd3genes_wide_matrix$region, dd3genes_wide_matrix$ind, sep="_"))
rownames(dd3genes_wide_matrix) <- dd3genes_wide_matrix$genoAPAregionInd  # set $genoAPAsession as rownames
names(dd3genes_wide_matrix)
dd3genes_wide_matrix <- dd3genes_wide_matrix[-c(1:4)] 
dd3genes_wide_matrix <- as.matrix(dd3genes_wide_matrix)
head(dd3genes_wide_matrix)
dd3genes_wide_matrix <- na.omit(dd3genes_wide_matrix)
head(dd3genes_wide_matrix)
tail(dd3genes_wide_matrix)


dd3genes_wide_matrix_cor <- round(cor(dd3genes_wide_matrix),2) 
head(dd3genes_wide_matrix_cor)
dev.off()
corrplot(dd3genes_wide_matrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)

## fancy ass correlation matrix
chart.Correlation(dd3genes_wide_matrix_cor, histogram=FALSE, pch=19, method = "spearman")

## saved as 2-heatmap-CA1CA3.png
col <- colorRampPalette(c("turquoise4", "white", "tan4"))(n = 299)
heatmap.2(dd3genes_wide_matrix_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none", 
          cellnote = dd3genes_wide_matrix_cor, notecol="black")
# See https://github.com/rasbt/R_snippets/blob/master/heatmaps/h3_categorizing.R



### useful plot code not used ----