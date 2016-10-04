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
names(qpcr)[10] <- "prkcz"


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

## Create "no homecage no WT" dataframe, anlayze with cq2counts function and naive model ----
head(qpcr)
nohomecage_nowt <- dplyr::filter(qpcr, APA != "home", genotype != "WT")
nohomecage_nowt <- droplevels(nohomecage_nowt)
str(nohomecage_nowt)

dd_nohomecage_nowt <- cq2counts(data=nohomecage_nowt, genecols=c(10:21), condcols=c(1:9), effic=eff)
head(dd_nohomecage_nowt)

naive_dd_nohomecage_nowt <- mcmc.qpcr(
  data=dd_nohomecage_nowt,
  fixed="APA",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_dd_nohomecage_nowt, col="grey50", cex=0.8)
HPDsummary(naive_dd_nohomecage_nowt, dd_nohomecage_nowt) -> summarynohomecage_nowt

# saved as 2-FMR1CA1qpcrdata.png
FMFR1Palette <- c('grey50','darkorange')
dd_nohomecage_nowt %>%  
  ggplot(aes(x=APA, y=count)) + 
  geom_boxplot(aes(fill=APA)) +  scale_y_log10(name="Gene Expression (Log10 Counts)") +
  facet_wrap(~gene, scales = "free_y") +
  scale_x_discrete(name="") +
  theme(legend.position="none", 
        strip.text = element_text(face = "italic")) + 
  scale_fill_manual(values = FMFR1Palette)

## Correlations!!! ----
head(dd_nohomecage_nowt)
dd_nohomecage_nowt_wide <- dcast(dd_nohomecage_nowt, ind + region.genotype + genoAPA + region ~ gene, value.var= "count", fun.aggregate=mean)
head(dd_nohomecage_nowt_wide)

## ploting one correlation at a time
ggplot(dd_nohomecage_nowt_wide, aes(x = cam2kd, y = creb, colour = region.genotype)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm)

## See http://www.sthda.com/english/wiki/correlation-matrix-a-quick-start-guide-to-analyze-format-and-visualize-a-correlation-matrix-using-r-software for code
## making a matrix to plot many correlations
dd_nohomecage_nowt_wide_matrix <- dd_nohomecage_nowt_wide
dd_nohomecage_nowt_wide_matrix$genoAPAregionInd <- as.factor(paste(dd_nohomecage_nowt_wide_matrix$genoAPA,dd_nohomecage_nowt_wide_matrix$region, dd_nohomecage_nowt_wide_matrix$ind, sep="_"))
rownames(dd_nohomecage_nowt_wide_matrix) <- dd_nohomecage_nowt_wide_matrix$genoAPAregionInd  # set $genoAPAsession as rownames
names(dd_nohomecage_nowt_wide_matrix)
dd_nohomecage_nowt_wide_matrix <- dd_nohomecage_nowt_wide_matrix[-c(1:4,17)] 
dd_nohomecage_nowt_wide_matrix <- as.matrix(dd_nohomecage_nowt_wide_matrix)
head(dd_nohomecage_nowt_wide_matrix)
dd_nohomecage_nowt_wide_matrix <- na.omit(dd_nohomecage_nowt_wide_matrix)
head(dd_nohomecage_nowt_wide_matrix)
tail(dd_nohomecage_nowt_wide_matrix)

dd_nohomecage_nowt_wide_matrix_cor <- round(cor(dd_nohomecage_nowt_wide_matrix),2) 
head(dd_nohomecage_nowt_wide_matrix_cor)
dev.off()
corrplot(dd_nohomecage_nowt_wide_matrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)

chart.Correlation(dd_nohomecage_nowt_wide_matrix_cor, histogram=FALSE, pch=19, method = "spearman")

## saved as 2-heatmap-CA1CA3.png
col <- colorRampPalette(c("turquoise4", "white", "tan4"))(n = 299)
heatmap.2(dd_nohomecage_nowt_wide_matrix_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none", 
          cellnote = dd_nohomecage_nowt_wide_matrix_cor, notecol="black")
# See https://github.com/rasbt/R_snippets/blob/master/heatmaps/h3_categorizing.R














## CA1 specific heatmap
head(dd_nohomecage_wide)
CA1onlymatrix <- dd_nohomecage_wide %>% 
  filter(region == "CA1")
names(CA1onlymatrix)
CA1onlymatrix$genoAPAregionInd <- as.factor(paste(CA1onlymatrix$genoAPA,CA1onlymatrix$region, CA1onlymatrix$ind, sep="_"))
rownames(CA1onlymatrix) <- CA1onlymatrix$genoAPAregionInd  # set $genoAPAsession as rownames
names(CA1onlymatrix)
CA1onlymatrix <- CA1onlymatrix[-c(1:4,17)] 
CA1onlymatrix <- as.matrix(CA1onlymatrix)
head(CA1onlymatrix)
CA1onlymatrix <- na.omit(CA1onlymatrix)
head(CA1onlymatrix)
tail(CA1onlymatrix)

CA1onlymatrix_cor <- round(cor(CA1onlymatrix),2) 
head(CA1onlymatrix_cor)
corrplot(CA1onlymatrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)

heatmap.2(CA1onlymatrix_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none", 
          cellnote = dd_nohomecage_wide_matrix_cor, notecol="black",
          main = "CA1") 

chart.Correlation(CA1onlymatrix_cor, histogram=TRUE, pch=19, method = "spearman")


## CA3 specific heatmap
head(dd_nohomecage_wide)
CA3onlymatrix <- dd_nohomecage_wide %>% 
  filter(region == "CA3")
names(CA3onlymatrix)
CA3onlymatrix$genoAPAregionInd <- as.factor(paste(CA3onlymatrix$genoAPA,CA3onlymatrix$region, CA3onlymatrix$ind, sep="_"))
rownames(CA3onlymatrix) <- CA3onlymatrix$genoAPAregionInd  # set $genoAPAsession as rownames
names(CA3onlymatrix)
CA3onlymatrix <- CA3onlymatrix[-c(1:4,17)] 
CA3onlymatrix <- as.matrix(CA3onlymatrix)
head(CA3onlymatrix)
CA3onlymatrix <- na.omit(CA3onlymatrix)
head(CA3onlymatrix)
tail(CA3onlymatrix)

CA3onlymatrix_cor <- round(cor(CA3onlymatrix),2) 
head(CA3onlymatrix_cor)
corrplot(CA3onlymatrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)

heatmap.2(CA3onlymatrix_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none", 
          cellnote = dd_nohomecage_wide_matrix_cor, notecol="black",
          main = "CA3") 

chart.Correlation(CA3onlymatrix_cor, histogram=FALSE, pch=19, method = "spearman")


## Create "all CA1 but no homecage" dataframe, anlayze with cq2counts function and naive model ----
CA1nohomecage <- filter(qpcr, APA != "home", region != "CA3")
CA1nohomecage <- droplevels(CA1nohomecage)

ddCA1nohomecage <- cq2counts(data=CA1nohomecage, genecols=c(10:21), condcols=c(1:9), effic=eff)
head(ddCA1nohomecage)

naive_CA1nohomecage <- mcmc.qpcr(
  data=ddCA1nohomecage,
  fixed="APA+region.genotype+APA:region.genotype",
  pr=T,pl=T, singular.ok=TRUE)
diagnostic.mcmc(model=naive_CA1nohomecage, col="grey50", cex=0.8)
HPDsummary(naive_CA1nohomecage, ddCA1nohomecage) -> summaryCA1nohomecage
trellisByGene(summaryCA1nohomecage,xFactor="region.genotype",groupFactor="APA")+xlab("group")


## Subset FMR1 data then anlyze with cq2counts function and naive model ----
FMR1KO <- filter(qpcr, genotype == "FMR1-KO")
FMR1KO <- droplevels(FMR1KO)
str(FMR1KO)

dd_FMR1KO <- cq2counts(data=FMR1KO, genecols=c(10:21), condcols=c(1:9), effic=eff)
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
nd_FMR1 <- melt(nd_FMR1, id.vars = c("ind", "time", "region", "APA", "genotype", "region.genotype", "genoAPA", "genoAPAregion"))
names(nd_FMR1)

# rename and relevel the genes
nd_FMR1 <- dplyr::rename(nd_FMR1, gene = variable) 
nd_FMR1$gene <- factor(nd_FMR1$gene, levels = c("cam2kd", "creb", "dlg4", "fmr1", "fos", "gria", "grim", "grin", "nsf", "pkmz", "rpl19", "rRNA18S"))
head(nd_FMR1)

ggplot(dd_FMR1KO, aes(x=genotype, y=count)) + 
  geom_boxplot(aes(fill=APA)) + 
  scale_y_log10() +
  facet_wrap(~gene)

## Subset WT-CA3-only data with cq2counts function and naive model ----
WT <- filter(qpcr, genotype == "WT", region == "CA3")
WT <- droplevels(WT)

ddwt <- cq2counts(data=WT, genecols=c(10:21), condcols=c(1:9), effic=eff)

naive_wt <- mcmc.qpcr(
  data=ddwt,
  fixed="APA",random="sample",
  pr=T,pl=T, geneSpecRes=FALSE)
diagnostic.mcmc(model=soft, col="grey50", cex=0.8)
HPDsummary(naive_wt, ddwt, relative=TRUE)

#get the normalize dataframes (using getNormalizedData funct.) and combine them into one 
nd_naive_wt <- getNormalizedData(naive_wt,data=ddwt) #export results
cbind(nd_naive_wt$conditions, nd_naive_wt$normData) -> nd_naive_wt
names(nd_naive_wt)

# melt the data and rename column to 
nd_naive_wt <- melt(nd_naive_wt, 
                    id.vars = c("ind", "time", "region", "APA", "genotype", "region.genotype", "genoAPA", "genoAPAregion")
)
# rename and relevel the genes
nd_naive_wt <- dplyr::rename(nd_naive_wt, gene = variable) 

nd_naive_wt$gene <- factor(nd_naive_wt$gene, 
                           levels = c("cam2kd", "creb", "dlg4", "fmr1", 
                                      "fos", "gria", "grim", "grin",
                                       "nsf", "pkmz", "rpl19", "rRNA18S")
)
head(nd_naive_wt)

ggplot(nd_naive_wt, aes(x=region.genotype, y=value)) + 
  geom_boxplot(aes(fill=APA)) + 
  facet_wrap(~gene)


## Suset 3 gene data then anlyzewith cq2counts function and naive model ----
CA1_3genes <- qpcr[c(1:12)] %>% filter( APA != "home")
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