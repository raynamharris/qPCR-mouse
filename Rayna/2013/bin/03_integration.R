## integrative analysis ----

########## gather three dataset and exampine ----
behavior <- wtfmr1
phy <- summary
qpcrcounts <- dds

head(behavior) #wide
head(phy) #wide
head(qpcrcounts) #long

## need to have lenghtened, with make a measure column and a value column that can be used to merge all -----

## Behavior: need to melt, remove some rows, make sesssionbehavior column, then widen
head(behavior)
behavior <- melt(wtfmr1, id=c("ind","genotype", "APA", "session", "genoAPA", "genoAPAsession", "genoAPAsessionInd", "filename"))
head(behavior)
behavior <- filter(behavior, !grepl("p.miss|TotalTime", variable))
behavior$measure <- as.factor(paste(behavior$session, behavior$variable, sep="_"))
behavior <- dplyr::select(behavior, ind, genotype, APA, measure, value)
head(behavior)

## Ephys: need to melt, remove some rows, make sesssionbehavior column, then widen
head(phy)
phy <- melt(summary, id=c("ind","genotype", "APA", "session", "genoAPA", "filename"))
head(phy)
phy <- dplyr::rename(phy, measure = variable)
phy <- dplyr::select(phy, ind, genotype, APA, measure, value)
head(phy)

## qpcr need to rename select variables
qpcrcounts <- dds %>%
  dplyr::rename(measure = gene, value=count) %>%
  dplyr::select(ind, genotype, APA, measure, value)
head(qpcrcounts)

## log transform count data
qpcrcounts$value[qpcrcounts$value == 0] <- NA
qpcrcounts$value <- log10(qpcrcounts$value)
head(qpcrcounts)

## widen to average, then melt
qpcrcounts <- dcast(qpcrcounts, ind + genotype + APA ~ measure, fun=mean)
qpcrcountslong <- melt(qpcrcounts, id=c("ind","genotype", "APA"))
qpcrcountslong <- dplyr::rename(qpcrcountslong, measure = variable)
head(qpcrcountslong)
str(qpcrcountslong)

## remove gens with >5 NA
qpcrcountslong <- qpcrcountslong %>% 
  filter(!grepl("creb|fos|fmr1|grim", measure))  
qpcrcountslong
  
############# merge all and widen!! -----

all <- dplyr::bind_rows(phy, qpcrcountslong)
all$ind <- as.factor(all$ind)
all$measure <- as.factor(all$measure)
all$APA <- as.factor(all$APA)
all$APA <- revalue(all$APA, c("control" = "untrained")) 
all$genoAPA <- as.factor(paste(all$genotype, all$APA, sep="_"))
all <- dplyr::select(all, ind, genotype, APA, genoAPA, measure, value)
all <- dcast(all, ind + genotype + APA + genoAPA ~ measure, value.var= "value")
str(all)

ggplot(all, aes(x=grin, y=IO_Max, color=genoAPA)) +
  geom_point(size = 8) +
  scale_color_manual(values=colorvalgenoAPA) + 
  geom_smooth(method = lm)

ggplot(all, aes(x=log10(TotalPunishment), y=log10(IO_Max), color=genoAPA)) +
  geom_point(size = 8) +
  scale_color_manual(values=colorvalgenoAPA) +
  geom_smooth(method = lm)

ggplot(all, aes(x=log10(TotalPunishment), y=log10(gria), color=genoAPA)) +
  geom_point(size = 8) +
  scale_color_manual(values=colorvalgenoAPA) +
  geom_smooth(method = lm)


ggplot(all, aes(x=log10(TotalPunishment), y=log10(grin), color=genoAPA)) +
  geom_point(size = 8) +
  scale_color_manual(values=colorvalgenoAPA) +
  geom_smooth(method = lm)

ggplot(all, aes(x=log10(TotalPunishment), y=log10(prkcz), color=genoAPA)) +
  geom_point(size = 8) +
  scale_color_manual(values=colorvalgenoAPA) +
  geom_smooth(method = lm)

all %>%
  filter(APA == "trained") %>%
ggplot(aes(x=log10(IO_Max), y=log10(grin), color=genoAPA)) +
  geom_point(size = 8) +
  scale_color_manual(values=colorvalgenoAPA) 
  geom_smooth(method = lm, alpha = 0.2)

correlation <- all %>%
  #filter(APA == "trained") %>%
  ggplot(aes(x=log10(grin), y=log10(gria), color=genoAPA)) +
  geom_point(size = 2) +
  scale_color_manual(values=colorvalgenoAPA) +
  geom_smooth(method = lm, alpha = 0.2) +
  theme_cowplot(font_size = 8, line_size = 0.25) +
  theme(legend.position="none",
        axis.title = element_blank()) 
correlation

pdf(file="../figures/3-correlation1.pdf", width=1.5, height=2)
plot(correlation)
dev.off()

correlation <- all %>%
  #filter(APA == "trained") %>%
  ggplot(aes(x=log10(grin), y=log10(IO_Max), color=genoAPA)) +
  geom_point(size = 2) +
  scale_color_manual(values=colorvalgenoAPA) +
  geom_smooth(method = lm, alpha = 0.2) +
  theme_cowplot(font_size = 8, line_size = 0.25) +
  theme(legend.position="none",
        axis.title = element_blank()) 
  
correlation

pdf(file="../figures/3-correlation2.pdf", width=1.5, height=2)
plot(correlation)
dev.off()
```


```{r correlationmatrix}
## then create a matrix of just gene expression data with NAs ommited
allmatrix <- all # prepare for matrix
allmatrix$genoAPAind <- as.factor(paste(allmatrix$genoAPA,allmatrix$ind, sep="_"))
rownames(allmatrix) <- allmatrix$genoAPAind  # set $genoAPAsession as rownames
names(allmatrix)
allmatrix <- allmatrix[-c(1:6,10,11,13,14,15,17,19,20,21,24)]  ## remove non-numeric columns
allmatrix <- log10(allmatrix + 1) #log transform all data
allmatrix <- as.matrix(allmatrix)
str(allmatrix)
allmatrix <- na.omit(allmatrix)
allmatrix_cor <- round(cor(allmatrix),2)  # then create a correlation matrix

corrplot(allmatrix_cor, type="lower", tl.col="black", tl.srt=45)
```