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

all <- dplyr::bind_rows(behavior, phy, qpcrcountslong)
head(all)
str(all)
all$ind <- as.factor(all$ind)
all$measure <- as.factor(all$measure)
all$genoAPA <- as.factor(paste(all$genotype, all$APA, sep="_"))
all <- dplyr::select(all, ind, genotype, APA, genoAPA, measure, value)
head(all)
all <- dcast(all, ind + genotype + APA +genoAPA ~ measure, value.var= "value")

head(all)
names(all)