## Data Wrangling for Student Projects

## All punches and tissues that Rayna has collected during the  2015 and 2016 NS&B course
## are described in the same csv file in BehavPhysRNAseq repo
## I'll use this script to create a file for the students and save it in the qPCR-mouse repo

library(plyr)
library(dplyr)
library(tidyr)
library(reshape2)

## 1. Read and join data
## 2. clean data
## 3. select samples for exercises
## 4. select huntington's mice
## 5. select fmr1 mice

## 1. read data
setwd("~/Github/qPCR-mouse/data")
setwd("Z:/NSB_2016/4_MouseMolecular/qPCR-mouse/data")
setwd("/Volumes/nsb/NSB_2016/4_MouseMolecular/qPCR-mouse/data")

punches <- read.csv("NSBpunches.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
animals <- read.csv("NSBanimals.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
punches$mouse <- as.factor(punches$Mouse)
animals$mouse <- as.factor(animals$Mouse)


# 2. melt slice data with atlas reference
animalsmelt<- melt(animals,id.vars=c("mouse","Genotype","Conflict","APA","Group"), measure.vars=c("Slice1","Slice2","Slice3","Slice4","Slice5","Slice6","Slice7","Slice8"),variable.name="Slice2",value.name="Atlas.location",na.rm=FALSE)
animalsmelt$Slice2 <- gsub("Slice","",animalsmelt$Slice2)
animalsmelt$Slice2 <- as.integer(animalsmelt$Slice2)


## join animals and punches using the mouse id
alldata <- punches %>% 
  left_join(animals, by="Mouse")
str(alldata)

## 2. clean data 
alldata$Mouse <- as.factor(alldata$Mouse)
alldata$Date <- as.Date(alldata$Date, "%m/%d/%y")
alldata$Slice <- as.factor(alldata$Slice)
alldata$Punch <- as.factor(alldata$Punch)
alldata$Tube <- as.factor(alldata$Tube)
alldata$Slice.collector <- as.factor(alldata$Slice.collector)
alldata$Group <- as.factor(alldata$Group)
alldata$APA <- as.factor(alldata$APA)
alldata$Genotype <- as.factor(alldata$Genotype)
alldata$storagebox <- as.factor(alldata$storagebox)
alldata$Conflict <- as.factor(alldata$Conflict)
alldata$mouse.x <- NULL
str(alldata)


## FMR1 mice
FMR1_mice <- alldata %>%
  filter(Genotype %in% c("FMR1", "WT"))   %>%
  distinct(Mouse, Genotype, Group, Punch)

WT_Htt_H_Tr_Un <- alldata %>%
  filter(Genotype != "FMR1", 
         Genotype != "WT", 
         Punch %in% c("CA1", "CA2", "DG"))%>%
  distinct(Mouse, Genotype, Group, Punch) %>%
  arrange(Group)

WT_H_Tr_Un <- alldata %>%
  filter(Genotype != "FMR1", 
         Genotype != "WT", Genotype != "Htt (YAC128)", 
         Punch %in% c("CA1", "CA2", "DG"))%>%
  distinct(Mouse, Genotype, Group, Punch) %>%
  arrange(Group)


##
ARmice <- alldata %>%
  filter(Genotype %in% c("Htt (YAC128)", "WT (FVB/NJ)")) %>%
  filter(APA %in% c("HomeCage", "NoShock")) %>%
  filter(Punch %in% c("CA2", "DG")) %>%
  distinct(Mouse, Genotype, APA, Group)

## 3. Group Project: FMR1 - conflict - slices
FMR1_conflict_slices <- alldata %>%
  filter(Punch == "slice", Genotype %in% c("WT", "FMR1"),  Group %in% c("FMR1 Conflict Trained", "WT Conflict Trained"),
         Conflict == "Conflict", Slice != "5", 
         Mouse != "16-119C", Mouse != "16-119D", Mouse != "16-118C", Mouse != "16-118D") %>%
  select(Mouse, Group, Slice,Punch, Tube, storagebox, Atlas.location) %>%
  arrange(Group, Mouse) 
View(FMR1_conflict_slices) 
#write.csv(FMR1_conflict_slices, "FMR1_conflict_slices.csv", quote=F, row.names = FALSE)


FMR1_conflict_animals <- dcast(FMR1_conflict_slices, Mouse + Group ~ Slice, value.var="Slice")
View(FMR1_conflict_animals)


## 5. just the huntington comparisons
HttMice <- alldata %>%
  filter(Punch != "slice", Genotype != "FMR1", Genotype != "WT") %>%
  distinct(Mouse, Genotype, Group) %>%
  arrange(Mouse)
View(HttMice)

HttGroups <- alldata %>%
  filter(Punch != "slice", Genotype != "FMR1", Genotype != "WT")  %>%
  distinct(Mouse, Genotype, Group) %>%
  arrange(Group, Genotype) %>%
  summarise(count(Group))
str(HttGroups)

HttTissues <- alldata %>%
  filter(Punch != "slice", Genotype != "FMR1", Genotype != "WT")  %>%
  select(Mouse, Genotype, Group, Punch) %>%
  arrange(Group, Genotype, Punch) 
HttTissues <- (count(HttTissues, vars=c("Group","Punch")))
HttTissues <- spread(HttTissues, Punch, freq)


CA2 <- alldata %>%
  filter(Punch == "CA2") 
CA2counts <- (count(CA2, vars=c("Group","Mouse", "Punch")))