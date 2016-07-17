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
punches <- read.csv("NSBpunches.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
animals <- read.csv("NSBanimals.csv", header=TRUE, sep=",", stringsAsFactors = FALSE)
punches$mouse <- as.factor(punches$Mouse)
animals$mouse <- as.factor(animals$Mouse)


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


## 3. Group Project: FMR1 - conflict - slices
FMR1_conflict_slices <- alldata %>%
  filter(Punch == "slice", Genotype %in% c("WT", "FMR1"), 
         Conflict == "Conflict", Slice != "5", 
         Mouse != "16-119C", Mouse != "16-119D", Mouse != "16-118C", Mouse != "16-118D") %>%
  select(Mouse, Group, Slice,Punch, Tube, storagebox) %>%
  arrange(Group, Mouse) 
View(FMR1_conflict_slices) 

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


