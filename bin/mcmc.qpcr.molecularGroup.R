##=============================================================================================================#
## Script created by Rayna Harris, rayna.harris@utexas.edu
## Script created in version R 3.3.1 
## This script is for analyzing data related to the STG qpcr projects
##=============================================================================================================#

##script stored in "Z:/NSB_2016/IntegrativeNeuroscience/STGsingleneuron2015/bin"
## set path to data dir
setwd("~/Desktop/Rachel/")
setwd("Z:/NSB_2016/4_MouseMolecular/qPCR-mouse/molecularGroup")


## The process:
## 1. Readin experimental prjoect files and create one big "rawdata" dataframe
## 2. Clean the data to make numbers numbers and rename important columns
## 3. Create dilutions dataframe with quantiry, target name, and ct.
## 4. Calculate primer efficiences with MCMC.qpcr PrimEFF function
## 5. create counts dataframe 
## 6. Join count and sample info, sort by sample, order in logical fashion
## 7. Turn cq into counts
## 8. Mixed model
## 9. Mixed model with diagnostic plots
## 10. Plot the data!!


#install.packages("dplyr")
#install.packages("plyr")
#install.packages("MCMC.qpcr")
#install.packages("reshape2")
library(dplyr)
library(plyr)
library(MCMC.qpcr)
library(reshape2)

# Load data -------------------

eff1  <- read.csv("2016-07-17_184859_ArctoFmr1.csv")
head(eff1)
summary(eff1) # For primer efficiencies, we care about Target.Name, CT, and Quantity

eff2  <- read.csv("2016-07-17_184859_FostoRpl19.csv")
head(eff2)
summary(eff1)

  
  ## Read data files
  ## for taqman, skip first 43 lines, 
  ## for sybr, skip first 44 lines
  
  RA <- read.delim("2016-07-20_112033_RA.txt", skip=43, header = TRUE, sep = "\t", stringsAsFactors = FALSE )[ ,1:24]
SN <- read.delim("2016-07-20_115248_ShaynNicole_CREB_PkcZ_Htt.txt", skip=44, header = TRUE, sep = "\t" , stringsAsFactors = FALSE)[ ,1:24]

## use rbind to bind data into a fill "raw data" dataframe
rawdata <- rbind(RA, SN)
names(rawdata)
str(rawdata)

## 2. Clean data
clean.eff1 <- rename(eff1, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 
clean.eff2 <- rename(eff2, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 

<<<<<<< HEAD:bin/mcmc.qpcr.molecularGroup.R
## 3. Create dilutions dataframe with quantity, target name, and ct. 
=======
  cleandata <- rawdata
cleandata <- rename(cleandata, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 
names(cleandata)
str(cleandata)
cleandata$dna <- as.numeric(cleandata$dna, na.rm = TRUE)
cleandata$cq <- as.numeric(cleandata$cq, na.rm = TRUE)
cleandata$gene <- as.factor(cleandata$gene)
str(cleandata)
>>>>>>> 09de634e82fe8d7ceaf84d44d5fde8000ab5efdf:bin/mcmc.qpcr.molecularGroup.R

dilutions.eff1 <- clean.eff1 %>%
  filter(Task == "STANDARD") %>%
  select(dna, cq, gene) 
str(dilutions.eff1)

dilutions.eff2 <- clean.eff2 %>%
  filter(Task == "STANDARD") %>%
  select(dna, cq, gene) 
str(dilutions.eff2)

dilutions = rbind(dilutions.eff1,dilutions.eff2)

## 4. Calculate primer efficiences with MCMC.qpcr PrimEFF function
quartz()
PrimEff(dilutions) # makes a plot with the primer efficiencies
amp.eff <- PrimEff(dilutions) #creates a table with the primer efficiencies

<<<<<<< HEAD:bin/mcmc.qpcr.molecularGroup.R
##### STOP HERE #####

=======
  >>>>>>> 09de634e82fe8d7ceaf84d44d5fde8000ab5efdf:bin/mcmc.qpcr.molecularGroup.R
## 5. create counts dataframe 
counts <- cleandata %>%
  filter(Task == "UNKNOWN") %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )

## 6. Join count and sample info, sort by sample, order in logical fashion
# read sample info, rename "Tube" to "sample", just read fir 16 fmr1 samples
samples <- read.csv("Z:/NSB_2016/4_MouseMolecular/qPCR-mouse/data/FMR1_conflict_slices.csv", header=TRUE, sep="," )
head(samples)
colnames(samples)[1] <- "sample"
samples <- samples[1:16,]

#use inner join to merge count sample dataframes
data <- inner_join(counts, samples) %>%
  arrange(sample) %>%
  select(sample, Mouse, Group, Slice, Atlas.location, Processor, 
         Ariane_Piota, NicolePkcZ, Raina_Fmr1, ShayneHtt, Raina_18S)

## 7. Turn cq into counts

dd=cq2counts(
  data=data,
  genecols=c(3:6), # where the Cq data are in the data table
  condcols=c(1:2), # which columns contain factors
  effic=amp.eff,
  Cq1=37
)
head(dd)

## 8. Mixed model
mm=mcmc.qpcr(
  fixed="condition",
  random="sample",
  data=dd
)
summary(mm)

## 9. Mixed model with diagnostic plots!
mmd=mcmc.qpcr(
  fixed="condition",
  random="sample",
  data=dd,
  pr=T,
  pl=T
)

diagnostic.mcmc(
  model=mmd,
  col="grey50",
  cex=0.8
)
dev.off()
## 10. Plot the data!!

HPDplot(
  model=mm,
  factors="conditionzymo",
  main="Maxwell vs Zymo"
)

S1=HPDsummary(model=mm,data=dd)
#png('Z:/NSB_2016/IntegrativeNeuroscience/qPCR-STG/ReproducibleExample/results/HPDsummary.png')
#plot(HPDsummary(model=mm,data=dd))
#dev.off()

#s0=HPDsummary(model=mm,data=dd,relative=TRUE)

## 2. Clean data

cleandata <- rawdata
cleandata <- rename(cleandata, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 
names(cleandata)
str(cleandata)
cleandata$dna <- as.numeric(cleandata$dna, na.rm = TRUE)
cleandata$cq <- as.numeric(cleandata$cq, na.rm = TRUE)
cleandata$gene <- as.factor(cleandata$gene)
str(cleandata)

## 3. Create dilutions dataframe with quantiry, target name, and ct. 

dilutions <- cleandata %>%
  filter(Task == "STANDARD") %>%
  select(dna, cq, gene) 
str(dilutions)


## 4. Calculate primer efficiences with MCMC.qpcr PrimEFF function
PrimEff(dilutions) # makes a plot with the primer efficiencies
amp.eff <- PrimEff(dilutions) #creates a table with the primer efficiencies

## 5. create counts dataframe 
counts <- cleandata %>%
  filter(Task == "UNKNOWN") %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )

## 6. Join count and sample info, sort by sample, order in logical fashion
# read sample info, rename "Tube" to "sample", just read fir 16 fmr1 samples
samples <- read.csv("Z:/NSB_2016/4_MouseMolecular/qPCR-mouse/data/FMR1_conflict_slices.csv", header=TRUE, sep="," )
head(samples)
colnames(samples)[1] <- "sample"
samples <- samples[1:16,]

#use inner join to merge count sample dataframes
data <- inner_join(counts, samples) %>%
  arrange(sample) %>%
  select(sample, Mouse, Group, Slice, Atlas.location, Processor, 
         Ariane_Piota, NicolePkcZ, Raina_Fmr1, ShayneHtt, Raina_18S)

## 7. Turn cq into counts

dd=cq2counts(
  data=data,
  genecols=c(3:6), # where the Cq data are in the data table
  condcols=c(1:2), # which columns contain factors
  effic=amp.eff,
  Cq1=37
)
head(dd)

## 8. Mixed model
mm=mcmc.qpcr(
  fixed="condition",
  random="sample",
  data=dd
)
summary(mm)

## 9. Mixed model with diagnostic plots!
mmd=mcmc.qpcr(
  fixed="condition",
  random="sample",
  data=dd,
  pr=T,
  pl=T
)

diagnostic.mcmc(
  model=mmd,
  col="grey50",
  cex=0.8
)
dev.off()
## 10. Plot the data!!

HPDplot(
  model=mm,
  factors="conditionzymo",
  main="Maxwell vs Zymo"
)

S1=HPDsummary(model=mm,data=dd)
#png('Z:/NSB_2016/IntegrativeNeuroscience/qPCR-STG/ReproducibleExample/results/HPDsummary.png')
#plot(HPDsummary(model=mm,data=dd))
#dev.off()

#s0=HPDsummary(model=mm,data=dd,relative=TRUE)
