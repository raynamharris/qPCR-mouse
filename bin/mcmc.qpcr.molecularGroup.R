##=============================================================================================================#
## Script created by Rayna Harris, rayna.harris@utexas.edu
## Script created in version R 3.3.1 
## This script is for analyzing data related to the STG qpcr projects
##=============================================================================================================#

##script stored in "Z:/NSB_2016/IntegrativeNeuroscience/STGsingleneuron2015/bin"
## set path to data dir
setwd("Z:/NSB_2016/4_MouseMolecular/qPCR-mouse/molecularGroup")

## The process:
## 1. Loop over all experimental prjoect files and create one big "rawdata" dataframe
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

## 1. Loop over all experimental prjoect files and create one big "rawdata" dataframe
## uses read.xlsx function to read 1 sheet, sharting at row 8, only first 20 columns, nothing imported as a factor

file_list <- list.files(pattern = ".xls") #creates a string will all the .xls in a diretory for us to loop through

rm(rawdata) # first removed any dataframe called rawdata, if any
rm(temp_dataset) # first removed any dataframe called rawdata, if any

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("rawdata")){
    rawdata <- read.xlsx(file, sheetIndex = 1, startRow=8, colIndex = (1:20), stringsAsFactors=FALSE)
  }
  
  # if the merged dataset does exist, append to it
  if (exists("rawdata")){
    temp_dataset <-read.xlsx(file, sheetIndex = 1, startRow=8, colIndex = (1:20), stringsAsFactors=FALSE)
    rawdata<-rbind(rawdata, temp_dataset)
    rm(temp_dataset)
  }
  
}


names(rawdata)
str(rawdata)

## 2. Clean data

cleandata <- rawdata
cleandata <- rename(cleandata, c("Sample.Name"="sample", "Target.Name"="gene",  "C?."="cq", "Quantity"="dna")) 
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
samples <- read.csv("sample_info_2015.csv", header=TRUE, sep="," )

data <- inner_join(counts, samples) %>%
  arrange(sample) %>%
  select(sample, condition, CbNaV, IH, inx1, inx2)

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
