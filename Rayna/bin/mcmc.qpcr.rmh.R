##=============================================================================================================#
## Script created by Rayna Harris, rayna.harris@utexas.edu
## Script created in version R 3.3.1 
## This script is for analyzing data related to the mouse qpcr projects
##=============================================================================================================#

##script stored in "Z:/NSB_2016/IntegrativeNeuroscience/qPCR-mouse/Rayna/bin"
## set path to data dir
setwd("Z:/NSB_2016/IntegrativeNeuroscience/qPCR-mouse/Rayna/data")

## The process:
## 1. Loop over all experimental prjoect files and create one big "rawdata" dataframe
## 2. Clean the data to make numbers numbers and rename important columns
## 3. Wrangle standard curve data with dplyr and plyr package into a dataframe called dilutions
## 4. Calculate primer efficiences with MCMC.qpcr PrimEFF function
## 5. Wrangle sample info a with dplyr packag
## 6. Analyze data with MCMC.qpcr package

#install.packages("xlsx", dependencies = TRUE) #note, must have Java installed on computer
#install.packages("plyr")
#install.packages("dplyr")
#install.packages("MCMC.qpcr")
#install.packages("reshape2")
library(xlsx)
library(dplyr)
library(plyr)
library(MCMC.qpcr)
library(reshape2)

## 1. read raw data starting with row 42, then make quanitity a real number

file_list <- list.files() #creates a string will all the files in a diretory for us to loop through

rm(rawdata) # first removed any dataframe called rawdata

for (file in file_list){
  
  # if the merged dataset doesn't exist, create it
  if (!exists("rawdata")){
    rawdata <- read.xlsx("2016-06-29_183933-mouse-RnRNA-tenfold-rmh.xls", sheetIndex = 1, startRow=42, stringsAsFactors=FALSE)
    
  }
  
  # if the merged dataset does exist, append to it
  if (exists("rawdata")){
    rawdata <- read.xlsx("2016-06-29_183933-mouse-RnRNA-tenfold-rmh.xls", sheetIndex = 1, startRow=42, stringsAsFactors=FALSE)
    rawdata<-rbind(rawdata, temp_dataset)
    rm(temp_dataset)
  }
  
}




rawdata <- read.xlsx("2016-06-29_183933-mouse-RnRNA-tenfold-rmh.xls", sheetIndex = 1, startRow=42, stringsAsFactors=FALSE)
names(rawdata)
str(rawdata)
rawdata$Quantity <- as.numeric(rawdata$Quantity)
str(rawdata$Quantity)



## 2. Create dilutions dataframe with quantiry, target name, and ct. 
##    Then rename for MCMC.qpcr
dilutions <- rawdata %>%
  filter(Task == "STANDARD") %>%
  select(Quantity, CT, Target.Name)
dilutions <- rename(dilutions, c("Quantity"="dna", "CT"="cq", "Target.Name"="gene")) 
str(dilutions)



## 3. Calculate primer efficiences with MCMC.qpcr PrimEFF function
PrimEff(dilutions) # makes a plot with the primer efficiencies
eff <- PrimEff(dilutions) #creates a table with the primer efficiencies


## 4. create counts dataframe 
counts <- rawdata %>%
  filter(Task == "UNKNOWN") %>%
  select(Sample.Name,Target.Name,CT)
str(counts)
counts <- rename(counts, c("Sample.Name"="sample", "Target.Name"="gene", "CT"="cq")) 


