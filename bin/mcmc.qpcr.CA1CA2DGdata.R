##=============================================================================================================#
## Script created by Rayna Harris, rayna.harris@utexas.edu
## Script created in version R 3.3.1 
## This script is for analyzing data related to the STG mouse projects
##=============================================================================================================#set path to data dir

#install.packages("dplyr")
#install.packages("plyr")
#install.packages("MCMC.qpcr")
#install.packages("reshape2")
#install.packages("car")
install.packages("wesanderson")
source("https://bioconductor.org/biocLite.R")
biocLite("ggPlot")
pallibrary(dplyr)
library(plyr)
library(MCMC.qpcr)
library(reshape2)
library(car)
library(wesanderson)
library(ggplot2)

names(wes_palettes)
library(ggplot2)
ggplot(iris, aes(Sepal.Length, Sepal.Width, color = Species)) + 
  geom_point(size = 3) + 
  scale_color_manual(values = wes_palette("GrandBudapest")) + 
  theme_gray()
wes_palette("Moonrise2")

setwd("/Volumes/nsb/NSB_2016/4_MouseMolecular/Final Project")
#setwd("Z:/NSB_2016/4_MouseMolecular/Final Project/CA1CA3DGsamplesdata")

# Run primer efficiency script first  ----------------------- 

## Go to Z:/NSB_2016/4_MouseMolecular/Final Project/standardcurvedata
## run mcmc.qpcr.standardcurvedata.R


# read sample info ----------------

samples <- read.csv("qPCR_95samples_072916.csv", header=TRUE)
head(samples)
colnames(samples)[1] <- "sample"
samples$Slice  <- as.factor(samples$Slice)
samples$Genotype  <- as.factor(samples$Genotype)
summary(samples)


# Load sample data with a for loop ----------------------- 
setwd("/Volumes/nsb/NSB_2016/4_MouseMolecular/Final Project/CA1CA3DGsamplesdata")
file_list <- list.files(pattern = ".csv") #creates a string will all the .xls in a diretory for us to loop through

rm(rawdata) # first removed any dataframe called rawdata, if any
rm(temp_dataset) # first removed any dataframe called rawdata, if any

for (file in file_list){
  # if the merged dataset doesn't exist, create it
  if (!exists("rawdata")){
    rawdata <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)[ ,1:24]
  }
  # if the merged dataset does exist, append to it
  else {
    temp_dataset <- read.csv(file, header = TRUE, stringsAsFactors = FALSE)[ ,1:24]
    rawdata<-rbind(rawdata, temp_dataset)
    rm(temp_dataset)
  }
}
str(rawdata)
names(rawdata)

rawdata$Target.Name  <- as.factor(rawdata$Target.Name) #make a factor
str(rawdata)

rawdata <- rawdata[c(1,4,5,9)] # subset for Well (to decast), Sample.Name, Target.Name, and CT
rawdata <- rename(rawdata, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq")) ## rename to fit with MCMC.qpcr
head(rawdata) # looks good

## rename a few misspelled genes
rawdata$gene <- recode(rawdata$gene, "'V1ar' = 'V1aR'")
rawdata$gene <- recode(rawdata$gene, "'Otxr' = 'OxtR'")
rawdata$gene <- recode(rawdata$gene, "'Nr3c1' = 'GR'")
rawdata$gene <- recode(rawdata$gene, "'V1br' = 'V1bR'")
rawdata$gene <- recode(rawdata$gene, "'RPL_19' = 'Rpl19'")


rawdata$gene <- recode(rawdata$gene, "'BDNF' = 'Bdnf'")
rawdata$gene <- recode(rawdata$gene, "'CamK2' = 'Camk2a'")

str(rawdata)
levels(rawdata$gene)

## Combine sample datasets as appropriate ----------------------- 

#### Shayn ----------
Shayn <- rawdata %>%
  filter(gene %in% c("PkmZ", "PkmIL"))
Shayn <- Shayn %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )  ## widen the dataframe
Shayn <- merge(Shayn,samples,by="sample")  ## join the cq and sample info
names(Shayn)
Shayn <- select(Shayn, sample, Well, Slice, Mouse, Category, 
                 Punch, Genotype, APA, Allen, Extractor, Pipettor,
                 PkmIL, PkmZ)  
head(Shayn)

Shayndd <- cq2counts(
  data=Shayn,
  genecols=c(12:13), # where the Cq data are in the data table
  condcols=c(1:11), # which columns contain factors
  effic=amp.eff
)
head(Shayndd)

Shayn$Genotype <- relevel(Shayn$Genotype,ref="WT") # relevel to make WT the reference
Shaynmm <- mcmc.qpcr(
  fixed="Genotype + Punch + Genotype:Punch",
  random="sample",
  data=Shayndd
)
summary(Shaynmm)

S1 <- HPDsummary(model=Shaynmm,data=Shayndd)
HPDsummary(Shaynmm,Shayndd,xgroup="Genotype")
s0 = HPDsummary(model=Shaynmm,data=Shayndd,relative=TRUE) # not super useful... it's just comparing the slices to slice1 WT

trellisByGene(S1,xFactor="Genotype",groupFactor="Punch")+xlab("Genotype")
trellisByGene(S1,xFactor="Genotype",groupFactor="Punch")+xlab("Genotype")





#### social ---------------

social <- rawdata %>%
  filter(gene %in% c("GR","OxtR", "CrH", "Rpl19"))
social <- social %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )  ## widen the dataframe
social <- merge(social,samples,by="sample")  ## join the cq and sample info
names(social)
str(social)

## select and reorder the columsn
social <- social %>%
  select(sample, Well, Slice, Mouse, Category, 
                Punch, Genotype, APA, Allen, Extractor, Pipettor,
                CrH, GR, OxtR, Rpl19)  %>%   
  filter(Punch != "CA1") %>% droplevels() %>%   
  filter(APA != "Trained") %>% droplevels() %>%   
  filter(Genotype != "Htt") %>% droplevels()
head(social)
#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
social$Punch <- relevel(social$Punch,ref="DG") # relevel to make WT the reference
social$Category <- relevel(social$Category,ref="NoStress") # relevel to make WT the reference

levels(social$Punch)

socialdd <- cq2counts(
  data=social,
  genecols=c(12:15), # where the Cq data are in the data table
  condcols=c(1:11), # which columns contain factors
  effic=amp.eff
)
head(socialdd)

#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
socialmm <- mcmc.qpcr(
  fixed= "Category + Punch + Punch:Category",
  random="sample", 
  data=socialdd,
  singular.ok=TRUE,
  pr=TRUE
)
summary(socialmm)

print1=getNormalizedData(socialmm,socialdd)
data<-data.frame(cbind(print1$normData,print1$conditions))

S1 <- HPDsummary(model=socialmm,data=socialdd)

HPDsummary(model=socialmm,data=socialdd,xgroup="Category")
s0 = HPDsummary(model=socialmm,data=socialdd,relative=TRUE) # not super useful... it's just comparing the slices to slice1 WT

Nagarajeffect <- trellisByGene(S1,xFactor="Punch",groupFactor="Category")+xlab("Genotype")
HPDsummary






#### social 2 ---------------

social <- rawdata %>%
  filter(gene %in% c("GR", "V1aR", "V1bR", "OxtR", "CrH", "Rpl19", "PkmZ", "PkmIL"))
social <- social %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )  ## widen the dataframe
social <- merge(social,samples,by="sample")  ## join the cq and sample info
names(social)

## select and reorder the columsn
social <- social %>%
  select(sample, Well, Slice, Mouse, Category, 
         Punch, Genotype, APA, Allen, Extractor, Pipettor,
         CrH, GR, OxtR, Rpl19)  %>%   
  filter(Punch != "CA1") %>% droplevels() %>%   
#  filter(APA != "Trained") %>% droplevels() %>%   
  filter(Genotype != "Htt") %>% droplevels()
head(social)
#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
social$Punch <- relevel(social$Punch,ref="DG") # relevel to make WT the reference
social$Category <- relevel(social$Category,ref="NoStress") # relevel to make WT the reference

levels(social$Punch)

socialdd <- cq2counts(
  data=social,
  genecols=c(12:17), # where the Cq data are in the data table
  condcols=c(1:11), # which columns contain factors
  effic=amp.eff
)
head(socialdd)
levels(socialdd$gene)

#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
socialmm <- mcmc.qpcr(
  fixed= "Category + Punch + Punch:Category",
  random="sample", 
  #control="Rpl19",
  data=socialdd,
  singular.ok=TRUE 
  #normalize=TRUE
)
summary(socialmm)


pal <- wes_palette("Moonrise2")
palette <- scale_color_manual(values=pal)
S1 <- HPDsummary(model=socialmm,data=socialdd)
S1$ggPlot + palette


p <- trellisByGene(S1,xFactor="Category",groupFactor="Punch")+xlab("Category")
p$ggPlot + palette


HPDsummary(model=socialmm,data=socialdd,xgroup="Category")

pal <- wes_palette("Moonrise2")
palette <- scale_color_manual(values=pal)
s0 <- HPDsummary(model=socialmm,data=socialdd,relative=TRUE) # not super useful... it's just comparing the slices to slice1 WT
S0$ggPlot + palette


Nagarajeffect <- trellisByGene(S1,xFactor="Punch",groupFactor="Category")+xlab("Genotype")




#### social 3 ---------------

social <- rawdata %>%
  filter(gene %in% c("GR", "V1aR", "V1bR", "OxtR", "CrH", "Rpl19"))
social <- social %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )  ## widen the dataframe
social <- merge(social,samples,by="sample")  ## join the cq and sample info
names(social)

## select and reorder the columsn
social <- social %>%
  select(sample, Well, Slice, Mouse, Category, 
         Punch, Genotype, APA, Allen, Extractor, Pipettor,
         CrH, GR, OxtR, Rpl19, V1bR, V1aR)  %>%   
  filter(Punch == "CA2") %>% droplevels()  
  #  filter(APA != "Trained") %>% droplevels() %>%   
  # filter(Genotype != "Htt") %>% droplevels()
head(social)
#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
#social$Punch <- relevel(social$Punch,ref="DG") # relevel to make WT the reference
social$Category <- relevel(social$Category,ref="NoStress") # relevel to make WT the reference

levels(social$Punch)

socialdd <- cq2counts(
  data=social,
  genecols=c(12:16), # where the Cq data are in the data table
  condcols=c(1:11), # which columns contain factors
  effic=amp.eff
)
head(socialdd)

#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
socialmm <- mcmc.qpcr(
  fixed= "Category + Genotype + Genotype:Category",
  random="sample", 
  control="Rpl19",
  data=socialdd,
  singular.ok=TRUE,
  normalize=TRUE
)
summary(socialmm)

S1 <- HPDsummary(model=socialmm,data=socialdd)
HPDsummary(model=socialmm,data=socialdd,xgroup="Genotype")
s0 = HPDsummary(model=socialmm,data=socialdd,relative=TRUE) # not super useful... it's just comparing the slices to slice1 WT

Nagarajeffect <- trellisByGene(S1,xFactor="Punch",groupFactor="Category")+xlab("Genotype")
??HPDsummary








#### social 4 ---------------

social <- rawdata %>%
  filter(gene %in% c("GR", "V1aR", "V1bR", "OxtR", "CrH"))
social <- social %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )  ## widen the dataframe
social <- merge(social,samples,by="sample")  ## join the cq and sample info
names(social)

## select and reorder the columsn
social <- social %>%
  select(sample, Well, Slice, Mouse, Category, 
         Punch, Genotype, APA, Allen, Extractor, Pipettor,
         CrH, GR, OxtR, Rpl19, V1bR, V1aR)  %>%   
  filter(Punch == "CA2") %>% droplevels()  
#  filter(APA != "Trained") %>% droplevels() %>%   
# filter(Genotype != "Htt") %>% droplevels()
head(social)
social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
#social$Punch <- relevel(social$Punch,ref="DG") # relevel to make WT the reference

levels(social$Genotype)

socialdd <- cq2counts(
  data=social,
  genecols=c(12:15), # where the Cq data are in the data table
  condcols=c(1:11), # which columns contain factors
  effic=amp.eff
)
head(socialdd)

#social$Genotype <- relevel(social$Genotype,ref="WT") # relevel to make WT the reference
socialmm <- mcmc.qpcr(
  fixed= "Category + Genotype + Genotype:Category",
  random="sample", 
  data=socialdd,
  singular.ok=TRUE
)
summary(socialmm)

S1 <- HPDsummary(model=socialmm,data=socialdd)
HPDsummary(model=socialmm,data=socialdd,xgroup="Genotype")
s0 = HPDsummary(model=socialmm,data=socialdd,relative=TRUE) # not super useful... it's just comparing the slices to slice1 WT

Nagarajeffect <- trellisByGene(S1,xFactor="Punch",groupFactor="Category")+xlab("Genotype")
??HPDsummary
