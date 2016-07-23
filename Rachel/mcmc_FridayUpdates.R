##=============================================================================================================#
## Script created by Rayna Harris, rayna.harris@utexas.edu
## Script created in version R 3.3.1 
## This script is for analyzing data related to the STG qpcr projects
##=============================================================================================================#set path to data dir
setwd("~/Desktop/Rachel")

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

# 1. Load primer efficiency data -------------------

eff1  <- read.csv("2016-07-17_184859_ArctoFmr1.csv")
head(eff1)
summary(eff1) # For primer efficiencies, we care about Target.Name, CT, and Quantity

eff2  <- read.csv("2016-07-17_184859_FostoRpl19.csv")
head(eff2)
summary(eff2)

# 2. Clean data ----------------------
clean.eff1 <- rename(eff1, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 
clean.eff2 <- rename(eff2, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 

# 3. Create dilutions dataframe with quantity, target name, and ct.  ---------------

dilutions.eff1 <- clean.eff1 %>%
  filter(Task == "STANDARD") %>%
  select(dna, cq, gene) 
str(dilutions.eff1)

dilutions.eff2 <- clean.eff2 %>%
  filter(Task == "STANDARD") %>%
  select(dna, cq, gene) 
str(dilutions.eff2)

dilutions = rbind(dilutions.eff1,dilutions.eff2)
summary(dilutions)

# Fos1 needs to be deleted because it only had one measurement after quality control
dilutions = dilutions[!dilutions$gene=="Fos1",]

# 4. Calculate primer efficiences with MCMC.qpcr PrimEFF function ----------------------
quartz() # or "windows()" if you're using a PC
PrimEff(dilutions) # makes a plot with the primer efficiencies
amp.eff <- PrimEff(dilutions) #creates a table with the primer efficiencies
amp.eff

# Above is an example of how you would create a primer efficiency table that would be used downstream in MCMC.qpcr analysis. Those primers don't match the primers we used in the 7/21/2016 runs, so we'll make fake primer efficiencies for now (just for practice)

# Load sample data ----------------------- 

# Load Raina + Ariane data
data1 = read.csv("2016-07-21_183800_Raina_Ariane.csv") # CHANGE THIS TO YOUR FILE
names(data1)
data1 = data1[c(1,4,5,9)] # subset for Well (to decast), Sample.Name, Target.Name, and CT
head(data1) # looks good
# load Shayn + Nicole data
data2 = read.csv("2016-07-21_183800_ShaynCREB_NicolePrkcZ.csv")
data2 = data2[c(1,4,5,9)]# subset for Well (to decast), Sample.Name, Target.Name, and CT
head(data2) # looks good
# load Rachel data
data3 = read.csv("2016-07-21_213928_Rachel_Rpl19_Fmr1_Prkci.csv")
head(data3) 
# my run also included some standards, so lets remove those
data3 = data3[data3$Task=="UNKNOWN",]
data3 = data3[c(1,4,5,9)]# subset for Well (to decast), Sample.Name, Target.Name, and CT
head(data3)

# combine all three
data = rbind(data1,data2,data3)
head(data)
summary(data)

cleandata <- rename(data, c("Sample.Name"="sample", "Target.Name"="gene",  "CT"="cq", "Quantity"="dna")) 
head(cleandata)

# 5. create counts dataframe  ----------------------

counts <- cleandata %>%
  select(Well, sample, gene,  cq) %>%
  dcast(Well + sample ~ gene )

head(counts)

# 6. Join count and sample info, sort by sample, order in logical fashion ----------------
# read sample info, rename "Tube" to "sample"
samples <- read.csv("FMR1_conflict_slices.csv", header=TRUE)
head(samples)
colnames(samples)[1] <- "sample"
samples$Slice  <- as.factor(samples$Slice)
summary(samples)
head(counts)

data = merge(counts,samples,by="sample")
names(data)

# Make fake primer efficiencies so we can continue ----------- 
# NORMALLY YOU WOULD NOT HAVE TO DO THIS. JUST USE YOUR "amp.eff" OBJECT

amp.eff2 = as.data.frame(cbind("gene"=c("FMR1", "PrkcI", "CREB", "PrkcZ", "Rpl19"),"eff"=rep(2,5))) #INSERT YOUR GENE NAMES HERE
amp.eff2$eff = as.numeric(levels(amp.eff2$eff))
str(amp.eff2)
amp.eff2

# 7. Turn cq into counts -------------
names(data)
dd = cq2counts(
  data=data,
  genecols=c(3:7), # where the Cq data are in the data table
  condcols=c(1,9:10), # which columns contain factors
  effic=amp.eff2
)
head(dd)
tail(dd)

# Subset data UNDO THIS IF YOU WANT !!!! ---------------------------
dd = dd[!dd$Slice=="4",] # Remove slice 4
dd$Slice = factor(dd$Slice) # Remove slice 4 from the levels of the dd object
summary(dd)

# 8. Mixed model ------------------------
dd$Group = relevel(dd$Group,ref="WT Conflict Trained") # relevel to make WT the reference
mm=mcmc.qpcr(
  fixed="Group + Slice + Group:Slice",
  random="sample",
  data=dd
)
summary(mm)

# 9. Mixed model with diagnostic plots! -------------
mmd = mcmc.qpcr(
  fixed="Group + Slice + Group:Slice",
  random="sample",
  data=dd,
  pr=T,
  pl=T
)

quartz() # or "windows()" if you're using a PC
diagnostic.mcmc(
  model=mmd,
  col="grey50",
  cex=0.8
)

# 10. Plot the data!! ---------------------

# HPDplot(
#   model=mm,
#   factors="Slice",
#   main="Slice"
# )

S1 = HPDsummary(model=mm,data=dd)
HPDsummary(mm,dd,xgroup="Slice")
s0 = HPDsummary(model=mm,data=dd,relative=TRUE) # not super useful... it's just comparing the slices to slice1 WT

trellisByGene(S1,xFactor="Slice",groupFactor="Group")+xlab("Slice")

# Group only model--------
mmGroup = mcmc.qpcr(
  fixed="Group",
  random="sample",
  data=dd
)
summary(mmGroup)
hpdGroup = HPDsummary(model=mmGroup,data=dd)

# Slice only model -------
mmSlice = mcmc.qpcr(
  fixed="Slice",
  random="sample",
  data=dd,
  singular.ok=T
)
summary(mmSlice)

hpdSlice = HPDsummary(model=mmSlice,data=dd)

# Mixed model (with control gene) ---------
mmCon = mcmc.qpcr(
  fixed = "Group + Slice + Group:Slice",
  random = "sample",
  data = dd,
  controls="Rpl19", # you can put multiple control genes here (ex: controls=c("con1", "con2")),
  pr=T,
  pl=T
)
summary(mmCon)

# Diagnostics
quartz() # or "windows()" if you're using a PC
diagnostic.mcmc(
  model=mmCon,
  col="grey50",
  cex=0.8
)

# Couple ways to plot the data
hpdCon = HPDsummary(mmCon,dd)
trellisByGene(hpdCon,xFactor="Slice",groupFactor="Group")+xlab("Slice")
hpdCon = HPDsummary(mmCon,dd, xgroup = "Slice")

