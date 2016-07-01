## set working director to source file locations
setwd()
## first time install.packages("MCMC.qpcr")
library(MCMC.qpcr)
data(beckham.data)
data(beckham.eff)
qs=cq2counts(data=beckham.data, effic=beckham.eff, genecols=c(4:13),condcols=c(1:3), Cq1=37)
qs$treatment.time=as.factor(paste(qs$tr,qs$time,sep="."))
qs$treatment.time=relevel(qs$treatment.time,ref="control.0h")
naive=mcmc.qpcr(data=qs, fixed="treatment.time")
s1=HPDsummary(model=naive,data=qs)
s0=HPDsummary(model=naive,data=qs,relative=TRUE)
summary(naive)
##clear workspace

###now with the NS&B 2013 data comparing response to active place avoidance (APA) training in CA1 and CA3 hippocampus regions
ct=read.table("ct.ca1.ca3.txt",header=T,sep="\t") 
dil=read.table("dil.ca1.ca3.txt",header=T,sep="\t") 
eff=PrimEff(dil)
qs=cq2counts(data=ct, effic=eff, genecols=c(6:8),condcols=c(1:5), Cq1=37)
qs$Region.APA=as.factor(paste(qs$Region,qs$APA,sep="."))
qs$Region.APA=relevel(qs$Region.APA,ref="CA3.control")
naive=mcmc.qpcr(data=qs, fixed="Region.APA")
s1=HPDsummary(model=naive,data=qs)
s0=HPDsummary(model=naive,data=qs,relative=TRUE)
summary(naive)

