##install MCMCglmm and MCMC.qpcr
library("MCMCglmm", lib.loc="C:/Program Files/R/R-3.0.1/library")
library("MCMC.qpcr", lib.loc="C:/Program Files/R/R-3.0.1/library")

setwd("Z:/Shared/NSB/Research/MouseMolBioData/Data Analysis")

##########Determine counts for Ca1 samples
eff1=read.table("effeciencies_ca1.txt",header=T,sep="\t") 
head(eff1)
PrimEff(eff1)
#add effeciency values to a table caleed cq1.txt then read the file with this command
amp.eff1=read.table("cq1_ca1.txt",header=T,sep="\t")
#########read in Ca1 data
ct1=read.table("ct_ca1.txt",header=T,sep="\t") 
head(ct1)
genec=c(5:7) 
conds=c(1:4)
#########turn cq into counts
dd1=cq2counts(data=ct1, genecols=c(5:7),condcols=c(1:4), effic=amp.eff1)
head(dd1)
write.table(dd1, "ca1_cq2counts.csv", col.names = NA, sep = ",")

##########Determine counts for Ca3 samples
eff3=read.table("effeciencies_ca3_conc.txt",header=T,sep="\t") 
head(eff3)
PrimEff(eff3)
#add effeciency values to a table caleed cq1_ca3.txt then read the file with this command
amp.eff3=read.table("cq1_ca3_conc.txt",header=T,sep="\t")
#########read in Ca3 data
ct3=read.table("ct_ca3_conc.txt",header=T,sep="\t") 
head(ct3)
genec=c(5:7) 
conds=c(1:4)
#############turn cq into counts
dd3=cq2counts(data=ct3, genecols=c(5:8),condcols=c(1:4), effic=amp.eff3)
head(dd3)
write.table(dd3, "ca3_cq2counts_conc.csv", col.names = NA, sep = ",")

##############open oup both cq2 counts documents and make a combo file
ddcombo.slim=read.table("ca1ca3_counts_slim.txt", header=T,sep="\t")

###linear mixed models 
mm.slim=mcmc.qpcr(data=ddcombo.slim, fixed="Region+APA+APA:Region",random=c("sample"))
summary(mm.slim)

CA1=list(control=list(factors=0), trained=list(factors=c(0,"APAtrain")))
CA3=list(control=list(factors=c(0,"RegionCA3")), trained=list(factors=c(0,"RegionCA3","APAtrain")))

plot.hpd.bygene.bygroup(model=mm.slim, gene="pkmz", group1=CA1, group2=CA3, ,jitter=.25, group3 = NULL,)
legend(0.45,5,"CA1",lty=1,pch=1,col="coral", bty="n")
legend(0.45,4.5,"CA3",lty=1,pch=1,col="cyan", bty="n")
legend(0.65,3,"a",bty="n")
legend(.9,4.4,"b",bty="n")
legend(1.65,3.4,"a",bty="n")
legend(1.9,5.2,"b",bty="n")

plot.hpd.bygene.bygroup(model=mm.slim, gene="rpl19", group1=CA1, group2=CA3, jitter=.25, group3 = NULL)
legend(2,8.5,"CA1",lty=1,pch=1,col="coral", bty="n")
legend(2,8,"CA3",lty=1,pch=1,col="cyan", bty="n")
legend(0.65,8,"a",bty="n")
legend(0.9,4.6,"b",bty="n")
legend(1.65,8.8,"a",bty="n")
legend(1.9,5.7,"b",bty="n")

plot.hpd.bygene.bygroup(model=mm.slim, gene="grim", group1=CA1, group2=CA3, group3 = NULL)
legend(0.45,4,"CA1",lty=1,pch=1,col="coral", bty="n")
legend(0.45,3.5,"CA3",lty=1,pch=1,col="cyan", bty="n")
legend(0.65,2.5,"a",bty="n")
legend(.9,2.8,"a",bty="n")
legend(1.7,3.5,"a",bty="n")
legend(1.9,4,"a",bty="n")

effects1=plot.hpd(model=mm.slim, factors="RegionCA3", plot=FALSE)
effects1

effects2=plot.hpd(model=mm.slim, factors="APAtrain", plot=FALSE)
effects2

effects3=plot.hpd(model=mm.slim, factors="APAtrain:RegionCA3", plot=FALSE)
effects3

all.effects=rbind(effects1,effects2,effects3)
padj.qpcr(all.effects)

#####plots
plot.hpd(model=mm.slim, factors="APAtrain", main="APAtrain", jitter=-0.15, ylim=c(-5,5))
plot.hpd(model=mm.slim, factors="RegionCA3", main="Region", jitter=-0.15, ylim=c(-5,5))
plot.hpd(model=mm.slim, factors="RegionCA3:APAtrain", main="APAtrain:RegionCA3", jitter=-0.15, ylim=c(-5,5))

APA=list(ddcomboslim=list(factors=c(6,"train")), control=list(factors=c(6,"control")))
plot.hpd.bygene(model=informed,gene="pkmz",conditions=APA,main="pkmz")

###infomred plots
informed=mcmc.qpcr(fixed="Time+APA+APA:Time", random="sample",data=ddcombo.slim,controls=c("rpl19"),m.fix=1.2)
summary(informed)
plot.hpd(model=informed, factors="APAtrain", main="APAtrain")
plot.hpd(model=informed, factors="TimePM", main="TimePM")
plot.hpd(model=informed, factors="TimePM:APAtrain", main="TimePM:APAtrain")

informed=mcmc.qpcr(fixed="Region+APA+APA:Region", random=c("sample"),data=ddcombo.slim,controls=c("rpl19"),m.fix=1.2)
summary(informed)
plot.hpd(model=informed, factors="APAtrain", main="Effect of training")
plot.hpd(model=informed, factors="RegionCA3", main="Effect of brain region, CA3 vs CA1")
plot.hpd(model=informed, factors="RegionCA3:APAtrain", main="Interactions")

plot.hpd(model=informed,factors=c("APAtrain"),main="training / region",jitter=-0.15,ylim=c(-5,5))
points.hpd(model=informed,factors=c("APAtrain","RegionCA3"),jitter=0.15,col="orange2",pch=17)
points.hpd(model=informed,factors=c("APAtrain:RegionCA3"),jitter=0.-15,col="grey50",pch=17)

#### nieve individual plots
plot.hpd(model=mm.slim, factors="RegionCA3", main="nieve RegionCA3")
plot.hpd(model=mm.slim, factors="APAtrain", main="APAtrain")
plot.hpd(model=mm.slim, factors="RegionCA3:APAtrain", main="Interactions")
plot.hpd(model=mm.slim, factors="RegionCA3","APAtrain", main="train and region")

plot.hpd(model=mm.slim,factors=c("APAtrain"),main="training / region",jitter=-0.15,ylim=c(-5,5))
points.hpd(model=mm.slim,factors=c("APAtrain","RegionCA3"),jitter=0.15,col="orange2",pch=17)

######determining Cq1

setwd("Z:/Shared/NSB/Research/MouseMolBioData/Data Analysis")
dilutions=read.csv("effeciencies_ca3_rpl19.csv", row.names=NULL)
head(dilutions)
rpl=dilutions

#primer efficiency
rpl=rpl[!is.na(rpl$cq),]
rpl=rpl[rpl$cq!=-1,]
plot(cq~log(dna,2),rpl)
ll=lm(cq~log(dna,2),rpl)
summary(ll)
coef(ll)
slope=coef(ll)[2]
eff=2^(-1/slope)
eff #1.88919
abline(ll)

#calculating cq1
max=max(rpl$cq)*log(eff,2)
c2=rpl$cq[rpl$cq>max]
plot(c2[order(c2)])
cq1=mean(c2,na.rm=T)
cq1 #37.4575
cq1*log(eff,2) #34.37729

# Boxplots
setwd("Z:/Shared/NSB/Research/MouseMolBioData/Data Analysis")
mouse=read.table("for box plot2.txt", header=T,sep="\t")

op <- par(mar = c(6,6,2,2)+ 0)
boxplot(pkmz_CA3~Group2,data=mouse, ylab="log expression", cex=1, cex.axis=1.8, cex.lab=2.5)

boxplot(grim.1_CA3~Group2,data=mouse,ylab="log expression", cex=2, cex.axis=1.8, cex.lab=2.5)
boxplot(dlg4_CA3~Group2,data=mouse, xlab="region.treatment", ylab="log expression", cex=2, cex.axis=2, cex.lab=1.8)

