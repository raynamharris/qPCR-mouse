setwd("Z:/rmharris/Course Developers/Hofmann Fenton Collaboration")
##install MCMCglmm and MCMC.qpcr
##load libraries

#install.packages("MCMCglmm")
#install.packages("MCMC.qpcr")
library("MCMCglmm")
library("MCMC.qpcr")


##read sample cts
ct=read.csv("samples with all gens only.csv") 
head(ct)

##read standard curves, calculate efficiencies
dil=read.table("effeciencies_ca3.txt",header=T,sep="\t") 
head(dil)
eff=PrimEff(dil)

##turn ct into molecule counts
genecolumns=c(7:18) 
conditions=c(1:6)
counts=cq2counts(data=ct, genecols=c(7:18),condcols=c(1:6),effic=eff)
head(counts)
write.table(counts, "counts for samples with all genes only.csv", col.names = NA, sep = ",")

counts=read.csv("counts for samples with all genes only. no home.v2.csv")
counts$strain.APA=factor(paste(counts$strain, counts$APA))


head(counts)

###linear mixed model 1 
model1=mcmc.qpcr(data=counts, fixed="strain+APA+APA:strain",random=c("sample","time", "year"),pr=TRUE,pl=TRUE)
summary(model1)
diagnostic.mcmc(model1)
plota=HPDsummary(model=model1,data=counts)
poltb=HPDsummary(model=model1,data=counts,relative=TRUE)

model2=mcmc.qpcr(data=counts, fixed="APA+strain+strain:APA",random=c("sample","time", "year"),pr=TRUE,pl=TRUE)
summary(model2)
diagnostic.mcmc(model2)
plota=HPDsummary(model=model2,data=counts)
poltb=HPDsummary(model=model2,data=counts,relative=TRUE)

model3=mcmc.qpcr(data=counts, fixed="strain.APA",random=c("sample","time"),pr=TRUE,pl=TRUE)
summary(model3)
diagnostic.mcmc(model3)
plota=HPDsummary(model=model3,data=counts)
poltb=HPDsummary(model=model3,data=counts,relative=TRUE)
summary(plota)


pd = position_dodge(0.3)
ggplot(
  plota$summary, 
  aes_string(x = "strain.APA", colour = "strain.APA",
             y = "mean")) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), lwd = 0.4, width = 0.7, position = pd) + geom_point(position = pd, size = 2.5) + 
  theme_bw() +facet_wrap(~gene,scales="free_y",nrow=2) +ylab("log(abundance)") +
  theme(legend.position="bottom") +
  theme(strip.text.x=element_text(face="italic",size = 12, family="serif")) +
  # edit the next line according to the names of your factors:
  xlab("Active Place Avoidance Training") + 
  # remark this line to see if defalt colors are good for you, if not, edit as needed
  scale_color_manual(values=c("#000000", "#FF0000", "#999999", "#FF9933")
                     + theme(axis.title.y = element_text(size=20))
  )

#install.packages("reshape2")
library(reshape2)
counts.slim=counts[,c(1,2,4,6,10)]


counts.l <- aggregate(count ~ ind + strain.APA+region +gene, data = counts.slim, mean)
counts.l <- aggregate(count ~ ind + strain.APA+region +gene, data = counts.slim)



attach(data)

#big group differences
par(mfrow=c(3,3))
boxplot(cam2kd~strain,data=data, main="cam2kd", ylab="counts")
boxplot(dlg4~strain,data=data, main="dlg4", ylab="counts")
boxplot(gria~strain,data=data, main="gria", ylab="counts")
boxplot(grim~strain,data=data, main="grim", ylab="counts")
boxplot(grin~strain,data=data, main="grin", ylab="counts")
boxplot(nsf~strain,data=data, main="nsf", ylab="counts")
boxplot(pkmz~strain,data=data, main="pkmz", ylab="counts")
boxplot(rpl19~strain,data=data, main="rpl19", ylab="counts")
boxplot(rRNA18S~strain,data=data, main="rRNA18S", ylab="counts")

#big group differences
boxplot(cam2kd~APA,data=data, main="cam2kd", ylab="counts")
boxplot(dlg4~APA,data=data, main="dlg4", ylab="counts")
boxplot(gria~APA,data=data, main="gria", ylab="counts")
boxplot(grim~APA,data=data, main="grim", ylab="counts")
boxplot(grin~APA,data=data, main="grin", ylab="counts")
boxplot(nsf~APA,data=data, main="nsf", ylab="counts")
boxplot(pkmz~APA,data=data, main="pkmz", ylab="counts")
boxplot(rpl19~APA,data=data, main="rpl19", ylab="counts")
boxplot(rRNA18S~APA,data=data, main="rRNA18S", ylab="counts")

#effects of training and strain
par(mfrow=c(2,2))
boxplot(cam2kd~strain*APA, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="cam2kd", xlab="Strain and Training")
boxplot(creb~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="creb", xlab="Strain and Training")
boxplot(dlg4~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="dlg4", xlab="Strain and Training")
boxplot(fmr1~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="fmr1", xlab="Strain and Training")
boxplot(fos~strain*APA, data=data, notch=FALSE,   col=(c("cyan","yellow", "deeppink")),
        main="fos", xlab="Strain and Training")
boxplot(gria~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="gria", xlab="Strain and Training")
boxplot(grim~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="grim", xlab="Strain and Training")
boxplot(grin~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="grin", xlab="Strain and Training")
boxplot(nsf~strain*APA, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="nsf", xlab="Strain and Training")
boxplot(pkmz~strain*APA, data=data, notch=FALSE, col=(c("cyan","yellow", "deeppink")),
        main="pkmz", xlab="Strain and Training")
boxplot(rpl19~strain*APA, data=data, notch=FALSE, col=(c("cyan","yellow", "deeppink")),
        main="rpl19", xlab="Strain and Training")
boxplot(rRNA18S~strain*APA, data=data, notch=FALSE, col=(c("cyan","yellow", "deeppink")),
        main="rRNA18S", xlab="Strain and Training")

#effects of training, starin, and region
par(mfrow=c(3,1))
boxplot(cam2kd~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="cam2kd", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(creb~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="creb", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(dlg4~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="dlg4", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(fmr1~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="fmr1", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(fos~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="fos", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(gria~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="gria", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(grim~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="grim", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(grin~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="grin", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(nsf~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="nsf", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(pkmz~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="pkmz", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(rpl19~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="rpl19", xlab="Strain * Training * Region", cex.axis=0.7)
boxplot(rRNA18S~APA*strained*region, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="rRNA18S", xlab="Strain * Training * Region", cex.axis=0.7)

#correlations
pairs(~cam2kd+creb+dlg4+fmr1+fos+gria+grim+grin+nsf+pkmz+rpl19+rRNA18S,data=data, 
      main="Simple Scatterplot Matrix")
pairs(~cam2kd+creb+dlg4+fmr1,data=data,main="Simple Scatterplot Matrix")
pairs(~fos+gria+grim+grin,data=data,main="Simple Scatterplot Matrix")
pairs(~nsf+pkmz+rpl19+rRNA18S,data=data, main="Simple Scatterplot Matrix")

install.packages("corrplot")
library(corrplot)
forcor=(data)
forcor$ind=NULL
forcor$time=NULL
forcor$region=NULL
forcor$APA=NULL
forcor$strained=NULL
forcor$sTrain=NULL
M=cor(forcor)
par(mfrow=c(1,1))
corrplot(M, method = "ellipse")
corrplot.mixed(M, lower = "ellipse", upper = "circle")

detach(data)


###linear mixed model 2
model2=mcmc.qpcr(data=counts, fixed="region+APA+APA:region",random=c("sample","time"),pr=TRUE,pl=TRUE)
summary(model2)
diagnostic.mcmc(model2)
plotc=HPDsummary(model=model2,data=counts)
plotd=HPDsummary(model=model2,data=counts,relative=TRUE)
print2=getNormalizedData(model2,counts)
datamodel2<-data.frame(cbind(print2$normData,print2$conditions))

###linear mixed model 3
model3=mcmc.qpcr(data=counts, fixed="region+APA+APA:region",random=c("sample","time"), controls=c("rpl19"),pr=TRUE,pl=TRUE)
summary(model3)
diagnostic.mcmc(model3)
plotc=HPDsummary(model=model3,data=counts)
plotd=HPDsummary(model=model3,data=counts,relative=TRUE)
print3=getNormalizedData(model3,counts)
datamodel3<-data.frame(cbind(print3$normData,print3$conditions))