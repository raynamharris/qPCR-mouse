##install MCMCglmm and MCMC.qpcr
##load libraries
library("MCMCglmm")
library("MCMC.qpcr")
##set working directory

##read sample cts
ct=read.csv("samples with all gens only.csv") 
head(ct)

##read standard curves, calculate efficiencies
dil=read.table("dilutions.txt",header=T,sep="\t") 
head(dil)
eff=PrimEff(dil)

##turn ct into molecule counts
genecolumns=c(7:18) 
conditions=c(1:6)
counts=cq2counts(data=ct, genecols=c(7:18),condcols=c(1:6),effic=eff)
head(counts)
write.table(counts, "counts for samples with all genes only.csv", col.names = NA, sep = ",")

###linear mixed model 1 
model1=mcmc.qpcr(data=counts, fixed="strained+APA+APA:strained",random=c("sample","time"),pr=TRUE,pl=TRUE)
summary(model1)
diagnostic.mcmc(model1)
#plota=HPDsummary(model=model1,data=counts)
#poltb=HPDsummary(model=model1,data=counts,relative=TRUE)
print1=getNormalizedData(model1,counts)
data<-data.frame(cbind(print1$normData,print1$conditions))
#remove negative numbers
data$fos[data$fos<0]=0
data$creb[data$creb<0]=0
data$fmr1[data$fmr1<0]=0
#remove NAs
data$fmr1[is.na(data$fmr1)] = 0
data$cam2kd[is.na(data$cam2kd)] = 0


attach(data)
#plot individual data
par(mfrow=c(4,3))
plot(ind, cam2kd, main="cam2kd", ylab="counts")
plot(ind, creb, main="creb", ylab="counts")
plot(ind, dlg4, main="dlg4", ylab="counts")
plot(ind, fmr1, main="fmr1", ylab="counts")
plot(ind, fos, main="fos", ylab="counts")
plot(ind, gria, main="gria", ylab="counts")
plot(ind, grim, main="grim", ylab="counts")
plot(ind, grin, main="grin", ylab="counts")
plot(ind, nsf, main="nsf", ylab="counts")
plot(ind, pkmz, main="pkmz", ylab="counts")
plot(ind, rpl19, main="rpl19", ylab="counts")
plot(ind, rRNA18S, main="rRNA18S", ylab="counts")


#big group differences
par(mfrow=c(4,3))
boxplot(cam2kd~strained,data=data, main="cam2kd", ylab="counts")
boxplot(creb~strained,data=data, main="creb", ylab="counts")
boxplot(dlg4~strained,data=data, main="dlg4", ylab="counts")
boxplot(fmr1~strained,data=data, main="fmr1", ylab="counts")
boxplot(fos~strained,data=data, main="fos", ylab="counts")
boxplot(gria~strained,data=data, main="gria", ylab="counts")
boxplot(grim~strained,data=data, main="grim", ylab="counts")
boxplot(grin~strained,data=data, main="grin", ylab="counts")
boxplot(nsf~strained,data=data, main="nsf", ylab="counts")
boxplot(pkmz~strained,data=data, main="pkmz", ylab="counts")
boxplot(rpl19~strained,data=data, main="rpl19", ylab="counts")
boxplot(rRNA18S~strained,data=data, main="rRNA18S", ylab="counts")
#big group differences
boxplot(cam2kd~APA,data=data, main="cam2kd", ylab="counts")
boxplot(creb~APA,data=data, main="creb", ylab="counts")
boxplot(dlg4~APA,data=data, main="dlg4", ylab="counts")
boxplot(fmr1~APA,data=data, main="fmr1", ylab="counts")
boxplot(fos~APA,data=data, main="fos", ylab="counts")
boxplot(gria~APA,data=data, main="gria", ylab="counts")
boxplot(grim~APA,data=data, main="grim", ylab="counts")
boxplot(grin~APA,data=data, main="grin", ylab="counts")
boxplot(nsf~APA,data=data, main="nsf", ylab="counts")
boxplot(pkmz~APA,data=data, main="pkmz", ylab="counts")
boxplot(rpl19~APA,data=data, main="rpl19", ylab="counts")
boxplot(rRNA18S~APA,data=data, main="rRNA18S", ylab="counts")

#effects of training and strain
par(mfrow=c(2,2))
boxplot(cam2kd~APA*strained, data=data, notch=FALSE,col=(c("cyan","yellow", "deeppink")),
        main="cam2kd", xlab="Strain and Training")
boxplot(creb~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="creb", xlab="Strain and Training")
boxplot(dlg4~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="dlg4", xlab="Strain and Training")
boxplot(fmr1~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="fmr1", xlab="Strain and Training")
boxplot(fos~APA*strained, data=data, notch=FALSE,   col=(c("cyan","yellow", "deeppink")),
        main="fos", xlab="Strain and Training")
boxplot(gria~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="gria", xlab="Strain and Training")
boxplot(grim~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="grim", xlab="Strain and Training")
boxplot(grin~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="grin", xlab="Strain and Training")
boxplot(nsf~APA*strained, data=data, notch=FALSE,  col=(c("cyan","yellow", "deeppink")),
        main="nsf", xlab="Strain and Training")
boxplot(pkmz~APA*strained, data=data, notch=FALSE, col=(c("cyan","yellow", "deeppink")),
        main="pkmz", xlab="Strain and Training")
boxplot(rpl19~APA*strained, data=data, notch=FALSE, col=(c("cyan","yellow", "deeppink")),
        main="rpl19", xlab="Strain and Training")
boxplot(rRNA18S~APA*strained, data=data, notch=FALSE, col=(c("cyan","yellow", "deeppink")),
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