##install MCMCglmm and MCMC.qpcr
##load libraries
library("MCMC.qpcr")

setwd("Y:/NSB_2014/3_Int Mol Neuro/R analysis/Rayna's examples")

##read sample cts
ct=read.csv("CA1 samples 3 genes only V2.csv") 
head(ct)

##read standard curves, calculate efficiencies
dil=read.table("dilutions.txt",header=T,sep="\t") 
head(dil)
eff=PrimEff(dil)

#turn ct into molecule counts
genecolumns=c(9:11) 
conditions=c(1:8)
counts=cq2counts(data=ct, genecols=genecolumns,condcols=conditions,effic=eff)
head(counts)
#control for dilution
if(counts$year=="2013"){counts$newcount=counts$count*2}
if(counts$year=="2014"){counts$newcount=counts$count*1}

hist(counts$count)
counts$logcount=log(counts$count)
hist(counts$logcount)

head(counts) 
counts$count=counts$newcount
counts$newcount=NULL
counts$logcount=NULL
head(counts) 

###uninformed linear mixed model 1 
model1=mcmc.qpcr(data=counts, fixed="strained+APA+strained:APA",pr=TRUE,pl=TRUE)
summary(model1)
diagnostic.mcmc(model1)
plota=HPDsummary(model=model1,data=counts)
poltb=HPDsummary(model=model1,data=counts,relative=TRUE)

model2=mcmc.qpcr(data=counts, fixed="APA+strained+APA:strained",pr=TRUE,pl=TRUE)
summary(model2)
diagnostic.mcmc(model2)
plota=HPDsummary(model=model2,data=counts)
poltb=HPDsummary(model=model2,data=counts,relative=TRUE)

model3=mcmc.qpcr(data=counts, fixed="APA+strained+APA:strained",pr=TRUE,pl=TRUE, controls = c("rpl19"), normalize=TRUE)
summary(model3)
diagnostic.mcmc(model3)
plota=HPDsummary(model=model3,data=counts)
poltb=HPDsummary(model=model3,data=counts,relative=TRUE)

model4=mcmc.qpcr(data=counts, fixed="strained+APA+strained:APA",pr=TRUE,pl=TRUE, controls = c("rpl19"), normalize=TRUE)
summary(model4)
diagnostic.mcmc(model4)
plota=HPDsummary(model=model4,data=counts)
poltb=HPDsummary(model=model4,data=counts,relative=TRUE)

print1=getNormalizedData(model1,counts)
data<-data.frame(cbind(print1$normData,print1$conditions))
#remove NAs
data$pkmz[is.na(data$pkmz)] = 0
#remove negative numbers
data$pkmz[data$pkmz<0]=0
data$grim[data$grim<0]=0


attach(data)
aov.grim.data = aov(grim~strained*APA,data=data)  
summary(aov.grim.data)
aov.pkmz.data = aov(pkmz~strained*APA,data=data)  
summary(aov.pkmz.data)
aov.rpl19.data = aov(rpl19~strained*APA,data=data)  
summary(aov.rpl19.data)

TukeyHSD(aov.grim.data)
TukeyHSD(aov.pkmz.data)
TukeyHSD(aov.rpl19.data)

par(mfrow=c(2,2))
boxplot(grim~APA*strained, data=data, notch=FALSE,  col=(c("gray30","gray90", "gray30","gray90")),
        main="grim", ylab="log2(abundance)", cex.axis=1.3, names=c("B6 untrained", "B6 trained", "fmr1- untrained","fmr1- trained"))
boxplot(pkmz~APA*strained, data=data, notch=FALSE, col=(c("gray30","gray90", "gray30","gray90")),
        main="pkmz", ylab="log2(abundance)", cex.axis=1.3, names=c("B6 untrained", "B6 trained", "fmr1- untrained","fmr1- trained"))
boxplot(rpl19~APA*strained, data=data, notch=FALSE, col=(c("gray30","gray90", "gray30","gray90")),
        main="rpl19", ylab="log2(abundance)",cex.axis=1.3, names=c("B6 untrained", "B6 trained", "fmr1- untrained","fmr1- trained"))
text(1, 15.8,"A")
text(2, 15.8, "A")
text(3, 15.8, "B")
text(4, 15.8, "B")
plot.new()
legend(0.1,1, c("untrained", "trained"), fill = c("gray30", "gray90"), title="APA", cex=1.5)

#examining correlations
pairs(~grim+pkmz+rpl19,data=data) 
pairs(~grim+pkmz,data=data) 
pairs(~rpl19+pkmz,data=data) 

cor.test(data$grim, data$pkmz, method = "spearman")
cor.test(data$grim, data$pkmz, method = "pearson")
cor.test(data$grim, data$rpl19, method = "spearman")
cor.test(data$grim, data$rpl19, method = "pearson")
cor.test(data$pkmz, data$rpl19, method = "spearman")
cor.test(data$pkmz, data$rpl19, method = "pearson")
detach(data)


print2=getNormalizedData(model2,counts)
data2<-data.frame(cbind(print1$normData,print2$conditions))
#remove NAs
data$pkmz[is.na(data$pkmz)] = 0
#remove negative numbers
data$pkmz[data$pkmz<0]=0
data$grim[data$grim<0]=0

attach(data2)
par(mfrow=c(2,2))
boxplot(grim~APA*strained, data=data2, notch=FALSE,  col=(c("gray30","gray90", "gray30","gray90")),
        main="grim", ylab="log2(abundance)", cex.axis=1.3, names=c("B6 untrained", "B6 trained", "fmr1- untrained","fmr1- trained"))
boxplot(pkmz~APA*strained, data=data2, notch=FALSE, col=(c("gray30","gray90", "gray30","gray90")),
        main="pkmz", ylab="log2(abundance)", cex.axis=1.3, names=c("B6 untrained", "B6 trained", "fmr1- untrained","fmr1- trained"))
boxplot(rpl19~APA*strained, data=data2, notch=FALSE, col=(c("gray30","gray90", "gray30","gray90")),
        main="rpl19", ylab="log2(abundance)",cex.axis=1.3, names=c("B6 untrained", "B6 trained", "fmr1- untrained","fmr1- trained"))
pairs(~grim+pkmz+rpl19,data=data2) 
pairs(~grim+pkmz,data=data2) 
pairs(~rpl19+pkmz,data=data2)
counts$logcount=log(counts$count)
detach(data2)

counts$logcount=log(counts$count)

write.csv(data2, "data2.10.1.14.csv" )
write.csv(data, "data1.10.1.14.csv" )

