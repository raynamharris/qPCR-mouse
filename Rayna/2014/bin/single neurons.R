##install MCMCglmm and MCMC.qpcr
##load libraries
library("MCMCglmm")
library("MCMC.qpcr")
setwd("Z:/NSB_2014/Course Developers/Hofmann Fenton Collaboration/single neurons")

##read sample cts
ct=read.table("single neurons no oriens.txt",header=T,sep="\t") 
head(ct)

##read standard curves, calculate efficiencies
dil=read.table("dilutions.txt",header=T,sep="\t") 
head(dil)
eff=PrimEff(dil)

##turn ct into molecule counts
genecolumns=c(8:10) 
conditions=c(1:7)
counts=cq2counts(data=ct, genecols=genecolumns,condcols=conditions,effic=eff)
head(counts)
counts$count=counts$count*18
write.table(counts, "counts for single neurons no oriens.csv", col.names = NA, sep = ",")

###linear mixed model 1 
model1=mcmc.qpcr(data=counts, fixed="neurons+regions+neurons:regions",random=c("ID","strain","mouse"),pr=TRUE,pl=TRUE)
summary(model1)
diagnostic.mcmc(model1)
plota=HPDsummary(model=model1,data=counts)
poltb=HPDsummary(model=model1,data=counts,relative=TRUE)


print1=getNormalizedData(model1,counts)
data<-data.frame(cbind(print1$normData,print1$conditions))
#remove NAs
data$grin[is.na(data$grin)] = 0
#remove negative numbers
data$fos[data$fos<0]=0
data$grin[data$grin<0]=0

data$inverse.fos=2^data$fos
data$inverse.gria=2^data$gria
data$inverse.grin=2^data$grin

attach(data)
#plot individual data
par(mfrow=c(2,2))
plot(ID, fos, main="fos", ylab="log2(abundance)")
plot(ID, gria, main="gria", ylab="log2(abundance)")
plot(ID, grin, main="grin", ylab="log2(abundance)")

#plot individual data
par(mfrow=c(2,2))
plot(neurons, fos, main="fos", ylab="log2(abundance)", xlab="# of pooled neurons", cex.lab=1.2, ylim=c(0,8))
plot(neurons, gria, main="gria", ylab="log2(abundance)", xlab="# of pooled neurons", cex.lab=1.2, ylim=c(0,8))
plot(neurons, grin, main="grin", ylab="log2(abundance)", xlab="# of pooled neurons", cex.lab=1.2, ylim=c(0,8))


#big group differences
par(mfrow=c(2,2))
boxplot(fos~regions,data=data, main="fos", ylab="log2(abundance)", xlab="neural population", cex.lab=1.2, ylim=c(0,8))
boxplot(gria~regions,data=data, main="gria", ylab="log2(abundance)", xlab="neural population",cex.lab=1.2, ylim=c(0,8))
boxplot(grin~regions,data=data, main="grin", ylab="log2(abundance)", xlab="neural population",cex.lab=1.2, ylim=c(0,8))



