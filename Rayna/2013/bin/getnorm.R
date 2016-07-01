library(MCMC.qpcr)
setwd("Z:/Shared/All projects/X. nigrensis mc4r/genotyping qpcr/data analysis")
dil=read.table("2014-02-14 analysis.dil.txt",header=T,sep="\t") 
eff=PrimEff(dil)
ct=read.table("2014-04-28 analysis.ct.txt",header=T,sep="\t") 
co=cq2counts(data=ct,genecols=c(8:12),condcols=c(1:7),effic=eff) 

##trying misha's soft model
soft=mcmc.qpcr(fixed="sex", random="ind", data=co, controls=c("xsrc", "ARa"), normalize=T)

#run the model
informed=mcmc.qpcr(fixed="sex", random="ind", data=co, controls=c("xsrc", "ARa"), m.fix=1, pr=TRUE)
summary(informed)
pp=getNormalizedData(informed,co)
data<-data.frame(cbind(pp$normData,pp$conditions))

#normalize by xsrc
data$a.norm=data$a*2 / (data$xsrc)  
data$b1.norm=data$b1*2 / (data$xsrc)
data$b2.norm=data$b2*2 / (data$xsrc)
data$ARa.norm=data$ARa / (data$xsrc)
#remove negative numbers
data$a.norm[data$a.norm<0]=0
data$b1.norm[data$b1.norm<0]=0
data$b2.norm[data$b2.norm<0]=0
#remove NAs
data$b1.norm[is.na(data$b1.norm)] = 0
data$b2.norm[is.na(data$b2.norm)] = 0
#calculate total b
data$total.b=data$b1.norm + data$b2.norm
data$b.over.a=data$total.b / data$a.norm 


#save the datafile
write.table(data, "data.04.28.14.csv", col.names = NA, sep = "," )

#make pretty graphs
attach(data)
par(mfrow=c(2,2))
plot(size.dis, total.b, main="group vs. total b allele", ylab="total b alleles")
plot(sex, total.b, main="group vs. total b allele", ylab="total b alleles")
plot(total.b, size.cont, xlim = NULL, ylim = NULL, main="body size vs. total b allele", ylab="body size mm", pch=19, )
detach(data)
