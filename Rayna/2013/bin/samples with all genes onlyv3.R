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

###linear mixed model 3 
model3=mcmc.qpcr(data=counts, fixed="strain.APA",random=c("sample","time"),pr=TRUE,pl=TRUE)
summary(model3)
diagnostic.mcmc(model3)
plota=HPDsummary(model=model3,data=counts)
poltb=HPDsummary(model=model3,data=counts,relative=TRUE)
summary(plota)

cbPalette <- c("#000000", "#FF0000", "#999999", "#FF9933")
pd = position_dodge(0.3)
ggplot(
  plota$summary, 
  aes_string(x = "strain.APA", colour = "strain.APA",
             y = "mean")) + 
  geom_errorbar(aes(ymin = lower, ymax = upper), lwd = 1, width = .8, position = pd) + geom_point(position = pd, size = 3) + 
  theme_bw() +facet_wrap(~gene,scales="free_y",nrow=2) +ylab("") +
  theme(legend.position="bottom") +
  theme(strip.text.x=element_text(face="italic",size = 14, family="serif")) +
  # edit the next line according to the names of your factors:
    # remark this line to see if defalt colors are good for you, if not, edit as needed
  scale_color_manual(values=cbPalette)+
  theme(axis.title.y = element_text(size=25), panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
        panel.background = element_blank()
  )

print3=getNormalizedData(model3,counts)
datamodel3<-data.frame(cbind(print3$normData,print3$conditions))
write.csv(datamodel3, "datamodel3.csv")


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



