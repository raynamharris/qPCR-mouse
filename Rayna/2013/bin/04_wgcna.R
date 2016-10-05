# WGCNA ----
# Install and Load WGCNA package
#source("https://bioconductor.org/biocLite.R")
#biocLite("WGCNA")
install.packages("flashClust")
library(WGCNA)
library(flashClust)
options(stringsAsFactors=FALSE)
allowWGCNAThreads()

########################################################    
#                     Load data
########################################################    
# Load behavior data, melt, remove some rows, make sesssionbehavior column, then widen ----
datExpr0 <- wtfmr1 
datExpr0 <- melt(datExpr0, id=c("ind","genotype", "APA", "session", "genoAPA", "genoAPAsession", "genoAPAsessionInd", "filename"))
head(datExpr0)
datExpr0 <- filter(datExpr0, !grepl("p.miss|TotalTime", variable))
datExpr0$sessionbeahvior <- as.factor(paste(datExpr0$session, datExpr0$variable, sep="_"))
datExpr0 <- dcast(datExpr0, ind + genotype+ APA ~ sessionbeahvior, value.var= "value")
rownames(datExpr0) <- datExpr0$ind     # set $genoAPAsessionInd as rownames
datExpr0 <- datExpr0[-c(1:3)] #delete all non-numeric columns 
head(datExpr0)

gsg=goodSamplesGenes(datExpr0, verbose = 1)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts

#-----Load trait data and make all the factor integers
datTraits <- wtfmr1
datTraits <- melt(datTraits, id=c("ind","genotype", "APA", "session", "genoAPA", "genoAPAsession", "genoAPAsessionInd", "filename"))
head(datTraits)
datTraits$sessionbeahvior <- as.factor(paste(datTraits$session, datTraits$variable, sep="_"))
datTraits <- dcast(datTraits, ind + genotype + APA + genoAPA ~ sessionbeahvior, value.var= "value")
rownames(datTraits) <- datTraits$ind     # set $genoAPAsessionInd as rownames
names(datTraits)
datTraits <- datTraits[c(2:4)] #keep only trait columns 
head(datTraits)
str(datTraits)

## making it a numeric
datTraits$genotype <- as.integer(factor(datTraits$genotype))
datTraits$APA <- as.integer(factor(datTraits$APA))
datTraits$genoAPA <- as.integer(factor(datTraits$genoAPA))
str(datTraits)
head(datTraits)


#######   #################    ################   #######    
#                 Call sample outliers
#######   #################    ################   #######   

#-----Sample dendrogram and traits
A=adjacency(t(datExpr0),type="signed")
#-----Calculate whole network connectivity
k=as.numeric(apply(A,2,sum))-1
#-----Standardized connectivity
Z.k=scale(k)
thresholdZ.k=3 
outlierColor=ifelse(Z.k<thresholdZ.k,"red","black")
sampleTree = flashClust(as.dist(1-A), method = "average")
#-----Convert traits to colors
traitColors=data.frame(numbers2colors(datTraits,signed=FALSE))
str(traitColors)
dimnames(traitColors)[[2]]=paste(names(datTraits))
datColors=data.frame(outlier=outlierColor,traitColors)

#-----Plot the sample dendrogram
#quartz()
plotDendroAndColors(sampleTree,groupLabels=names(datColors),
                    colors=datColors,main="Sample dendrogram and trait heatmap")

#-----Remove outlying samples 
#remove.samples= Z.k<thresholdZ.k | is.na(Z.k)
#datExpr0=datExpr0[!remove.samples,]
#datTraits=datTraits[!remove.samples,]
#A=adjacency(t(datExpr0),type="distance")
#k=as.numeric(apply(A,2,sum))-1
#Z.k=scale(k)


#######   #################    ################   #######    
#                     Choose soft threshold
#######   #################    ################   #######     

dim(datExpr0)
dim(datTraits)
powers= c(seq(1,10,by=1), seq(from =12, to=20, by=2)) #choosing a set of soft-thresholding powers
sft = pickSoftThreshold(datExpr0, powerVector=powers, verbose =5,networkType="signed") #call network topology analysis function

sft <- pickSoftThreshold(
  datExpr0, 
  dataIsExpr = TRUE,
  RsquaredCut = 0.90, 
  powerVector = c(seq(1, 10, by = 1), seq(12, 20, by = 2)), 
  removeFirst = FALSE, nBreaks = 10, blockSize = NULL, 
  corFnc = cor, corOptions = list(use = 'p'), 
  networkType = "unsigned",
  moreNetworkConcepts = FALSE,
  verbose = 0, indent = 0)
sft

#quartz()
par(mfrow= c(1,2))
cex1=0.9
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], xlab= "Soft Threshold (power)", ylab="Scale Free Topology Model Fit, signed", type= "n", main= paste("Scale independence"))
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2], labels=powers, cex=cex1, col="red")
abline(h=0.90, col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5], xlab= "Soft Threshold (power)", ylab="Mean Connectivity", type="n", main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1, col="red")
dev.off()
#softPower=24

#######   #################    ################   #######    
#                    Construct network
#######   #################    ################   #######     

adjacency=adjacency(datExpr0, type="signed") 
TOM= TOMsimilarity(adjacency, TOMType="signed")
dissTOM= 1-TOM

geneTree= flashClust(as.dist(dissTOM), method="average")

#quartz()
dev.off()
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE, hang=0.04)

#######   #################    ################   #######    
#                    Make modules
#######   #################    ################   ####### 

minModuleSize=9 
dynamicMods= cutreeDynamic(dendro= geneTree, distM= dissTOM, deepSplit=2, pamRespectsDendro= FALSE, minClusterSize= minModuleSize)
table(dynamicMods)

dynamicColors= labels2colors(dynamicMods)

#quartz()
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut", dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang= 0.05, main= "Gene dendrogram and module colors")

#-----Merge modules whose expression profiles are very similar
MEList= moduleEigengenes(datExpr0, colors= dynamicColors)
MEs= MEList$eigengenes
#Calculate dissimilarity of module eigenegenes
MEDiss= 1-cor(MEs)
#Cluster module eigengenes
METree= flashClust(as.dist(MEDiss), method= "average")

#quartz()
plot(METree, main= "Clustering of module eigengenes", xlab= "", sub= "")
MEDissThres= 0.3
abline(h=MEDissThres, col="red")
merge= mergeCloseModules(datExpr0, dynamicColors, cutHeight= MEDissThres, verbose =3)

mergedColors= merge$colors
mergedMEs= merge$newMEs

#quartz()
plotDendroAndColors(geneTree, cbind(dynamicColors, mergedColors), c("Dynamic Tree Cut", "Merged dynamic"), dendroLabels= FALSE, hang=0.03, addGuide= TRUE, guideHang=0.05)


moduleColors= mergedColors
colorOrder= c("grey", standardColors(50))
moduleLabels= match(moduleColors, colorOrder)-1
MEs=mergedMEs

#save(MEs, moduleLabels, moduleColors, geneTree, file= "SamplesAndColors_thresh24merge42_signed.RData")



#######   #################    ################   #######    
#                Relate modules to traits
#######   #################    ################   ####### 

datt=datExpr0

#-----Define numbers of genes and samples
nGenes = ncol(datt);
nSamples = nrow(datt);
#-----Recalculate MEs with color labels
MEs0 = moduleEigengenes(datt, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#-----Correlations of genes with eigengenes
moduleGeneCor=cor(MEs,datt)
moduleGenePvalue = corPvalueStudent(moduleGeneCor, nSamples);

moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);


#---------------------Module-trait heatmap

#quartz()
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.7,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
######--------------------end--------------------#######



#---------------------Gene significance by Module membership scatterplots
whichTrait="APA" #Replace this with the trait of interest

#quartz()
nGenes = ncol(datt);
nSamples = nrow(datt);
selTrait = as.data.frame(datTraits[,whichTrait]);
names(selTrait) = whichTrait
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(signedKME(datt, MEs));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(datt, selTrait, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(selTrait), sep="");
names(GSPvalue) = paste("p.GS.", names(selTrait), sep="");
par(mfrow=c(2,3))
counter=0
for(module in modNames[1:length(modNames)]){
  counter=counter+1
  if (counter>6) {
    quartz()
    par(mfrow=c(2,3))
    counter=1
  }
  column = match(module, modNames);
  moduleGenes = moduleColors==module;
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                     abs(geneTraitSignificance[moduleGenes, 1]),
                     xlab = paste(module,"module membership"),
                     ylab = paste("GS for", whichTrait),
                     col = module,mgp=c(2.3,1,0))
}
## saved as 4-MM-APA for APA modules
######--------------------end--------------------#######



#---------------------Eigengene heatmap
which.module="blue" #replace with module of interest
datME=MEs
datExpr=datt
#quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", names.arg=c(row.names(datt)), cex.names=0.5, cex.main=2,
        ylab="eigengene expression",xlab="sample")
######--------------------end--------------------#######


#### output "gene" list ----
## see https://github.com/ClaireMGreen/TDP-43_Code/blob/afb43cddb8ec1a940fbcfa106a1cc3cf77568b7e/WGCNA2.R

#find out what the IDs are of the genes that are contained within a module. 

blue <- as.data.frame(colnames(datExpr0)[moduleColors=='blue'])
blue$module <- "blue"
colnames(blue)[1] <- "sessionbeahvior"
red <- as.data.frame(colnames(datExpr0)[moduleColors=='red'])
red$module <- "red"
colnames(red)[1] <- "sessionbeahvior"
green <- as.data.frame(colnames(datExpr0)[moduleColors=='green'])
green$module <- "green"
colnames(green)[1] <- "sessionbeahvior"
yellow <- as.data.frame(colnames(datExpr0)[moduleColors=='yellow'])
yellow$module <- "yellow"
colnames(yellow)[1] <- "sessionbeahvior"
brown <- as.data.frame(colnames(datExpr0)[moduleColors=='brown'])
brown$module <- "brown"
colnames(brown)[1] <- "sessionbeahvior"
turquoise <- as.data.frame(colnames(datExpr0)[moduleColors=='turquoise'])
turquoise$module <- "turquoise"
colnames(turquoise)[1] <- "sessionbeahvior"



MM <- dplyr::bind_rows(blue,red,yellow,green,brown, turquoise)

MM$session <- ifelse(grepl("pretraining", MM$sessionbeahvior), "pretraining", 
                         ifelse(grepl("training1", MM$sessionbeahvior), "training1",
                                ifelse(grepl("training2", MM$sessionbeahvior), "training2",
                                       ifelse(grepl("training3", MM$sessionbeahvior), "training3",
                                              ifelse(grepl("retention", MM$sessionbeahvior), "retention", "NA")))))
MM$session  ## check that all names good with no NAs                                 
MM$session <- as.factor(MM$session)  
MM$session <- factor(MM$session, levels = c("pretraining", "training1", "training2", "training3", "retention"))

head(MM)
MM_session <- dcast(MM, module ~ session, value.var = 'session')

