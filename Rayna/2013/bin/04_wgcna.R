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
datExpr0 <- all 
rownames(datExpr0) <- datExpr0$ind     # set $genoAPAsessionInd as rownames
datExpr0 <- datExpr0[-c(1:4)] #delete all non-numeric columns 
head(datExpr0)
names(datExpr0)

gsg=goodSamplesGenes(datExpr0, verbose = 1)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts
head(gsg)

## removing bad measures
gsg
datExpr0 <- datExpr0[-c(9:10)]  #removes LTP measures with excessive NAs
gsg=goodSamplesGenes(datExpr0, verbose = 1)
gsg$allOK #If the last statement returns TRUE, all genes have passed the cuts



#-----Make a trait data frame
datTraits <- all
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
MEDissThres= 0.4
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
dev.off()
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
##saved as 4-Moduletrait-alldata.png
######--------------------end--------------------#######



#---------------------Eigengene heatmap
which.module="red" #replace with module of interest
datME=MEs
datExpr=datt
#quartz()
ME=datME[, paste("ME",which.module, sep="")]
par(mfrow=c(2,1), mar=c(0.3, 5.5, 3, 2))
plotMat(t(scale(datExpr[,moduleColors==which.module ]) ),
        nrgcols=30,rlabels=F,rcols=which.module,
        main=which.module, cex.main=2)
par(mar=c(5, 4.2, 0, 0.7))
barplot(ME, col=which.module, main="", names.arg=(datTraits$genoAPA), cex.names=0.5, cex.main=2,
        ylab="eigengene expression",xlab="sample")

## saved as 4-<color>-alldata 

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
black <- as.data.frame(colnames(datExpr0)[moduleColors=='black'])
black$module <- "black"
colnames(black)[1] <- "sessionbeahvior"

## merged data frame 
MM <- dplyr::bind_rows(blue,turquoise, brown,yellow,red, black)

## new column to pull out major data structures
MM$session <- ifelse(grepl("pretraining", MM$sessionbeahvior), "pretraining", 
                         ifelse(grepl("training1", MM$sessionbeahvior), "training1",
                                ifelse(grepl("training2", MM$sessionbeahvior), "training2",
                                       ifelse(grepl("training3", MM$sessionbeahvior), "training3",
                                              ifelse(grepl("retention", MM$sessionbeahvior), "retention", 
                                                     ifelse(grepl("rRNA18S|rpl19|prkcz|nsf|grin|gria|dlg4|cam2kd", MM$sessionbeahvior), "genes", 
                                                            ifelse(grepl("IO_Max", MM$sessionbeahvior), "physiol",
                                                                   ifelse(grepl("Total", MM$sessionbeahvior), "Totals",
                                                                          ifelse(grepl("Last|Punishment|Retention", MM$sessionbeahvior), "retention","other")))))))))
MM$session  ## check that all names good with no NAs                                 
MM$session <- as.factor(MM$session)  
MM_session <- dcast(MM, module ~ session, value.var = 'session')
head(MM_session)
blue
turquoise
brown
yellow
red
black

#write.csv(MM_session, "MM_session.csv" , row.names = F)
