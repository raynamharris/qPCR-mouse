## integrative analysis ----

# check col names of datasets or merge
names(wtfmr1_long)
names(dd_nohomecage)

# make grouping column that is genoAPAind
wtfmr1_long$grouping <- as.factor(paste(wtfmr1_long$genoAPA,wtfmr1_long$ind, sep="_"))
dd_nohomecage$grouping <- as.factor(paste(dd_nohomecage$genoAPA,dd_nohomecage$ind, sep="_"))


# start with qpcr data, slim to just grouping, value, and variable
integtive <- dd_nohomecage %>%
  dplyr::rename(variable = gene, value=count) %>%
  dplyr::select(grouping, genoAPA, variable, value)
names(integtive) 
str(integtive)

# then with behav/ephys data, slim to just grouping, value, and variable
wtfmr1_long_forintegtive <- wtfmr1_long %>%
  dplyr::select(grouping, genoAPA, variable, value)
names(wtfmr1_long_forintegtive)
str(wtfmr1_long_forintegtive)

## joing data, then rename factors
integtive <-  dplyr::bind_rows(integtive, wtfmr1_long_forintegtive)
str(integtive)
integtive$variable <- as.factor(integtive$variable)

## widen the data, make it a matrix, make it a cor matrix
head(integtive)
integtive_wide <- dcast(integtive, grouping +genoAPA ~ variable, value.var= "value", fun.aggregate=mean)
names(integtive_wide)
rownames(integtive_wide) <- integtive_wide$grouping  # set $gropuing as rownames
integtive_wide_matrix <- integtive_wide[-c(1,2,58)] # remove grouping, genoAPA, and TotalTime
integtive_wide_matrix <- as.matrix(integtive_wide_matrix)
integtive_wide_matrix <- na.omit(integtive_wide_matrix)
head(integtive_wide_matrix)
integtive_wide_matrix_cor <- round(cor(integtive_wide_matrix),2) 
head(integtive_wide_matrix_cor)

## plot the cor heat map
col <- colorRampPalette(c("turquoise4", "white", "tan4"))(n = 299)
heatmap.2(integtive_wide_matrix_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none")
corrplot(integtive_wide_matrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)

png('3_integrative_corrplot.png')
setwd("~/Github/qPCR-mouse/Rayna/2013/results")
corrplot(integtive_wide_matrix_cor, type="upper",
         order="hclust", tl.col="black", tl.srt=45)
setwd("~/Github/qPCR-mouse/Rayna/2013/bin")
dev.off()

## plot with signficicance

## first define a new function
cor.mtest <- function(mat, conf.level = 0.60){
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat <- lowCI.mat <- uppCI.mat <- matrix(NA, n, n)
  diag(p.mat) <- 0
  diag(lowCI.mat) <- diag(uppCI.mat) <- 1
  for(i in 1:(n-1)){
    for(j in (i+1):n){
      tmp <- cor.test(mat[,i], mat[,j], conf.level = conf.level)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
      lowCI.mat[i,j] <- lowCI.mat[j,i] <- tmp$conf.int[1]
      uppCI.mat[i,j] <- uppCI.mat[j,i] <- tmp$conf.int[2]
    }
  }
  return(list(p.mat, lowCI.mat, uppCI.mat))
}

## caluculate residguals
res1 <- cor.mtest(integtive_wide_matrix,0.60)
## specialized the insignificant value according to the significant level
corrplot(integtive_wide_matrix_cor, type="upper",
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.5,
         p.mat = res1[[1]], insig = "blank")


## ploting one correlation at a time ----
head(integtive)
integtive_wide <- dcast(integtive, grouping + genoAPA~ variable, value.var= "value", fun.aggregate=mean)
head(integtive_wide)
ggplot(integtive_wide, aes(x = pkmz, y = grim, colour = genoAPA)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm)

## igraph ----
library(igraph)
head(integtive_wide_matrix_cor)

cor_mat<-integtive_wide_matrix_cor 
diag(cor_mat)<-0
graph<-graph.adjacency(cor_mat, weighted=TRUE, mode="upper") 
E(graph)[ weight>0.9 ]$color <- "orangered4" 
E(graph)[ weight>0.80 & weight < 0.899]$color <- "orangered2" 
#E(graph)[ weight<-0.9 ]$color <- "blue" 
#E(graph)[ weight>-0.80 & weight < -0.899]$color <- "blue" 
E(graph)[ weight<-0.8 ]$color <- "black" 

# removed edges and nodes below a certain value
graph <- delete.edges(graph, E(graph)[ weight < 0.8 & weight > -0.8])
graph <- delete.vertices(graph,which(degree(graph)<1))

## color codes nodes
V(graph)$name  #prints column/vertice names
V(graph)[c(6:8,11:15,20,25,33:34)]$color <- "tan2"
V(graph)[c(9:10,16:17,22:24,36,38:42)]$color <- "turquoise2"

## default plot
plot(graph)

## gives coords to plot in a circle then plots
coords <- layout_(graph, in_circle())
plot(graph, layout = coords)

reingold <- layout.reingold.tilford(graph, circular=T)
plot(graph, layout = reingold)



g <- make_ring(10) %>%
  set_vertex_attr("name", value = LETTERS[1:10])
g
V(g)

g2 <- delete_vertices(graph, c(1:5,18:19,21,26:32,37)) 
V(g2)$name  #prints column/vertice names
V(g2)[c(1:3,6:10,13,17:19)]$color <- "tan2"
V(g2)[c(4:5,11:12,14:16,20:26)]$color <- "turquoise2"

g2 <- delete.edges(g2, E(g2)[ weight < 0.8 & weight > -0.8])
g2 <- delete.vertices(g2,which(degree(g2)<1))

E(g2)[ weight > 0 ]$color <- "red" 
E(g2)[ weight < 0 ]$color <- "blue" 

plot(g2, layout = reingold)
plot(g2, layout = coords)

