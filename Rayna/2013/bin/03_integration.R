## integrative analysis ----

# check col names of datasets or merge
names(wtfmr1_long)
names(dd_nohomecage)
names(summary_long)

# make grouping column that is genoAPAind
wtfmr1_long$grouping <- as.factor(paste(wtfmr1_long$genoAPA,wtfmr1_long$ind, sep="_"))
dd_nohomecage$grouping <- as.factor(paste(dd_nohomecage$genoAPA,dd_nohomecage$ind, sep="_"))
summary_long$grouping <- as.factor(paste(summary_long$genoAPA,summary_long$ind, sep="_"))

# start with qpcr data, slim to just grouping, value, and variable
integrative <- dd_nohomecage %>%
  dplyr::rename(variable = gene, value=count) %>%
  dplyr::select(grouping, genoAPA, variable, value)
names(integrative) 
str(integrative)

# then with behav data, slim to just grouping, value, and variable
wtfmr1_long_forintegrative <- wtfmr1_long %>%
  dplyr::select(grouping, genoAPA, variable, value)
names(wtfmr1_long_forintegrative)
str(wtfmr1_long_forintegrative)


## then with the summary data
summary_long_forintegrative <- summary_long %>%
  filter(grepl("Total|IO_Max", variable))%>% 
  dplyr::select(grouping, genoAPA, variable, value)
names(summary_long_forintegrative)
str(summary_long_forintegrative)
tail(summary_long_forintegrative)

## joing data, then rename factors
integrative <-  dplyr::bind_rows(integrative, wtfmr1_long_forintegrative, summary_long_forintegrative)
str(integrative)
integrative$variable <- as.factor(integrative$variable)

## widen the data, make it a matrix, make it a cor matrix
head(integrative)
integrative_wide <- dcast(integrative, grouping +genoAPA ~ variable, value.var= "value", fun.aggregate=mean)
names(integrative_wide)
rownames(integrative_wide) <- integrative_wide$grouping  # set $gropuing as rownames
integrative_wide_matrix <- integrative_wide[-c(1,2,62)] # remove grouping, genoAPA, and TotalTime
integrative_wide_matrix <- as.matrix(integrative_wide_matrix)
integrative_wide_matrix <- na.omit(integrative_wide_matrix)
head(integrative_wide_matrix)
integrative_wide_matrix_cor <- round(cor(integrative_wide_matrix),2) 

names(integrative_wide)
integrative_wide_matrix_slim <- integrative_wide[-c(1:13,16:18,20,23:30,36:47,49,62)]
integrative_wide_matrix_slim <- na.omit(integrative_wide_matrix_slim)
integrative_wide_matrix_slim_cor <- round(cor(integrative_wide_matrix_slim),2) 

## plot the cor heat map
col <- colorRampPalette(c("turquoise4", "white", "tan4"))(n = 299)
heatmap.2(integrative_wide_matrix_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none")
dev.off()
corrplot(integrative_wide_matrix_cor, type="upper", order="hclust", tl.col="black", tl.srt=45)


corrplot(integrative_wide_matrix_slim_cor, type="lower", order="hclust", tl.col="black", tl.srt=45)
heatmap.2(integrative_wide_matrix_slim_cor, col=col, 
          density.info="none", trace="none", dendrogram=c("row"), 
          symm=F,symkey=T,symbreaks=T, scale="none")


png('3_integrative_corrplot.png')
setwd("~/Github/qPCR-mouse/Rayna/2013/results")
corrplot(integrative_wide_matrix_cor, type="upper",
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
res1 <- cor.mtest(integrative_wide_matrix,0.60)
## specialized the insignificant value according to the significant level
corrplot(integrative_wide_matrix_cor, type="upper",
         order="hclust", tl.col="black", tl.srt=45, tl.cex = 0.5,
         p.mat = res1[[1]], insig = "blank")


## ploting one correlation at a time ----
head(integrative)
integrative_wide <- dcast(integrative, grouping + genoAPA~ variable, value.var= "value", fun.aggregate=mean)
head(integrative_wide)
ggplot(integrative_wide, aes(x = pkmz, y = grim, colour = genoAPA)) + 
  geom_point() +
  scale_x_log10() + scale_y_log10() +   
  geom_smooth(method=lm)

## igraph ----
library(igraph)
head(integrative_wide_matrix_cor)

cor_mat<-integrative_wide_matrix_cor 
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

