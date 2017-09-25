## integrative script


# animals details from sequenencing file
head(TraitsTrained)

# beahvior from pca
head(scoresdf) 
scoresdf$T4Path2Entr <- forpca$T4_C1_Path2ndEntr

#Combine the two, drop all buy first 5 pcs
FactorsBeahavPCA <- inner_join(TraitsTrained, scoresdf, by="ID")
FactorsBeahavPCA <- FactorsBeahavPCA[-c(13:24)]

ggplot(FactorsBeahavPCA, aes(x = T4Path2Entr, y = PC5, colour = factor(APA.y))) + 
  geom_point(size = 4)






