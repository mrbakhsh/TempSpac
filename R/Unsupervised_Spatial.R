##################################################################
######unsupervised learning & visualization of spatial data#######
##################################################################




library(ggfortify)
library(gplots)
library(ComplexHeatmap)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(caret)
library(factoextra)
library(cluster)
library(dendextend)
library(Rtsne)



#load the data
df <-
    fread(file.choose())

##################################################################
#################### 1. k-means clustering
##################################################################

set.seed(1)
#define the optimal number of cluster (e.g., WCC)
fviz_nbclust(df[,-c(1,8)], kmeans, method = "wss")
#perform k-mean clustering & visualization on PCA
autoplot(kmeans(df[,-c(1,8)], 3), data = df[,-c(1,8)])+
    theme_bw()

##################################################################
#################### 2. Fuzzy clustering
##################################################################

#perfrom fuzzy clustering & visualization on PCA
set.seed(1)
autoplot(fanny(df[,-c(1,8)], 3)) +
    theme_bw()

##################################################################
#################### 3. Hierarchical clustering
##################################################################


#Hierarchical clustering & data visulaization on heatmap
Heatmap(as.matrix(df[,-c(1,8)]),
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    cluster_columns = TRUE,
    clustering_method_columns = "complete",
    show_row_names = F,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "color Key"))

#or visualize the data as cluster dendogram
set.seed(2)
hr <- # compute dissimilarity matrix and cluster rows
    hclust(dist(df[,-c(1,8)], method="euclidean"),
        method="complete")

treeR <- #cluster dendogram
    as.dendrogram(hr, method = "average")
hclustCutree =
    cutree(hr, h=1) #cut tree at height of 1.0

#plot the dendogram
plot(treeR,
    leaflab = "none",
    main = "Gene Clustering",
    ylab = "Height")

#color the clusters
colored_bars(hclustCutree, treeR,
    sort_by_labels_order = T, y_shift=-0.1,
    rowLabels = c("h=3"),cex.rowLabels=0.7)
#this will add lines showing the cut heights
abline(h=1.0, lty = 2, col="grey")



######### Overlay marker on PCA plot ############
pca_res <- #keep the numeric data
    prcomp(df[,-c(1,8)], scale. = TRUE)
autoplot(pca_res, data = df[,-1], colour = "Marker") +
    theme_bw()


##################################################################
#################### 4. t-sne plot 
##################################################################

## Executing the algorithm on curated data
tsne_out <-
    Rtsne(df[,-c(1,8)],
        dims = 2, perplexity=30,#can be adjusted
        verbose=TRUE, max_iter = 500)

tsne_output <-
    as.data.frame(tsne_out$Y)

ggplot(tsne_output, aes(V1, V2, color=df$Marker)) +
    geom_point() +
    theme_bw() +
    labs(x= "Dim1", y="Dim2", title = "Preplexity = 30")





