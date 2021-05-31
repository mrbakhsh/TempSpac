##################################################################
######unsupervised learning & visualization of temporal data######
##################################################################


library(gplots)
library(ComplexHeatmap)
library(circlize)
library(Mfuzz)
library(ggplot2)
library(data.table)
library(dplyr)
library(tidyr)
library(reshape2)
library(dendextend)
library(kohonen)
library(magrittr)
library(factoextra)



df <-
    read.table(file.choose("Temporal_SampleData.txt"),
        sep="\t", header =T, row.names = 1)
scaledata <-
    t(scale(t(df))) #center and scale data input (optional)

##################################################################
#################### 1. Hierarchical Clustering
##################################################################

# Hierarchical clustering & visualizing the result as heatmap 
col_fun = #define the color
    colorRamp2(c(-2,-1, 0,1,2),
        c("#08306b","#08519c","black","#f4e52a","#f2e216"))
Heatmap(as.matrix(scaledata),
    cluster_rows = TRUE,
    clustering_distance_rows = "euclidean",
    clustering_method_rows = "complete",
    cluster_columns = TRUE,
    clustering_method_columns = "complete",
    show_row_names = F,
    show_heatmap_legend = TRUE,
    heatmap_legend_param = list(title = "color Key"),
    col = col_fun)

# Hierarchical clustering using basic functions
set.seed(222)
hr <- # compute dissimilarity matrix and cluster rows
    hclust(dist(scaledata, method="euclidean"),
        method="complete")

treeR <- #cluster dendogram
    as.dendrogram(hr, method = "average")
hclustCutree =
    cutree(hr, h=3.9) #cut tree at height of 1.0

#plot the dendogram
plot(treeR,
    leaflab = "none",
    main = "Gene Clustering",
    ylab = "Height")

#color the clusters
colored_bars(hclustCutree, treeR,
    sort_by_labels_order = T, y_shift=-0.1,
    rowLabels = c("h=3.9"),cex.rowLabels=0.7)
#add line showing the cut heights
abline(h=3.9, lty = 2, col="grey")


#plotting the centroids

extClust <- #extract the cluster
    hclustCutree


# function to find centroid in cluster i
ClustCent = function(i, dat, clusters) {
    ind = (clusters == i)
    colMeans(dat[ind,])
}
kClustcentroids <-
    sapply(levels(factor(extClust)),
        ClustCent, scaledata, extClust)

hc_longdf <-
    reshape2::melt(kClustcentroids)
colnames(hc_longdf) <-
    c('condition','cluster','value')

#plot
hc_longdf$condition <-
    factor(hc_longdf$condition,
        levels=unique(hc_longdf$condition))

ggplot(hc_longdf,
    aes(x=condition,y=value,
        group=cluster,
        colour=as.factor(cluster))) +
    geom_point() +
    geom_line() +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
    ylab("value")+
    labs(x = "condition", y = "value(i.e., Relative abundace)",
        title = "cluster centroids")



#using a cluster score to identify core genes

clust <- #Core genes for any cluster(e.g., cluster=3)
    hc_longdf[hc_longdf$cluster==3,]


hc <-#extract the selected cluster from the scaled data
    (scaledata[extClust==3,])

#compute correlation with cluster centroid
corScore <-
    function(x){cor(x,clust$value)}
cor.comp <-
    as.data.frame(apply(hc, 1, corScore))
colnames(cor.comp)[1] <-
    "cor.score"
cor.comp <-
    tibble::rownames_to_column(cor.comp, "gene")
#convert to long df
melt_hc <-
    melt(hc)
colnames(melt_hc) <-
    c('gene','condition','value')
melt_hc$condition <-
    factor(melt_hc$condition, levels=unique(melt_hc$condition))
#add the correlation score
melt_hc <-
    left_join(melt_hc, cor.comp, by='gene')

# plot
ggplot(melt_hc, aes(x=condition,y=value)) +
    geom_line(aes(colour=cor.score, group=gene)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
    scale_colour_gradientn(colours=c('black','yellow')) +
    #this adds the core
    geom_line(data=clust,
        aes(condition,value,
            group=cluster), color="red",inherit.aes=FALSE, size = 2) +
    labs(x = "condition", y = "value(i.e., Relative abundace)",
        title = "Cluster score to identify core genes")


# Add cluster number to the original df
df_Cluster <-
    df %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(cluster = hclustCutree)




##################################################################
#################### 2. k-means Clustering
##################################################################


#define the optimal number of clusters
fviz_nbclust(scaledata, kmeans, method = "wss")

set.seed(222)
kClust <-
    kmeans(scaledata, centers=4, nstart = 1000, iter.max = 20)



#Plotting the centroids
extClust <- #extract the cluster
    kClust$cluster



# function to find centroid in cluster i
ClustCent = function(i, dat, clusters) {
    ind = (clusters == i)
    colMeans(dat[ind,])
}
kClustcentroids <-
    sapply(levels(factor(extClust)), ClustCent, scaledata, extClust)

k_longdf <-
    melt(kClustcentroids)
colnames(k_longdf) <-
    c('condition','cluster','value')

#plot
k_longdf$condition <-
    factor(k_longdf$condition, levels=unique(k_longdf$condition))

print(ggplot(k_longdf,
    aes(x=condition,y=value,
        group=cluster,
        colour=as.factor(cluster))) +
        geom_point() +
        geom_line() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
        ylab("value")+
        labs(x = "condition", y = "value(i.e., Relative abundace)",
            title = "cluster centroids"))


#add cluster to df
df_Cluster <-
    df %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(cluster = kClust$cluster)

#using a cluster score to identify core genes

clust <- #Core genes for any cluster (e.g., cluster = 2)
    k_longdf[k_longdf$cluster==2,]


K <-###extract the selected cluster from the scaled data
    (scaledata[extClust==2,])

#compute correlation with cluster centroid
corScore <-
    function(x){cor(x,clust$value)}
cor.comp <-
    as.data.frame(apply(K, 1, corScore))
colnames(cor.comp)[1] <-
    "cor.score"
cor.comp <-
    tibble::rownames_to_column(cor.comp, "gene")
#convert to long df
meltK <-
    melt(K)
colnames(meltK) <-
    c('gene','condition','value')
meltK$condition <-
    factor(meltK$condition, levels=unique(meltK$condition))
#add the correlation score
meltK <-
    left_join(meltK, cor.comp, by='gene')

# plot
ggplot(meltK, aes(x=condition,y=value)) +
    geom_line(aes(colour=cor.score, group=gene)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
    scale_colour_gradientn(colours=c('black','yellow')) +
    #this adds the core
    geom_line(data=clust,
        aes(condition,value, group=cluster),
        color="red",inherit.aes=FALSE, size = 2) +
    labs(x = "condition", y = "value(i.e., Relative abundace)",
        title = "Cluster score to identify core genes")




##################################################################
#################### 3. FCM clustering
##################################################################

#add time point to the data
time <- c(1,2,3,4,5)
df.t <- rbind(time, df)
row.names(df.t)[1]<-"time"


#save it to a temp file
tmp <- tempfile()
write.table(df.t,file=tmp, sep='\t', quote = F, col.names=NA)

#read the file as an expression set
df.t <-
    table2eset(file=tmp)

#scale the data
scaleddf <-
    standardise(df.t)

#estimate the fuzzier
fuzz <-
    mestimate(scaleddf)

set.seed(222)
fuzzy_clust <- #cluster the data using FCM approach
    mfuzz(scaleddf,c=4,m=fuzz)


mfuzz.plot(scaleddf,cl=fuzzy_clust,mfrow=c(1,1),
    new.window=FALSE,
    time.labels=c(0,10,20,30,40))


#extract membership score
memScore <-
    acore(scaleddf,fuzzy_clust,min.acore=0)


memScore_df <-
    do.call(rbind, lapply(seq_along(memScore),
        function(i){ data.frame(CLUSTER=i, memScore[[i]])}))




##################################################################
#################### 4. SOM clustering
##################################################################

set.seed(222)
g <-
    somgrid(xdim = 3, #dimensions of the grid
        ydim = 3,#dimensions of the grid
        topo = "hexagonal")# topology of the grid.
map <-
    som(scaledata,
        grid = g,
        alpha = c(0.05,0.01),
        radius = 1)
plot(map,
    type = 'codes',
    palette.name = rainbow)
#how many data points are in each node
plot(map, type = "count")
#how many data points are in each node
plot(map, type = "mapping")

#add cluster number to the original df
df_Cluster <-
    cbind(df,map$unit.classif) %>%
    dplyr::rename(cluster = `map$unit.classif`)%>%
    tibble::rownames_to_column("Gene")





#Plotting the centroids
extClust <- #extract the cluster
    df_Cluster$cluster



# function to find centroid in cluster i
ClustCent = function(i, dat, clusters) {
    ind = (clusters == i)
    colMeans(dat[ind,])
}
kClustcentroids <-
    sapply(levels(factor(extClust)), ClustCent, scaledata, extClust)

k_longdf <-
    melt(kClustcentroids)
colnames(k_longdf) <-
    c('condition','cluster','value')

#plot
k_longdf$condition <-
    factor(k_longdf$condition, levels=unique(k_longdf$condition))

print(ggplot(k_longdf,
    aes(x=condition,y=value,
        group=cluster,
        colour=as.factor(cluster))) +
        geom_point() +
        geom_line() +
        theme_bw() +
        theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            axis.line = element_line(colour = "black")) +
        ylab("value")+
        labs(x = "condition", y = "value(i.e., Relative abundace)",
            title = "cluster centroids"))


#add cluster to df
df_Cluster <-
    df %>%
    as.data.frame(.) %>%
    tibble::rownames_to_column("Gene") %>%
    mutate(cluster = kClust$cluster)

#using a cluster score to identify core genes

clust <- #Core genes for any cluster (e.g., cluster = 2)
    k_longdf[k_longdf$cluster==2,]


K <-###extract the selected cluster from the scaled data
    (scaledata[extClust==2,])

#compute correlation with cluster centroid
corScore <-
    function(x){cor(x,clust$value)}
cor.comp <-
    as.data.frame(apply(K, 1, corScore))
colnames(cor.comp)[1] <-
    "cor.score"
cor.comp <-
    tibble::rownames_to_column(cor.comp, "gene")
#convert to long df
meltK <-
    melt(K)
colnames(meltK) <-
    c('gene','condition','value')
meltK$condition <-
    factor(meltK$condition, levels=unique(meltK$condition))
#add the correlation score
meltK <-
    left_join(meltK, cor.comp, by='gene')

# plot
ggplot(meltK, aes(x=condition,y=value)) +
    geom_line(aes(colour=cor.score, group=gene)) +
    theme_bw()+
    theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.line = element_line(colour = "black")) +
    scale_colour_gradientn(colours=c('black','yellow')) +
    #this adds the core
    geom_line(data=clust,
        aes(condition,value, group=cluster),
        color="red",inherit.aes=FALSE, size = 2) +
    labs(x = "condition", y = "value(i.e., Relative abundace)",
        title = "Cluster score to identify core genes")



