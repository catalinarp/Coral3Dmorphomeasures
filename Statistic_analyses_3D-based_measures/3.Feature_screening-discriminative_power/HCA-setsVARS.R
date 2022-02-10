###-------------------------------------------------------------------
### Hierarchical Clustering Analysis (HCA)
###-------------------------------------------------------------------

### Aim: to inspect clustering using the preliminary selected subset of variables 

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load packages
library(readr)      # load data 
library(dplyr)      # transform data if missing values
library(cluster)    # dissimilarity matrix
library(purrr)      # apply function to several items
library(nomclust)   # clustering algorithms
library(mclust)     # evaluating optimal number of clusters
library(NbClust)
library(FactoMineR)
library(factoextra)
library(cowplot)
library(pvclust)
library(reshape2)
library(tidyverse)
library(dendextend)
library(arsenal)    # to compare datasets 

## B) Set working directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Load the data sets

# All variables:
ALLvarsTR <- as.data.frame(read_delim("CompleteDataSet/ALLvarsTR.txt",
                                      "\t", escape_double = FALSE,
                                      na = "NA", trim_ws = TRUE))

summary(ALLvarsTR) #Check the data set characteristics

# Preliminary selected variables from list:
prelimselVARS <- as.data.frame(read_delim("CompleteDataSet/PrelselectecVARS.txt",
                                          delim="\t",col_names = F))

prelimselVARSlist <- as.factor(prelimselVARS$X1)
prelimselVARSdata <- dplyr::select(ALLvarsTR,
                                   all_of(prelimselVARSlist))
prelimselVARSdf <-cbind(ID=ALLvarsTR$ID,
                        Year=ALLvarsTR$Collection_year,
                        Species=ALLvarsTR$Species,
                        prelimselVARSdata)

##-------------------------------------------------------------------------------
## 2. Scale data and compute dissimilarity matrix
##-------------------------------------------------------------------------------

## A) Scale data to make the variables comparable and compute Euclidean distance matrix

# All variables
scaledata1 <- cbind(ALLvarsTR[ ,-c(7:59)],
                    scale(ALLvarsTR[,7:59], center=TRUE)) 
summary(scaledata1)

row.names(scaledata1) <- paste0(scaledata1$ID,sep="_",scaledata1$Species)

disimatrix1 <- dist(scaledata1[,7:59], method = "euclidean")

# Preliminary selected variables
scaledataPrelim <- cbind(prelimselVARSdf[ ,-c(4:32)],
                    scale(prelimselVARSdf[,4:32], center=TRUE)) 
summary(scaledataPrelim)

row.names(scaledataPrelim) <- paste0(scaledataPrelim$ID,sep="_",scaledataPrelim$Species)

disimatrix2 <- dist(scaledataPrelim[,4:32],
                        method = "euclidean")


## B) Preliminary distance heatmaps

plot1 <- fviz_dist(disimatrix1, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
plot1 #All variables


plot2 <- fviz_dist(disimatrix2, gradient = list(low = "#00AFBB", mid = "white", high = "#FC4E07"))
plot2 #Preliminary selected variables

##-------------------------------------------------------------------------------
## 3. Compare clustering methods using the agglomerative coefficient
##-------------------------------------------------------------------------------

## A)  List clustering methods to be compared
m <- c( "average", "single", "complete", "ward")
names(m) <- c( "average", "single", "complete", "ward")

## B) Create a function for comparison
# All variables
ac1 <- function(x) {
  agnes(disimatrix1, method = x)$ac # Change according to the dissimilarity matrix used
}

# Preliminary selected variables
ac2 <- function(x) {
  agnes(disimatrix2, method = x)$ac # Change according to the dissimilarity matrix used
}

## C) Check the result of the agglomerative coefficient (values closer to 1 suggest stronger clustering)
map_dbl(m, ac1) # All variables
map_dbl(m, ac2) # Preliminary selected variables

##-------------------------------------------------------------------------------
## 4. Generate the dendrograms with different clustering methods 
##-------------------------------------------------------------------------------
## NOTE: Compare those ones that gave best agglomerative coefficient values

## A) Ward:
# All variables
dev.off() # to clean plot window
hc1.1 <- agnes(disimatrix1, method = "ward") # Calculates clustering
pltree(hc1, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix1))

# Preliminary selected variables
dev.off() # to clean plot window
hc1.2 <- agnes(disimatrix2, method = "ward") # Calculates clustering
pltree(hc1.2, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix2))

## B) Average:
# All variables
hc2.1 <- agnes(disimatrix1, method = "average") 
pltree(hc2.1, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix1))

# Preliminary selected variables
hc2.2 <- agnes(disimatrix2, method = "average")
pltree(hc2.2, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix2))

## C) Complete:
# All variables
hc3.1 <- agnes(disimatrix1, method = "complete") 
pltree(hc3.1, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix1))

# Preliminary selected variables
hc3.2 <- agnes(disimatrix2, method = "complete")
pltree(hc3.2, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix2))

##-------------------------------------------------------------------------------
## 5. Determining number of optimal clusters
##-------------------------------------------------------------------------------

## A) Elbow method: function to compute total within-cluster sum of square 
wss <- function(k) {
  kmeans(disimatrix2, k, nstart = 50, iter.max = 20)$tot.withinss
} # Change the dissimilarity matrix each time

k.values <- 1:15
wss_values <- map_dbl(k.values, wss)

plot(k.values, wss_values,
     type="b", pch = 19, frame = FALSE, 
     xlab="Number of clusters K",
     ylab="Total within-clusters sum of squares") #Check inflection point

## B) Average Silhouette Method
## function to compute average silhouette for k clusters
avg_sil <- function(k) {
  km.res <- kmeans(disimatrix2, centers = k, nstart = 25)
  ss <- silhouette(km.res$cluster, dist(disimatrix1))
  mean(ss[, 3]) # Change the dissimilarity matrix each time
}
## Compute and plot wss for k = 2 to k = 10
k.values <- 2:10
## extract avg silhouette for 2-15 clusters
avg_sil_values <- map_dbl(k.values, avg_sil)

plot(k.values, avg_sil_values,
     type = "b", pch = 19, frame = FALSE, 
     xlab = "Number of clusters K",
     ylab = "Average Silhouettes") # Check for the peak

## C) Compute several methods at once
# All variables
nbclust_out <- NbClust(
  data = scaledata1[,7:59],
  distance = "euclidean",
  min.nc = 2, # minimum number of clusters
  max.nc = 7, # maximum number of clusters
  method = "ward.D" # one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans"
)

# Preliminary selected variables
nbclust_out <- NbClust(
  data = scaledataPrelim[,4:32],
  distance = "euclidean",
  min.nc = 2, # minimum number of clusters
  max.nc = 7, # maximum number of clusters
  method = "ward.D" # one of: "ward.D", "ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid", "kmeans"
)

##-------------------------------------------------------------------------------
## 6. Hierarchical Clustering with Bootstrapped p-values
##-------------------------------------------------------------------------------

# All variables
mydt <- as.data.frame(t(as.matrix(scaledata1[,7:59])))

fitaverage <- pvclust(mydt, method.hclust="average",
               method.dist="euclidean")
plot(fitaverage) # dendrogram with p values

fitWard <- pvclust(mydt, method.hclust="ward.D2",
                      method.dist="euclidean")
plot(fitWard) # dendrogram with p values

fitcomplete <- pvclust(mydt, method.hclust="complete",
                   method.dist="euclidean")
plot(fitcomplete) # dendrogram with p values

# Preliminary selected variables
mydt2 <- as.data.frame(t(as.matrix(scaledataPrelim[,4:32])))

fitaverage2 <- pvclust(mydt2, method.hclust="average",
                      method.dist="euclidean")
plot(fitaverage2) # dendrogram with p values

fitWard2 <- pvclust(mydt2, method.hclust="ward.D2",
                   method.dist="euclidean")
plot(fitWard2) # dendrogram with p values

fitcomplete2 <- pvclust(mydt2, method.hclust="complete",
                       method.dist="euclidean")
plot(fitcomplete2) # dendrogram with p values


# OPTIONAL: add rectangles around groups highly supported by the data
pvrect(fitaverage, alpha=.90) 

##-------------------------------------------------------------------------------
## 7. Compare two HCA directly
##-------------------------------------------------------------------------------

## A) Calculate dendrograms using the best clustering method
# All variables
dev.off()
hcA <- agnes(disimatrix1, method = "ward") 
  pltree(hcA, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix1))

# Preliminary selected variables
dev.off()
hcB <- agnes(disimatrix2, method = "ward") 
pltree(hcB, cex = 0.6, hang = -1, main = "Dendrogram", labels = rownames(disimatrix2))


## B) Plot a tanglegram to compare groups obtained 

dend1 <- as.dendrogram(hcA)
dend2 <- as.dendrogram(hcB)
tanglegram1 <- tanglegram(dend1, dend2)
dend_list <- dendlist(dend1, dend2)

scaledata1mod <- scaledata1 %>% 
  mutate(Species = case_when(Species == "bifurcata" ~ "#114260",
                             Species == "cytherea"  ~ "#4D724D",
                             Species == "hyacinthus"  ~ "#D19240"))

scaledataPrelimod <- scaledataPrelim %>% 
  mutate(Species = case_when(Species == "bifurcata" ~ "#114260",
                             Species == "cytherea"  ~ "#4D724D",
                             Species == "hyacinthus"  ~ "#D19240"))

  
# let's add some color:
colors_to_use1 <- scaledata1mod$Species
colors_to_use2 <- scaledataPrelimod$Species
# But sort them based on their order in dend:
colors_to_use1 <- colors_to_use1[order.dendrogram(dend1)]
colors_to_use2 <- colors_to_use2[order.dendrogram(dend2)]
# Now we can use them
labels_colors(dend1) <- colors_to_use1
labels_colors(dend2) <- colors_to_use2
# Now each species has a color
labels_colors(dend1) 
labels_colors(dend2) 


## C)  Align and plot two dendrograms side by side
tanglegram2 <- dendlist(dend1,dend2) %>%
  untangle(method = "step2side") %>% # Find the best alignment layout
  set("labels_cex", 0.5) %>%
  set("branches_k_color", value = c("#a7a7a7", "#424242"), k = 2) %>%
tanglegram(common_subtrees_color_lines = FALSE, 
           highlight_distinct_edges  = FALSE)

# Just to check the entanglement degree
tanglegram3 <- tanglegram(dend1, dend2,
                          highlight_distinct_edges = FALSE, # Turn-off dashed lines
                          common_subtrees_color_lines = FALSE, # Turn-off line colors
                          common_subtrees_color_branches = TRUE, # Color common branches 
                          main = paste("entanglement =", round(entanglement(tanglegram2), 2)))

