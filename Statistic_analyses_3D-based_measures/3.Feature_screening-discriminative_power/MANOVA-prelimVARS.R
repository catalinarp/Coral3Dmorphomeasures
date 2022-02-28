###------------------------------------------------------------------------------
### Multivariate analysis of variance - MANOVA
###------------------------------------------------------------------------------

### Aim: to test for significant interspecific differences using the preliminary selected VARS

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load required libraries/packages
library(readr)
library(rstatix)

## B) Set working directory: 
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Upload data set:
# All variables:
ALLvarsTR <- read_delim("CompleteDataSet/ALLvarsTR.txt", delim="\t",
                        col_types = cols(`Collection_year`= col_factor(levels = c("2015", "2018", "2019")),
                                         Species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                         Dataset = col_factor(levels = c("Testing", "Training"))))

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
                        Dataset=ALLvarsTR$Dataset,
                        prelimselVARSdata) 


##-------------------------------------------------------------------------------
## 2. MANOVA
##-------------------------------------------------------------------------------

## A) Check for outliers (MANOVA is sensitive to them):

mahalanobis_distance(data = prelimselVARSdf[,4:32])$is.outlier

## B) Prepare data set:

prelimselVARS <- as.matrix(cbind(prelimselVARSdf[,4:32]))

## C) Prepare data set:
MANOVAtest <- manova(prelimselVARS ~ Species,
                     data = prelimselVARSdf)

## D) Inspect results:
summary(MANOVAtest, intercept=TRUE)
summary(MANOVAtest, intercept=TRUE,test="Wilks")
summary(MANOVAtest, intercept=TRUE,test="Roy")

