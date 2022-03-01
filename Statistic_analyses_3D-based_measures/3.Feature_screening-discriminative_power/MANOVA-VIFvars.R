###------------------------------------------------------------------------------
### Multivariate analysis of variance - MANOVA
###------------------------------------------------------------------------------

### Aim: to test for significant multivariate interspecific differences

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load required libraries/packages
library(readr)
library(rstatix)

## B) Set working directory: 
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Upload data set:
# Variables that do not present multicollinearity:
VARSVIF_data <- read_delim("CompleteDataSet/allVARSVIF_data.txt", delim="\t",
                        col_types = cols(year = col_factor(levels = c("2015", "2018", "2019")),
                                         species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                         dataset = col_factor(levels = c("Testing", "Training"))))

summary(VARSVIF_data) #Check the data set characteristics

##-------------------------------------------------------------------------------
## 2. MANOVA
##-------------------------------------------------------------------------------

## A) Check for outliers (MANOVA is sensitive to them):

mahalanobis_distance(data = VARSVIF_data[,4:24])$is.outlier

## B) Prepare data set:

VARSVIF <- as.matrix(cbind(VARSVIF_data[,4:24]))

## C) Prepare data set:
MANOVAtest <- manova(VARSVIF ~ species,
                     data = VARSVIF_data)

## D) Inspect results:
summary(MANOVAtest, intercept=TRUE)
summary(MANOVAtest, intercept=TRUE,test="Wilks")
summary(MANOVAtest, intercept=TRUE,test="Roy")

