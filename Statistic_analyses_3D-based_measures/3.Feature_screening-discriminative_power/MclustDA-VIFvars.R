###------------------------------------------------------------------------------
### MclustDA - Discriminant analysis based on Gaussian finite mixture modelling
###------------------------------------------------------------------------------

### Aim: to assess the discriminating power of variables that do not present multicollinearity

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load packages
library(readr)
library(mclust)
library(dplyr)
library(tidyverse)
library(car)

## B) Set working directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")


## C) Load the data set
# All variables:
ALLvarsTR <- read_delim("CompleteDataSet/ALLvarsTR.txt", delim="\t",
                        col_types = cols(`Collection_year`= col_factor(levels = c("2015", "2018", "2019")),
                                         Species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                         Dataset = col_factor(levels = c("Testing", "Training"))))

summary(ALLvarsTR) #Check the data set characteristics

# Variables that do not present multicollinearity:
VIFvars <- as.data.frame(read_delim("CompleteDataSet/VIFvars.txt",
                                             delim="\t",col_names = F))
VIFvarslist <- as.factor(VIFvars$X1)
VIFvarsdata <- dplyr::select(ALLvarsTR,
                                      all_of(VIFvarslist))
VIFvarsdatadf <-cbind(ID=ALLvarsTR$ID,
                      Year=ALLvarsTR$Collection_year,
                      Species=ALLvarsTR$Species,
                      VIFvarsdata)

str(VIFvarsdatadf)


# Set factor for species 
species <- factor(VIFvarsdatadf$species,
                    levels= c("bifurcata",
                              "cytherea",
                              "hyacinthus"),
                    labels= c("A","B","C"))

X <- VIFvarsdatadf[,4:14] # Select only numeric variables

##-------------------------------------------------------------------------------
## 2. Discriminant analysis based on Gaussian finite mixture modelling
##-------------------------------------------------------------------------------

## A) Create training and testing dataset:
set.seed(36734)

train <- sample(1:nrow(X),
                size = round(nrow(X)*0.68),
                replace = FALSE)

X.train <- X[train,] 
Class <-VIFvarsdatadf[,3]
Class.train <- Class[train]
table(Class.train)

X.test <- X[-train,]
Class.test <- Class[-train]
table(Class.test)

## B) Discriminant analysis based on EDDA model 
# Common covariance structure selected by BIC: imposes a single mixture component for each group
modEDDA <- MclustDA(X.train, Class.train, modelType = "EDDA")
summary(modEDDA, newdata = X.test, newclass = Class.test)
#plot(modEDDA, what = "scatterplot")

## C) EDDA model cross-validation (10-fold):
cvEDDA <- cvMclustDA(modEDDA, nfold = 10)
unlist(cvEDDA[c("error", "se")])

## D) Discriminant analysis based on MclustDA model
# General covariance structure selected by BIC: a finite mixture of Gaussian distributions is used within each class
modMclustDA <- MclustDA(X.train, Class.train, modelType = "MclustDA")
summary(modMclustDA, newdata = X.test, newclass = Class.test)

## C) MclustDA model cross-validation (10-fold):
cvMclustDA <- cvMclustDA(modMclustDA)
unlist(cvMclustDA[c("error", "se")])

##-------------------------------------------------------------------------------
## 3. Exploratory plots
##-------------------------------------------------------------------------------

## A) Discriminant analysis based on EDDA model 
plot(modEDDA)
plot(modEDDA, dimens = 3:6) # Plot bivariate

## B) Discriminant analysis based on MclustDA model
plot(modMclustDA)
