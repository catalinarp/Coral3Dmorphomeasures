###------------------------------------------------------------------------------
### Correlation between preliminary selected variables
###------------------------------------------------------------------------------

### Aim: Evaluate correlation between variables

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------
## A) Load required libraries/packages
library(Hmisc)
library(readr)
library(psych)
library(RColorBrewer)

## B) Set working directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Upload data set:
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
## 2. Assess correlation using Pearson correlation coefficients
##-------------------------------------------------------------------------------

## A) Create a function: 
flattenCorrMatrix <- function(cormat, pmat) {
  ut <- upper.tri(cormat)
  data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
}

## B) Apply function to the dataset: 
corrALLvarsTR<-rcorr(as.matrix(prelimselVARSdf[,4:32]))

corrALLvarsTR.tab <-flattenCorrMatrix(corrALLvarsTR$r, corrALLvarsTR$P)

# Save output
write.table(corrALLvarsTR.tab, "CompleteDataSet/Tests/corrPrelimselVARSTR.tab.txt",
            sep="\t", row.names = F)

## C) Plot graphical output:

spp <- factor(prelimselVARSdf$Species,
              levels= c("bifurcata",
                        "cytherea",
                        "hyacinthus"),
              labels= c("A",
                        "B",
                        "C"))

pairs.panels(prelimselVARSdf[,4:32],
             gap = 0,
             hist.col="grey",
             bg = c("#114260","#4D724D",
                    "#D19240")[spp],
             pch = 21,
             ellipses=FALSE,
             cex.cor=4,
             cex=1,
             stars = TRUE)
