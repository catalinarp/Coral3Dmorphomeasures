###------------------------------------------------------------------------------
### Select variables according to their variance inflation factor (VIF)
###------------------------------------------------------------------------------

### Aim: to subset variables the do not present multicollinearity

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load packages
library(readr)      
library(dplyr)
library(tidyverse)
library(car)

## B) Set working directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")


## C) Load and prepare the data
ALLvarsTR <- read_delim("CompleteDataSet/ALLvarsTR.txt", delim="\t",
                        col_types = cols(`Collection_year`= col_factor(levels = c("2015", "2018", "2019")),
                                         Species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                         Dataset = col_factor(levels = c("Testing", "Training"))))

summary(ALLvarsTR) #Check the data set characteristics

ALLvarsTRmod <- clean_names(ALLvarsTR)
ALLvarsTRdf <-as.data.frame(dplyr::select(ALLvarsTRmod,
                            -(k_tip_mean_x_t),
                            -(k_tip_skew_x_t),
                            -(k2_tip_skew),
                            -(k2_tip_mean_x_t)))
summary(ALLvarsTRdf)

##-------------------------------------------------------------------------------
## 2. Variance inflation factor (VIF)
##-------------------------------------------------------------------------------

## A) Check multicollinearity in the dataset:

VIFmulticoltest <- vif(ALLvarsTRdf[,7:55]) 

# Save output
write.table(VIFmulticoltest, "CompleteDataSet/Tests/VIF-allVARS.txt",
            sep="\t", row.names = F) # To export results 

## B) Identify collinear variables that should be excluded:
multicoltest <- vifcor(ALLvarsTRdf[,7:55], th=0.9)
multicoltest@results

## C) Reduce the dataset accordingly:
multicolexcl <- vifstep(ALLvarsTRdf[,7:55], th=10) 
multicolexcl@results

# Save output
write.table(multicolexcl@results, "CompleteDataSet/Tests/VIF-allVARSincluded.txt",
            sep="\t", row.names = F) # To export results 

## D) Subset a list of the variables with VIF <10:

VARSVIF_list <- as.factor(multicolexcl@results[[1]])
VARSVIF_data <- dplyr::select(ALLvarsTRdf,
                              all_of(VARSVIF_list))
VARSVIF_data <- cbind(species=ALLvarsTRdf$species,
                      dataset=ALLvarsTRdf$dataset,
                      year=ALLvarsTRdf$collection_year,
                      VARSVIF_data)
# Save output
write.table(VARSVIF_data, "CompleteDataSet/allVARSVIF_data.txt",
sep="\t", row.names = F) # To export results 
