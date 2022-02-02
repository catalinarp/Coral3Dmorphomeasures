###------------------------------------------------------------------------------
### Exploratory stats for 3D-based coral measures
###------------------------------------------------------------------------------

### Aim: to test assumptions of normality, homocedasticity and transform variables 

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------
## A) Load required libraries/packages

library(readr)
library(ggpubr)
library(stats)
library(car)
library(bestNormalize)
library(MVN)

## B) Set working directory: 
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")


## C) Upload data set:
ALLvars <- read_csv("Data/3Dmorpho-data.csv",
                    col_types = cols(`Collection year`= col_factor(levels = c("2015", "2018", "2019")),
                                     species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                     Dataset = col_factor(levels = c("Testing", "Training"))))

summary(ALLvars) # Check column types

##-------------------------------------------------------------------------------
## 2. Evaluate normality assumption in the variables
##-------------------------------------------------------------------------------

## A) Load extra packages to manipulate the data:
library(tidyr)
library(tidyselect)
library(tidyverse)
library(dplyr)
library(broom)

## B) Shapiro-Wilk test:
SWtest <- ALLvars %>% 
  gather(key = "variable_name", value = "value", 7:59) %>% 
  group_by(variable_name)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtest, "CompleteDataSet/Tests/SWtest.txt",
            sep="\t",row.names = F)

## C) Shapiro-Wilk test per species:
SWtest_spp <- ALLvars %>% 
  gather(key = "variable_name", value = "value", 7:59) %>% 
  group_by(variable_name, Species)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtest_spp,
            "CompleteDataSet/Tests/SWtest_spp.txt",
            sep="\t",row.names = F)

##-------------------------------------------------------------------------------
## 3. Transform (TR) variables that violate the normality assumption
##-------------------------------------------------------------------------------
## For each variable with p-value <0.05: 
# A) Plot q-q plot and test normality using Shapiro-Wilk 
# B) Run bestNormalize for x5 and choose the most frequently suggested option
# C) Perform the suggested transformation
# D) Plot q-q plot again
# E) Plot q-q plot
# F) Perform Shapiro-Wilk test to check normality again

#b_angle_var
ggqqplot(ALLvars$b_angle_var)
shapiro.test(ALLvars$b_angle_var)
bestNormalize(ALLvars$b_angle_var)
b_angle_var<-log_x(ALLvars$b_angle_var)
ggqqplot(b_angle_var$x.t)
shapiro.test(b_angle_var$x.t)

#br_length_var
ggqqplot(ALLvars$br_length_var)
shapiro.test(ALLvars$br_length_var)
bestNormalize(ALLvars$br_length_var)
br_length_var<-boxcox(ALLvars$br_length_var)
ggqqplot(br_length_var$x.t)
shapiro.test(br_length_var$x.t)

#br_rate_var
ggqqplot(ALLvars$br_rate_var)
shapiro.test(ALLvars$br_rate_var)
bestNormalize(ALLvars$br_rate_var)
br_rate_var<-log_x(ALLvars$br_rate_var)
ggqqplot(br_rate_var$x.t)
shapiro.test(br_rate_var$x.t)

#br_spacing_v1_mean
ggqqplot(ALLvars$br_spacing_v1_mean)
shapiro.test(ALLvars$br_spacing_v1_mean)
bestNormalize(ALLvars$br_spacing_v1_mean)
br_spacing_v1_mean<-orderNorm(ALLvars$br_spacing_v1_mean)
ggqqplot(br_spacing_v1_mean$x.t)
shapiro.test(br_spacing_v1_mean$x.t)

#br_spacing_v1_var
ggqqplot(ALLvars$br_spacing_v1_var)
shapiro.test(ALLvars$br_spacing_v1_var)
bestNormalize(ALLvars$br_spacing_v1_var)
br_spacing_v1_var<-log_x(ALLvars$br_spacing_v1_var)
ggqqplot(br_spacing_v1_var$x.t)
shapiro.test(br_spacing_v1_var$x.t)

#br_spacing_v2_var
ggqqplot(ALLvars$br_spacing_v2_var)
shapiro.test(ALLvars$br_spacing_v2_var)
bestNormalize(ALLvars$br_spacing_v2_var)
br_spacing_v2_var<-log_x(ALLvars$br_spacing_v2_var)
ggqqplot(br_spacing_v2_var$x.t)
shapiro.test(br_spacing_v2_var$x.t)

#d_avg_var
ggqqplot(ALLvars$d_avg_var)
shapiro.test(ALLvars$d_avg_var)
bestNormalize(ALLvars$d_avg_var)
d_avg_var<-boxcox(ALLvars$d_avg_var)
ggqqplot(d_avg_var$x.t)
shapiro.test(d_avg_var$x.t)

#da_avg_mean
ggqqplot(ALLvars$da_avg_mean)
shapiro.test(ALLvars$da_avg_mean)
bestNormalize(ALLvars$da_avg_mean)
da_avg_mean<-boxcox(ALLvars$da_avg_mean)
ggqqplot(da_avg_mean$x.t)
shapiro.test(da_avg_mean$x.t)

#da_mean
ggqqplot(ALLvars$da_mean)
shapiro.test(ALLvars$da_mean)
bestNormalize(ALLvars$da_mean)
da_mean<-boxcox(ALLvars$da_mean)
ggqqplot(da_mean$x.t)
shapiro.test(da_mean$x.t)

#da_var
ggqqplot(ALLvars$da_var)
shapiro.test(ALLvars$da_var)
bestNormalize(ALLvars$da_var)
da_var<-boxcox(ALLvars$da_var)
ggqqplot(da_var$x.t)
shapiro.test(da_var$x.t)

#db_mean
ggqqplot(ALLvars$db_mean)
shapiro.test(ALLvars$db_mean)
bestNormalize(ALLvars$db_mean)
db_mean<-orderNorm(ALLvars$db_mean)
ggqqplot(db_mean$x.t)
shapiro.test(db_mean$x.t)

#db_var
ggqqplot(ALLvars$db_var)
shapiro.test(ALLvars$db_var)
bestNormalize(ALLvars$db_var)
db_var<-orderNorm(ALLvars$db_var)
ggqqplot(db_var$x.t)
shapiro.test(db_var$x.t)

#dc_mean
ggqqplot(ALLvars$dc_mean)
shapiro.test(ALLvars$dc_mean)
bestNormalize(ALLvars$dc_mean)
dc_mean<-orderNorm(ALLvars$dc_mean)
ggqqplot(dc_mean$x.t)
shapiro.test(dc_mean$x.t)

#dc_var
ggqqplot(ALLvars$dc_var)
shapiro.test(ALLvars$dc_var)
bestNormalize(ALLvars$dc_var)
dc_var<-log_x(ALLvars$dc_var)
ggqqplot(dc_var$x.t)
shapiro.test(dc_var$x.t)

#H_tip_skew
ggqqplot(ALLvars$H_tip_skew)
shapiro.test(ALLvars$H_tip_skew)
bestNormalize(ALLvars$H_tip_skew)
H_tip_skew<-orderNorm(ALLvars$H_tip_skew)
ggqqplot(H_tip_skew$x.t)
shapiro.test(H_tip_skew$x.t)

#H_tip_var
ggqqplot(ALLvars$H_tip_var)
shapiro.test(ALLvars$H_tip_var)
bestNormalize(ALLvars$H_tip_var)
H_tip_var<-orderNorm(ALLvars$H_tip_var)
ggqqplot(H_tip_var$x.t)
shapiro.test(H_tip_var$x.t)

#H_var
ggqqplot(ALLvars$H_var)
shapiro.test(ALLvars$H_var)
bestNormalize(ALLvars$H_var)
H_var<-orderNorm(ALLvars$H_var)
ggqqplot(H_var$x.t)
shapiro.test(H_var$x.t)

#K_mean
ggqqplot(ALLvars$K_mean)
shapiro.test(ALLvars$K_mean)
bestNormalize(ALLvars$K_mean)
K_mean<-orderNorm(ALLvars$K_mean)
ggqqplot(K_mean$x.t)
shapiro.test(K_mean$x.t)

#K_skew
ggqqplot(ALLvars$K_skew)
shapiro.test(ALLvars$K_skew)
bestNormalize(ALLvars$K_skew)
K_skew<-orderNorm(ALLvars$K_skew)
ggqqplot(K_skew$x.t)
shapiro.test(K_skew$x.t)

#K_tip_var
ggqqplot(ALLvars$K_tip_var)
shapiro.test(ALLvars$K_tip_var)
bestNormalize(ALLvars$K_tip_var)
K_tip_var<-orderNorm(ALLvars$K_tip_var)
ggqqplot(K_tip_var$x.t)
shapiro.test(K_tip_var$x.t)

#K_var
ggqqplot(ALLvars$K_var)
shapiro.test(ALLvars$K_var)
bestNormalize(ALLvars$K_var)
K_var<-orderNorm(ALLvars$K_var)
ggqqplot(K_var$x.t)
shapiro.test(K_var$x.t)

#k1_kurt
ggqqplot(ALLvars$k1_kurt)
shapiro.test(ALLvars$k1_kurt)
bestNormalize(ALLvars$k1_mean)
k1_kurt<-orderNorm(ALLvars$k1_kurt)
ggqqplot(k1_kurt$x.t)
shapiro.test(k1_kurt$x.t)

#k1_mean
ggqqplot(ALLvars$k1_mean)
shapiro.test(ALLvars$k1_mean)
bestNormalize(ALLvars$k1_mean)
k1_mean<-orderNorm(ALLvars$k1_mean)
ggqqplot(k1_mean$x.t)
shapiro.test(k1_mean$x.t)

#k1_skew
ggqqplot(ALLvars$k1_skew)
shapiro.test(ALLvars$k1_skew)
bestNormalize(ALLvars$k1_skew)
k1_skew<-boxcox(ALLvars$k1_skew)
ggqqplot(k1_skew$x.t)
shapiro.test(k1_skew$x.t)

#k1_tip_kurt
ggqqplot(ALLvars$k1_tip_kurt)
shapiro.test(ALLvars$k1_tip_kurt)
bestNormalize(ALLvars$k1_tip_kurt)
k1_tip_kurt<-orderNorm(ALLvars$k1_tip_kurt)
ggqqplot(k1_tip_kurt$x.t)
shapiro.test(k1_tip_kurt$x.t)

#k1_tip_mean
ggqqplot(ALLvars$k1_tip_mean)
shapiro.test(ALLvars$k1_tip_mean)
bestNormalize(ALLvars$k1_tip_mean)
k1_tip_mean<-orderNorm(ALLvars$k1_tip_mean)
ggqqplot(k1_tip_mean$x.t)
shapiro.test(k1_tip_mean$x.t)

#k1_tip_var
ggqqplot(ALLvars$k1_tip_var)
shapiro.test(ALLvars$k1_tip_var)
bestNormalize(ALLvars$k1_tip_var)
k1_tip_var<-orderNorm(ALLvars$k1_tip_var)
ggqqplot(k1_tip_var$x.t)
shapiro.test(k1_tip_var$x.t)

#k1_var
ggqqplot(ALLvars$k1_var)
shapiro.test(ALLvars$k1_var)
bestNormalize(ALLvars$k1_var)
k1_var<-orderNorm(ALLvars$k1_var)
ggqqplot(k1_var$x.t)
shapiro.test(k1_var$x.t)

#k2_mean
ggqqplot(ALLvars$k2_mean)
shapiro.test(ALLvars$k2_mean)
bestNormalize(ALLvars$k2_mean)
k2_mean<-orderNorm(ALLvars$k2_mean)
ggqqplot(k2_mean$x.t)
shapiro.test(k2_mean$x.t)

#k2_skew
ggqqplot(ALLvars$k2_skew)
shapiro.test(ALLvars$k2_skew)
bestNormalize(ALLvars$k2_skew)
k2_skew<-orderNorm(ALLvars$k2_skew)
ggqqplot(k2_skew$x.t)
shapiro.test(k2_skew$x.t)

#k2_tip_kurt
ggqqplot(ALLvars$k2_tip_kurt)
shapiro.test(ALLvars$k2_tip_kurt)
bestNormalize(ALLvars$k2_tip_kurt)
k2_tip_kurt<-yeojohnson(ALLvars$k2_tip_kurt)
ggqqplot(k2_tip_kurt$x.t)
shapiro.test(k2_tip_kurt$x.t)

#k2_tip_var
ggqqplot(ALLvars$k2_tip_var)
shapiro.test(ALLvars$k2_tip_var)
bestNormalize(ALLvars$k2_tip_var)
k2_tip_var<-orderNorm(ALLvars$k2_tip_var)
ggqqplot(k2_tip_var$x.t)
shapiro.test(k2_tip_var$x.t)

#k2_var
ggqqplot(ALLvars$k2_var)
shapiro.test(ALLvars$k2_var)
bestNormalize(ALLvars$k2_var)
k2_var<-orderNorm(ALLvars$k2_var)
ggqqplot(k2_var$x.t)
shapiro.test(k2_var$x.t)

#S
ggqqplot(ALLvars$S)
shapiro.test(ALLvars$S)
bestNormalize(ALLvars$S)
S<-boxcox(ALLvars$S)
ggqqplot(S$x.t)
shapiro.test(S$x.t)

##-------------------------------------------------------------------------------
## 4. Transform marginally normal overall and not normal when examined among species 
##-------------------------------------------------------------------------------

#br_spacing_v2_mean
ggqqplot(ALLvars$br_spacing_v2_mean)
shapiro.test(ALLvars$br_spacing_v2_mean)
bestNormalize(ALLvars$br_spacing_v2_mean)
br_spacing_v2_mean<-orderNorm(ALLvars$br_spacing_v2_mean)
ggqqplot(br_spacing_v2_mean$x.t)
shapiro.test(br_spacing_v2_mean$x.t)

#FD
ggqqplot(ALLvars$FD)
shapiro.test(ALLvars$FD)
bestNormalize(ALLvars$FD)
FD<-orderNorm(ALLvars$FD)
ggqqplot(FD$x.t)
shapiro.test(FD$x.t)

#H_kurt
ggqqplot(ALLvars$H_kurt)
shapiro.test(ALLvars$H_kurt)
bestNormalize(ALLvars$H_kurt)
H_kurt<-yeojohnson(ALLvars$H_kurt)
ggqqplot(H_kurt$x.t)
shapiro.test(H_kurt$x.t)

#H_tip_mean
ggqqplot(ALLvars$H_tip_mean)
shapiro.test(ALLvars$H_tip_mean)
bestNormalize(ALLvars$H_tip_mean)
H_tip_mean<-orderNorm(ALLvars$H_tip_mean)
ggqqplot(H_tip_mean$x.t)
shapiro.test(H_tip_mean$x.t)

#K_tip_mean
ggqqplot(ALLvars$K_tip_mean)
shapiro.test(ALLvars$K_tip_mean)
bestNormalize(ALLvars$K_tip_mean)
K_tip_mean<-center_scale(ALLvars$K_tip_mean)
ggqqplot(K_tip_mean$x.t)
shapiro.test(K_tip_mean$x.t)

#K_tip_skew
ggqqplot(ALLvars$K_tip_skew)
shapiro.test(ALLvars$K_tip_skew)
bestNormalize(ALLvars$K_tip_skew)
K_tip_skew<-log_x(ALLvars$K_tip_skew)
ggqqplot(K_tip_skew$x.t)
shapiro.test(K_tip_skew$x.t)

#k2_kurt
ggqqplot(ALLvars$k2_kurt)
shapiro.test(ALLvars$k2_kurt)
bestNormalize(ALLvars$k2_kurt)
k2_kurt<-boxcox(ALLvars$k2_kurt)
ggqqplot(k2_kurt$x.t)
shapiro.test(k2_kurt$x.t)

#k2_tip_mean
ggqqplot(ALLvars$k2_tip_mean)
shapiro.test(ALLvars$k2_tip_mean)
bestNormalize(ALLvars$k2_tip_mean)
k2_tip_mean<-orderNorm(ALLvars$k2_tip_mean)
ggqqplot(k2_tip_mean$x.t)
shapiro.test(k2_tip_mean$x.t)

##-------------------------------------------------------------------------------
## 5. Create new data set with transformed variables (TR)
##-------------------------------------------------------------------------------

varsNOTnormal <- subset(SWtest, p.value < 0.05)
varsNOTn_list <- as.factor(varsNOTnormal$variable_name) #List of TR variables
#varsNOTn_data <- dplyr::select(ALLvars, all_of(varsNOTn_list))

ALLvarsTR <- dplyr::select(ALLvars,
                    -all_of(varsNOTn_list),
                    -(br_spacing_v2_mean),
                    -(H_kurt),
                    -(FD),
                    -(H_tip_mean),
                    -(K_tip_mean),
                    -(K_tip_skew),
                    -(k2_kurt),
                    -(k2_tip_mean))

ALLvarsTR <- cbind(ALLvarsTR,
                   b_angle_var$x.t, br_length_var$x.t,br_rate_var$x.t,
                   br_spacing_v1_mean$x.t, br_spacing_v1_var$x.t,
                   br_spacing_v2_mean$x.t,br_spacing_v2_var$x.t,d_avg_var$x.t,
                   da_avg_mean$x.t,da_mean$x.t,da_var$x.t,
                   db_mean$x.t,db_var$x.t,dc_mean$x.t,dc_var$x.t,FD$x.t,
                   H_tip_mean$x.t,H_tip_skew$x.t,H_tip_var$x.t,H_var$x.t,H_kurt$x.t,
                   K_mean$x.t,K_skew$x.t,K_tip_mean$x.t,K_tip_skew$x.t,
                   K_tip_var$x.t,K_var$x.t,k1_kurt$x.t,k1_mean$x.t,k1_skew$x.t,
                   k1_tip_kurt$x.t,k1_tip_mean$x.t,k1_tip_var$x.t,k1_var$x.t,
                   k2_kurt$x.t,k2_mean$x.t,k2_skew$x.t,k2_tip_kurt$x.t,
                   k2_tip_mean$x.t,k2_tip_var$x.t,k2_var$x.t,S$x.t)

# Save output
write.table(ALLvarsTR, "CompleteDataSet/ALLvarsTR.txt",
            sep="\t", row.names = FALSE) # To export results 

##-------------------------------------------------------------------------------
## 6. Evaluate normality assumption in the transformed (TR) variables
##-------------------------------------------------------------------------------

SWtestTR <- ALLvarsTR %>% 
  gather(key = "variable_name", value = "value", 7:59) %>% 
  group_by(variable_name)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtestTR, "CompleteDataSet/Tests/SWtestTR.txt",
            sep="\t", row.names = F)

SWtest_sppTR <- ALLvarsTR %>% 
  gather(key = "variable_name", value = "value", 7:59) %>% 
  group_by(variable_name, Species)  %>% 
  do(tidy(shapiro.test(.$value))) %>% 
  ungroup() %>% 
  dplyr::select(-method)

# Save output
write.table(SWtest_sppTR, "CompleteDataSet/Tests/SWtest_sppTR.txt",
            sep="\t",row.names = F)

##-------------------------------------------------------------------------------
## 7. Evaluate multivariate normality in the transformed (TR) variables
##-------------------------------------------------------------------------------

MVNtest <- mvn(ALLvarsTR[7:59], subset = NULL,
               mvnTest = c("hz"), # Change the test type each time to have all values:"mardia","hz","royston","dh","energy"
               covariance = TRUE, tol = 1e-30, alpha = 0.5,
               scale = FALSE, desc = TRUE, transform = "none", R = 1000,
               univariateTest = c("SW", "CVM", "Lillie", "SF", "AD"),
               univariatePlot = "none", multivariatePlot = "none",
               multivariateOutlierMethod = "none", bc = FALSE,
               bcType = "rounded",
               showOutliers = FALSE, showNewData = FALSE)

#To inspect the results of the MVN test:
MVNtest$multivariateNormality

##-------------------------------------------------------------------------------
## 8. Evaluate homocedasticity assumption in the transformed (TR) variables
##-------------------------------------------------------------------------------

## A) Load file in case it is not in the environment:

ALLvarsTR <- read_delim("CompleteDataSet/ALLvarsTR.txt", delim="\t",
                        col_types = cols(`Collection_year`= col_factor(levels = c("2015", "2018", "2019")),
                                         Species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                         Dataset = col_factor(levels = c("Testing", "Training"))))

summary(ALLvarsTR) # Check column types

## B) Perform Levene test of homogeneity of variance:

levene_test_resultsTR <-lapply(ALLvarsTR[,-c(1:6)], leveneTest, group = ALLvarsTR$Species)

# Save output
write.table(levene_test_resultsTR, "CompleteDataSet/Tests/Levenetest-TR.txt",
            sep="\t", row.names = F)

