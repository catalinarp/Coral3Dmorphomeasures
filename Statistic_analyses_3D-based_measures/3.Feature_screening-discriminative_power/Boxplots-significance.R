###------------------------------------------------------------------------------
### Box Plots and t-tests
###------------------------------------------------------------------------------

### Aim: to show significant interspecific differences using the preliminary selected VARS

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load required libraries/packages
library(readr)
library(tidyverse)
library(rstatix)
library(ggpubr)
library(RColorBrewer)

## B) Set working directory: 
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

cols <- c("#114260","#4D724D","#D19240") # set species color scheme

## D) Prepare data frame:
PrelimselectedvarsTR <- prelimselVARSdf %>%
  mutate(Species = factor(Species,
                                levels = c("bifurcata", "cytherea", "hyacinthus"),
                                labels = c("A","B","C"))) %>%
  dplyr::select(-ID,
                -`Year`)

PrelimvarsTR_tibb <-  as_tibble(PrelimselectedvarsTR)
PrelimvarsTR_tibb %>% sample_n(6) # Check

mydata.long <- PrelimvarsTR_tibb %>%
  tidyr::pivot_longer(-Species,
                      names_to = "variables",values_to = "value")

mydata.long %>% sample_n(6) # Check

mydata.long %>%
  group_by(variables, Species) %>%
  dplyr::summarise(
    n = n(),
    mean = mean(value),
    sd = sd(value)
  ) %>%
  ungroup()

##-------------------------------------------------------------------------------
## 2. Calculate significant pairwise differences using the Bonferroni adjusted t-test
##-------------------------------------------------------------------------------

## A) t-test
stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ Species, p.adjust.method = "bonferroni")

## B) Remove unnecessary columns and display the outputs
ttest <- stat.test %>%
  dplyr::select(-.y., -statistic, -df) 


##-------------------------------------------------------------------------------
## 3. Create multi-panel Boxplots with t-test p-values
##-------------------------------------------------------------------------------

## A) Create the plot
myplot <- ggboxplot(
  mydata.long, x = "Species", y = "value",
  error.plot = "crossbar", add = "jitter", notch=FALSE,
  fill = "Species", palette = cols, legend = "none",
  ggtheme = theme_pubr(border = TRUE)) +
  facet_wrap(~variables, scales="free_y") +
  scale_fill_manual(name="Species",
                    values = c("#114260","#4D724D","#D19240"),
                    labels = c("A. cf. bifurcata (A)",
                               "A. cf. cytherea (B)",
                               "A. aff. hyacinthus (C)"))

## B) Add statistical test p-values
stat.test <- stat.test %>% add_xy_position(x = "Species")
myplot + stat_pvalue_manual(stat.test, label = "p.adj.signif")

## C) Save output
ggsave("CompleteDataSet/BoxplotsSigPrelimselvarsTR.pdf",
       width = 40, height = 40, units = "cm")
