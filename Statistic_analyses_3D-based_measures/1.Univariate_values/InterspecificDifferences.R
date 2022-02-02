###-------------------------------------------------------------------------------
### Interspecific differences 
###-------------------------------------------------------------------------------

### Aim: Evaluate significant differences between species

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------
## A) Load required libraries/packages

library(readr)
library(stats)
library(rstatix)
library(reshape2)

## B) Set working directory: 
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Upload data set:
ALLvarsTRdf <- as.data.frame(read_delim("CompleteDataSet/ALLvarsTR.txt", delim="\t",
                                        col_types = cols(`Collection_year`= col_factor(levels = c("2015", "2018", "2019")),
                                                         Species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                                         Dataset = col_factor(levels = c("Testing", "Training")))))
##-------------------------------------------------------------------------------
## 2. Assessing interspecific differences
##-------------------------------------------------------------------------------

## A) Analysis of variance (ANOVA)

# Store all formulae in a list
formulae <- lapply(colnames(ALLvarsTRdf)[7:ncol(ALLvarsTRdf)], function(x) as.formula(paste0("`",x,"`", " ~ Species")))
ANOVA_resultsTR <- lapply(formulae, function(x) summary(aov(x, data = ALLvarsTRdf)))
names(ANOVA_resultsTR) <- format(formulae)

# Extract just p-values
p <- unlist(lapply(ANOVA_resultsTR, function(x) x[[1]]$"Pr(>F)"[1]))
dfpval <- cbind(variable=colnames(ALLvarsTRdf[,7:59]),p)

# Save output
write.table(dfpval, "CompleteDataSet/Tests/ANOVA_resultsTR.txt", sep="\t",
            row.names = FALSE)


# NOTE: List variables with p-values < 0.05 and create a subset
ANOVATRS_list <- as.list(read_delim("Data/ANOVA-TRS_list.txt",
                                    delim="\t",
                                    col_names = FALSE))


varsANOVA_list <- as.factor(ANOVATRS_list$X1)

varsANOVA_data <- dplyr::select(ALLvarsTR, varsANOVA_list)

write.table(varsANOVA_data, "Data/varsANOVA_data.txt", sep="\t",
            row.names = FALSE)


## B)  Post hoc Tukey test:

an <- lapply(ALLvarsTRdf[,7:59], function(x) aov(x~Species, data = ALLvarsTRdf))
Posthoc <- sapply(an, TukeyHSD , simplify=FALSE)
pval <- sapply(Posthoc, function(x) tail(unlist(x[1]), 3))
comparison <- as.factor(c("cytherea-bifurcata",
                          "hyacinthus-bifurcata",
                          "hyacinthus-cytherea"))
row.names(pval)<- comparison

pvaluedf.m <- melt(pval, value.name="p.value")

# Save output
write.table(pvaluedf.m, "CompleteDataSet/Tests/PostHocTR-Tukey.txt",
            sep="\t", row.names = F)

## C) Two-sample t-tests (Bonferroni adjusted)

ALLvarsTR <- ALLvarsTR %>%
  mutate(Species = factor(Species,
                          levels = c("bifurcata", "cytherea", "hyacinthus"),
                          labels = c("A","B","C"))) %>%
  dplyr::select(-ID,
                -`Collection_year`,
                -Genus,
                -`Open_nomenclature`,
                -`Dataset`) 

ALLvarsTR_tibb <-  as_tibble(ALLvarsTR)
ALLvarsTR_tibb %>% sample_n(6) # Check

mydata.long <- ALLvarsTR_tibb %>%
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

stat.test <- mydata.long %>%
  group_by(variables) %>%
  t_test(value ~ Species, p.adjust.method = "bonferroni")


# Remove unnecessary columns and display the outputs
ttest <- stat.test %>%
  dplyr::select(-.y., -statistic, -df) 

# Save output
write.table(ttest, "CompleteDataSet/Tests/t-test_Bonf-p-adj.txt",
            sep="\t", row.names = F)

## D) OPTIONAL: Post hoc test with letters
library(agricolae)

List <- names(ALLvarsTRdf)[7:59] # select just the variables

model1 <- lapply(List, function(x) {
  lm(substitute(i~Species, list(i = as.name(x))), data = ALLvarsTRdf)})

lapply(model1, summary)

letters = lapply(model1, function(m) HSD.test((m), "Species", group = TRUE, console = TRUE))

