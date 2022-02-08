###------------------------------------------------------------------------------
### Morphodiagrams: Bi-variate scatter plots and density plots 
##-------------------------------------------------------------------------------

### Aim: to assess graphically the morphospaces overlapping between species

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------
## A) Load required libraries/packages

library(readr)
library(reshape2)
library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(RColorBrewer)
library(patchwork)

## B) Set working directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Load data :
ALLvarsTR <- as.data.frame(read_delim("CompleteDataSet/ALLvarsTR.txt",
                                      "\t", escape_double = FALSE,
                                      na = "NA", trim_ws = TRUE))

##-------------------------------------------------------------------------------
## 2. Select variables with significant differences between all variables
##-------------------------------------------------------------------------------
##
an <- lapply(ALLvarsTR[,7:59], function(x) aov(x~Species, data = ALLvarsTR))
Posthoc <- sapply(an, TukeyHSD , simplify=FALSE)
pval <- sapply(Posthoc, function(x) tail(unlist(x[1]), 3))
comparison <- as.factor(c("cytherea-bifurcata",
                          "hyacinthus-bifurcata",
                          "hyacinthus-cytherea"))
row.names(pval)<- comparison
pvaluedf.m <- melt(pval, value.name="p.value")


varsInterDif <- subset(pvaluedf.m, p.value < 0.05)

varsInterDif3spp <- varsInterDif %>%
  group_by(Var2) %>%
  filter(n() > 2) %>%
  droplevels %>%
  arrange(Var2) %>%
  slice(1) %>% 
  ungroup() %>% 
  dplyr::select(Var2)

vars_list <- as.factor(varsInterDif3spp$Var2)

varsInterDif_df <- ALLvarsTR %>%
  dplyr::select(Species,all_of(vars_list))


##-------------------------------------------------------------------------------
## 3. Plot the diagrams
##-------------------------------------------------------------------------------
## A) Assess colors to species

species <- factor(varsInterDif_df$Species,
                  levels = c("bifurcata","cytherea","hyacinthus"),
                  labels=c("A","B","C"))
cols <- c("#114260","#4D724D","#D19240")

## B) Check correlation and bi-plots

library(psych)

pairs.panels(varsInterDif_df[,2:11],
             gap = 0,
             hist.col="grey",
             bg = cols,
             pch = 21,
             ellipses=FALSE,
             density = TRUE,
             cex.cor=1,
             cex=1,
             stars = TRUE,
             smooth=FALSE)

## C) Plot bi-variate scatter plots and density plots

# Example: H_tip_var vs. H_var
scatter1 <- ggplot(varsInterDif_df,
                   aes(x=`H_tip_var$x.t`,
                       y=`H_var$x.t`,
                       color=species, fill=species)) +
  geom_point(size=2) + stat_ellipse(geom = "polygon",
                                    linetype = "blank",
                                    alpha=0.1) +
  scale_colour_manual(name = "Species",
                      labels = c("A","B","C"),
                      values = cols) +
  scale_fill_manual(name = "Species",
                    labels = c("A","B","C"),
                    values = cols) +
  theme_bw() +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank()) +
  geom_hline(yintercept = 0, color = "grey", linetype = "dashed",
             size=0.8)+
  geom_vline(xintercept = 0.3, color = "grey", linetype = "dashed",
             size=0.8)+
  scale_x_continuous(limits = c(-3, 3))+
  scale_y_continuous(limits = c(-3, 3))+
  theme(legend.position = c(.9, .2))+
  labs(y= "Mean curvature variance (H_var)",
       x = "Mean curvature variance at the tip (H_tip_var)")

scatter1

xdensity1 <- ggplot(varsInterDif_df,
                    aes(`H_tip_var$x.t`, fill=species)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(-3, 3)) + theme_classic() +
  theme(legend.position = "none") +
  labs(y= "Density",
       x = "H_tip_var") + 
  theme(axis.title.x=element_blank())

xdensity1

ydensity1 <- ggplot(varsInterDif_df,
                    aes(`H_var$x.t`, fill=species)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = cols) + theme_classic() +
  theme(legend.position = "none") + 
  scale_x_continuous(limits = c(-3, 3)) +
  labs(y= "Density",
       x = "H_var") +
  theme(axis.title.y=element_blank()) +
  coord_flip()

ydensity1

# Complete plot:
#(xdensity1 + plot_spacer())/ (scatter1 | ydensity1) 

xdensity1 + plot_spacer() + scatter1 + ydensity1 +
  plot_layout(ncol = 2, nrow = 2, 
                heights = c(2, 7),
                widths = c(6, 1))

## D) Density plots

varsInterDif_df %>%
 gather(key = "measure", value = "value", -Species) %>%
  ggplot(aes(value, fill=Species, color=Species)) +
  scale_colour_manual(name = "Species",
                      labels = c("A","B","C"),
                      values = cols) +
  scale_fill_manual(name = "Species",
                    labels = c("A","B","C"),
                    values = cols) +
  labs(y= "Density", x = "3D-based measure") + theme_bw() + 
  theme(legend.position = c(.9, .1))+
  facet_wrap(~ measure, scales = "free") +
  geom_density(alpha=.5) +
  coord_cartesian(xlim = c(-2.3, 2.3), ylim = c(0, 0.65)) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

