###------------------------------------------------------------------------------
### Linear discriminant analiysis (LDA)
###------------------------------------------------------------------------------

### Aim: to assess the discriminating power of variables that do not present multicollinearity

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load packages
library(readr)      
library(dplyr)
library(plyr)
library(reshape2)
library(MASS)
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggord)

## B) Set working directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses")

## C) Load the data sets
# Variables that do not present multicollinearity:
VARSVIF_data <- read_delim("CompleteDataSet/allVARSVIF_data.txt", delim="\t",
                           col_types = cols(year = col_factor(levels = c("2015", "2018", "2019")),
                                            species = col_factor(levels = c("bifurcata", "cytherea", "hyacinthus")),
                                            dataset = col_factor(levels = c("Testing", "Training"))))

summary(VARSVIF_data) # Check the data set characteristics

##-------------------------------------------------------------------------------
## 2. Perform linear discriminant analysis (LDA)
##-------------------------------------------------------------------------------

## A) LDA
ldasppMASS <- lda(x=VARSVIF_data[,c(4:24)],
                  grouping=VARSVIF_data$species,
                  method= "mle")

ldasppMASS

## B) Exploratory plots 
plot(ldasppMASS) 

ggord(ldasppMASS, VARSVIF_data$species,
      ylim = c(-5, 5),
      xlim = c(-6, 5))

##-------------------------------------------------------------------------------
## 3. Graphical representation
##-------------------------------------------------------------------------------

## A) Prepare data
lda.data <- cbind(VARSVIF_data, predict(ldasppMASS)$x)
df <- na.omit(lda.data)

cols <- c("#114260","#4D724D","#D19240")

## B) Biplot
a <- ggplot(df,
            aes(LD1, LD2, colour=species,
                fill = species)) +
  stat_ellipse(geom = "polygon",
               linetype = "blank",
               alpha=0.1,
               level = 0.95) +
  geom_point() + xlab(label = "LD1 (69.58%)") + 
  ylab(label = "LD2 (30.42%)") +
  scale_colour_manual(name = "Species",
                      labels = c("A", "B","C"),
                      values = cols) +
  scale_fill_manual(name = "Species",
                    labels = c("A", "B","C" ), 
                    values = cols) +
  theme_bw() + theme(legend.position     = "bottom") +
  theme(legend.text = element_text(face = "italic",size = 12)) +
  theme(legend.title = element_text(size = 14))
a

## C) Density plots
xdensity <- ggplot(df, aes(LD1, fill=species)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = cols) +
  scale_x_continuous(limits = c(-6, 5)) +
  theme(legend.position = "none") + theme_classic()
xdensity

ydensity <- ggplot(df, aes(LD2, fill=species)) + 
  geom_density(alpha=.5) + 
  scale_fill_manual(values = cols) + 
  scale_x_continuous(limits = c(-4, 5)) +
  theme(legend.position = "none") + theme_classic() +
  coord_flip()
ydensity


## C) LDA scaling information: Bar plots using coefficients of the linear discriminants

LDAscaling <- as.data.frame(ldasppMASS$scaling)

#LD1
ggplot(data=LDAscaling,
       aes(x=reorder(row.names(LDAscaling),desc(LD1)), y=(LD1), fill = LD1)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "#eae9e9", high = "#353434") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "LD1 coefficients")

#LD2
ggplot(data=LDAscaling,
       aes(x=reorder(row.names(LDAscaling),desc(LD2)), y=LD2, fill = LD2)) +
  geom_bar(stat="identity") +
  scale_fill_gradient(low = "#eae9e9", high = "#353434") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank()) +
  labs(y = "LD2 coefficients") 


##-------------------------------------------------------------------------------
## 4. Prediction accuracy tables
##-------------------------------------------------------------------------------

## A) Complete LDA model

predictALL <- predict(ldasppMASS,
                      VARSVIF_data[,c(4:24)])$class
tableALL <- table(Predicted = predictALL,
                Actual = VARSVIF_data$species)
tableALL # Confusion matrix 
sum(diag(tableALL))/sum(tableALL) # Prediction accuracy

## B) Validate model

# Create random subsets:
set.seed(123)
ind <- sample(2, nrow(VARSVIF_data),
              replace = TRUE,
              prob = c(0.7, 0.3))
training <- VARSVIF_data[ind==1,]
testing <- VARSVIF_data[ind==2,]

# LDA:
linear <- lda(species~., training) 

# Confusion matrix - training subset:
p1 <- predict(linear, training)$class
tab1 <- table(Predicted = p1, Actual = training$species)
tab1
sum(diag(tab1))/sum(tab1)

# Confusion matrix - testing subset:
p2 <- predict(linear, testing)$class
tab2 <- table(Predicted = p2, Actual = testing$species)
tab2
sum(diag(tab2))/sum(tab2)
