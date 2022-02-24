###------------------------------------------------------------------------------
### Heatmaps: Distribution comparisons using the Mann-Whitney test results
###------------------------------------------------------------------------------

### Aim: To test for interspecific differences between the distributions of the variables

##-------------------------------------------------------------------------------
## 1. Preliminary steps
##-------------------------------------------------------------------------------

## A) Load packages
library(readr)
library(heatmap3)
library(data.table)

## B) Move to the files directory
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses/Data/manwhit_test_full/")

## C) Load and prepare the data
# Set sample IDs
IDs <- as.data.frame(read_csv("IDs.csv", col_names = FALSE))
IDsf <- as.factor(IDs$X1)

# Create a function to sum matrices
add_matrices <- function(...) {
  a <- (...)
  cols <- sort(unique(unlist(lapply(a, colnames))))
  rows <- sort(unique(unlist(lapply(a, rownames))))
  out <- array(0, dim = c(length(rows), length(cols)), dimnames = list(rows,cols))
  for (m in a) out[rownames(m), colnames(m)] <- out[rownames(m), colnames(m)] + m
  out
}

# Set colors for heatmap gradient
gradientcolors<- c("white","#ecebeb","#4d4d4d")


## D) Set working directory according to the p-value files location
setwd("C:/Users/catal/Dropbox/Research/1.ULB_PhD/Thesis/Chapter 3 - Morphometric species delimitation Acropora/Manuscript_MethodsEcolEvol2021/Analyses/Data/manwhit_test_full/mwpval/")

##-------------------------------------------------------------------------------
## 2. Distribution comparison score (DCS): Summarize replicates of the Mann-Whitney U test per each variable
##-------------------------------------------------------------------------------
# p-values obtained from each of the pairwise comparisons were transformed
# into integers according to an alpha (??) of 0.05 of significance: 
# if p-value > 0.05, then = 1 (i.e., samples come from similar distributions)
# if p-value ??? 0.05, then = -1 (i.e., samples do not come from similar distributions)
# The integer values of the ten replicates were then added cumulatively to obtain 
# a final value or distribution comparison score (DCS)

## For each variable-------------------------------------------------------------

## A) Create list of file names

## B) list_matrices: Set function to summarize all matrices per each variable
# read each file into a matrix, assign the rownames and colnames as with the
# single file case, put everything in a list.
# Now each matrix is in a component of the list

## C) Combine all matrices for the same variable in one

## D) Format the matrices

## E) Plot the heatmap of the sumamrized matrices for each variable

#--------------------------------------------------------------------------------
# 2.1 Gaussian curvature (K)
#--------------------------------------------------------------------------------

## A)
files1 = list.files(path=".", pattern="manwhit_curvatures_new_Gauss_size1000_")
names1 <- substr(files1,1,39) #Create list of file names without the ".csv" part 

## B)
list_matrices <-
  lapply(files1, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy value where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix1 <- add_matrices(list_matrices)
fullmatrix1 <- fullmatrix1[!rownames(fullmatrix1) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix1) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples1 <- as.factor(rownames(fullmatrix1))
IDcol1 <-as.data.frame(Samples1, colnames="X1")
colnames(IDs)<- c("Samples1", "Morphospecies")
IDcol1 <-merge(IDcol1,IDs, by  = "Samples1") 
Morphospecies1 <- as.factor(IDcol1$Morphospecies)

cols1 <- c("#114260","#4D724D","#D19240")[Morphospecies1]

## E)
pdf(file='Figures/MWtest-pvalbinary_K.pdf')
heatmap3(fullmatrix1, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="average", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols1, 
         ColSideLabs = "Morphospecies",
         labCol=Samples1,labRow=Samples1,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.2 Gaussian curvature at the tip of the branch (K_tip)
#--------------------------------------------------------------------------------

## A)
files2 = list.files(path=".", pattern="manwhit_curvatures_new_Gauss_tips_size1000_")
names2 <- substr(files2,1,44) #Create list of file names without the ".csv" part 

## B)
list_matrices2 <-
  lapply(files2, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix2 <- add_matrices(list_matrices2)
fullmatrix2 <- fullmatrix2[!rownames(fullmatrix2) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix2) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples2 <- as.factor(rownames(fullmatrix2))
IDcol2 <-as.data.frame(Samples2, colnames="X1")
colnames(IDs)<- c("Samples2", "Morphospecies")
IDcol2 <-merge(IDcol2,IDs, by  = "Samples2") 
Morphospecies2 <- as.factor(IDcol2$Morphospecies)

cols2 <- c("#114260","#4D724D","#D19240")[Morphospecies2]

## E)
pdf(file='Figures/MWtest-pvalbinary_K-tip.pdf')
heatmap3(fullmatrix2, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols2, 
         ColSideLabs = "Morphospecies",
         labCol=Samples2,labRow=Samples2,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.3 Maximum curvature (k1)
#--------------------------------------------------------------------------------

## A)
files3 = list.files(path=".", pattern="manwhit_curvatures_new_Maximum_size1000_")
names3 <- substr(files3,1,41) #Create list of file names without the ".csv" part 

## B)
list_matrices3 <-
  lapply(files3, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })


## C)
fullmatrix3 <- add_matrices(list_matrices3)
fullmatrix3 <- fullmatrix3[!rownames(fullmatrix3) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix3) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples3 <- as.factor(rownames(fullmatrix3))
IDcol3 <-as.data.frame(Samples3, colnames="X1")
colnames(IDs)<- c("Samples3", "Morphospecies")
IDcol3 <-merge(IDcol3,IDs, by  = "Samples3") 
Morphospecies3 <- as.factor(IDcol3$Morphospecies)

cols3 <- c("#114260","#4D724D","#D19240")[Morphospecies3]

## E)
pdf(file='Figures/MWtest-pvalbinary_k1.pdf')
heatmap3(fullmatrix3, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols3, 
         ColSideLabs = "Morphospecies",
         labCol=Samples3,labRow=Samples3,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.4 Maximum curvature at branch tip (k1_tip)
#--------------------------------------------------------------------------------

## A)
files4 = list.files(path=".", pattern="manwhit_curvatures_new_Maximum_tips_size1000_")
names4 <- substr(files4,1,46) #Create list of file names without the ".csv" part 

## B)
list_matrices4 <-
  lapply(files4, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix4 <- add_matrices(list_matrices4)
fullmatrix4 <- fullmatrix4[!rownames(fullmatrix4) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix4) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples4 <- as.factor(rownames(fullmatrix4))
IDcol4 <-as.data.frame(Samples4, colnames="X1")
colnames(IDs)<- c("Samples4", "Morphospecies")
IDcol4 <-merge(IDcol4,IDs, by  = "Samples4") 
Morphospecies4 <- as.factor(IDcol4$Morphospecies)

cols4 <- c("#114260","#4D724D","#D19240")[Morphospecies4]

## E)
pdf(file='Figures/MWtest-pvalbinary_k1-tip.pdf')
heatmap3(fullmatrix4, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols4, 
         ColSideLabs = "Morphospecies",
         labCol=Samples4,labRow=Samples4,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.5 Mean curvature (H)
#--------------------------------------------------------------------------------

## A)
files5 = list.files(path=".", pattern="manwhit_curvatures_new_Mean_size1000_")
names5 <- substr(files5,1,38) #Create list of file names without the ".csv" part 

## B)
list_matrices5 <-
  lapply(files5, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix5 <- add_matrices(list_matrices5)
fullmatrix5 <- fullmatrix5[!rownames(fullmatrix5) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix5) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples5 <- as.factor(rownames(fullmatrix5))
IDcol5 <-as.data.frame(Samples5, colnames="X1")
colnames(IDs)<- c("Samples5", "Morphospecies")
IDcol5 <-merge(IDcol5,IDs, by  = "Samples5") 
Morphospecies5 <- as.factor(IDcol5$Morphospecies)

cols5 <- c("#114260","#4D724D","#D19240")[Morphospecies5]

## E)
pdf(file='Figures/MWtest-pvalbinary_H.pdf')
heatmap3(fullmatrix5, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols5, 
         ColSideLabs = "Morphospecies",
         labCol=Samples5,labRow=Samples5,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.6 Mean curvature at the branch tips (H_tip)
#--------------------------------------------------------------------------------

## A)
files6 = list.files(path=".", pattern="manwhit_curvatures_new_Mean_tips_size1000_")
names6 <- substr(files6,1,43) #Create list of file names without the ".csv" part 

## B)
list_matrices6 <-
  lapply(files6, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix6 <- add_matrices(list_matrices6)
fullmatrix6 <- fullmatrix6[!rownames(fullmatrix6) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix6) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples6 <- as.factor(rownames(fullmatrix6))
IDcol6 <-as.data.frame(Samples6, colnames="X1")
colnames(IDs)<- c("Samples6", "Morphospecies")
IDcol6 <-merge(IDcol6,IDs, by  = "Samples6") 
Morphospecies6 <- as.factor(IDcol6$Morphospecies)

cols6 <- c("#114260","#4D724D","#D19240")[Morphospecies6]

## E)
pdf(file='Figures/MWtest-pvalbinary_H_tip.pdf')
heatmap3(fullmatrix6, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols6, 
         ColSideLabs = "Morphospecies",
         labCol=Samples6,labRow=Samples6,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.7 Minimum curvature at the branch tips (k2)
#--------------------------------------------------------------------------------

## A)
files7 = list.files(path=".", pattern="manwhit_curvatures_new_Minimum_size1000_")
names7 <- substr(files7,1,41) #Create list of file names without the ".csv" part 

## B)
list_matrices7 <-
  lapply(files7, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix7 <- add_matrices(list_matrices7)
fullmatrix7 <- fullmatrix7[!rownames(fullmatrix7) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix7) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples7 <- as.factor(rownames(fullmatrix7))
IDcol7 <-as.data.frame(Samples7, colnames="X1")
colnames(IDs)<- c("Samples7", "Morphospecies")
IDcol7 <-merge(IDcol7,IDs, by  = "Samples7") 
Morphospecies7 <- as.factor(IDcol7$Morphospecies)

cols7 <- c("#114260","#4D724D","#D19240")[Morphospecies7]

## E)
pdf(file='Figures/MWtest-pvalbinary_k2.pdf')
heatmap3(fullmatrix7, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols7, 
         ColSideLabs = "Morphospecies",
         labCol=Samples7,labRow=Samples7,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.8 Minimum curvature at the branch tips (k2_tip)
#--------------------------------------------------------------------------------

## A)
files8 = list.files(path=".", pattern="manwhit_curvatures_new_Minimum_tips_size1000_")
names8 <- substr(files8,1,46) #Create list of file names without the ".csv" part 

## B)
list_matrices8 <-
  lapply(files8, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix8 <- add_matrices(list_matrices8)
fullmatrix8 <- fullmatrix8[!rownames(fullmatrix8) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix8) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples8 <- as.factor(rownames(fullmatrix8))
IDcol8 <-as.data.frame(Samples8, colnames="X1")
colnames(IDs)<- c("Samples8", "Morphospecies")
IDcol8 <-merge(IDcol8,IDs, by  = "Samples8") 
Morphospecies8 <- as.factor(IDcol8$Morphospecies)

cols8 <- c("#114260","#4D724D","#D19240")[Morphospecies8]

## E)
pdf(file='Figures/MWtest-pvalbinary_k2_tip.pdf')
heatmap3(fullmatrix8, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols8, 
         ColSideLabs = "Morphospecies",
         labCol=Samples8,labRow=Samples8,
         xlab="Samples")
dev.off()


#--------------------------------------------------------------------------------
# 2.9 Branch length of branches that have a minimum of 4 vertices (br_rate1_long)
#--------------------------------------------------------------------------------

## A)
files9 = list.files(path=".", pattern="manwhit_skel_distances_br_rate1_long_size1000_")
names9 <- substr(files9,1,47) #Create list of file names without the ".csv" part 

## B)
list_matrices9 <-
  lapply(files9, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix9 <- add_matrices(list_matrices9)
fullmatrix9 <- fullmatrix9[!rownames(fullmatrix9) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                           !colnames(fullmatrix9) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples9 <- as.factor(rownames(fullmatrix9))
IDcol9 <-as.data.frame(Samples9, colnames="X1")
colnames(IDs)<- c("Samples9", "Morphospecies")
IDcol9 <-merge(IDcol9,IDs, by  = "Samples9") 
Morphospecies9 <- as.factor(IDcol9$Morphospecies)

cols9 <- c("#114260","#4D724D","#D19240")[Morphospecies9]

## E)
pdf(file='Figures/MWtest-pvalbinary_brlength-long.pdf')
heatmap3(fullmatrix9, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols9, 
         ColSideLabs = "Morphospecies",
         labCol=Samples9,labRow=Samples9,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.10 Branch spacing 1 (brspacing_v1)
#--------------------------------------------------------------------------------

## A)
files10 = list.files(path=".", pattern="manwhit_skel_distances_br_spacing1_size1000_")
names10 <- substr(files10,1,45) #Create list of file names without the ".csv" part 
names10

## B)
list_matrices10 <-
  lapply(files10, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix10 <- add_matrices(list_matrices10)
fullmatrix10 <- fullmatrix10[!rownames(fullmatrix10) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                             !colnames(fullmatrix10) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples10 <- as.factor(rownames(fullmatrix10))
IDcol10 <-as.data.frame(Samples10, colnames="X1")
colnames(IDs)<- c("Samples10", "Morphospecies")
IDcol10 <-merge(IDcol10,IDs, by  = "Samples10") 
Morphospecies10 <- as.factor(IDcol10$Morphospecies)

cols10 <- c("#114260","#4D724D","#D19240")[Morphospecies10]

## E)
pdf(file='Figures/MWtest-pvalbinary_brspacing1.pdf')
heatmap3(fullmatrix10, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols10, 
         ColSideLabs = "Morphospecies",
         labCol=Samples10,labRow=Samples10,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.11 Branch spacing 2 (brspacing_v2)
#--------------------------------------------------------------------------------

## A)
files11 = list.files(path=".", pattern="manwhit_skel_distances_br_spacing2_size1000_")
names11 <- substr(files11,1,45) #Create list of file names without the ".csv" part 
names11

## B)
list_matrices11 <-
  lapply(files11, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix11 <- add_matrices(list_matrices11)
fullmatrix11 <- fullmatrix11[!rownames(fullmatrix11) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                             !colnames(fullmatrix11) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples11 <- as.factor(rownames(fullmatrix11))
IDcol11 <-as.data.frame(Samples11, colnames="X1")
colnames(IDs)<- c("Samples11", "Morphospecies")
IDcol11 <-merge(IDcol11,IDs, by  = "Samples11") 
Morphospecies11 <- as.factor(IDcol11$Morphospecies)

cols11 <- c("#114260","#4D724D","#D19240")[Morphospecies11]

## E)
pdf(file='Figures/MWtest-pvalbinary_brspacing2.pdf')
heatmap3(fullmatrix11, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols11, 
         ColSideLabs = "Morphospecies",
         labCol=Samples11,labRow=Samples11,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.12 Branch angle (br_angle)
#--------------------------------------------------------------------------------

## A)
files12 = list.files(path=".", pattern="manwhit_spheres_angles_end_angle_size1000_")
names12 <- substr(files12,1,43) #Create list of file names without the ".csv" part 
names12

## B)
list_matrices12 <-
  lapply(files12, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix12 <- add_matrices(list_matrices12)
fullmatrix12 <- fullmatrix12[!rownames(fullmatrix12) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                             !colnames(fullmatrix12) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples12 <- as.factor(rownames(fullmatrix12))
IDcol12 <-as.data.frame(Samples12, colnames="X1")
colnames(IDs)<- c("Samples12", "Morphospecies")
IDcol12 <-merge(IDcol12,IDs, by  = "Samples12") 
Morphospecies12<- as.factor(IDcol12$Morphospecies)

cols12 <- c("#114260","#4D724D","#D19240")[Morphospecies12]

## E)
pdf(file='Figures/MWtest-pvalbinary_brangle.pdf')
heatmap3(fullmatrix12, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols12, 
         ColSideLabs = "Morphospecies",
         labCol=Samples12,labRow=Samples12,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.13. Branch thickness at bifurcation/base (da)
#--------------------------------------------------------------------------------

## A)
files13 = list.files(path=".", pattern="manwhit_spheres_angles_end_da_size1000_")
names13 <- substr(files13,1,40) #Create list of file names without the ".csv" part 
names13

## B)
list_matrices13 <-
  lapply(files13, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix13 <- add_matrices(list_matrices13)
fullmatrix13 <- fullmatrix13[!rownames(fullmatrix13) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                             !colnames(fullmatrix13) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples13 <- as.factor(rownames(fullmatrix13))
IDcol13 <-as.data.frame(Samples13, colnames="X1")
colnames(IDs)<- c("Samples13", "Morphospecies")
IDcol13 <-merge(IDcol13,IDs, by  = "Samples13") 
Morphospecies13 <- as.factor(IDcol13$Morphospecies)

cols13 <- c("#114260","#4D724D","#D19240")[Morphospecies13]

## E)
pdf(file='Figures/MWtest-pvalbinary_da.pdf')
heatmap3(fullmatrix13, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols13, 
         ColSideLabs = "Morphospecies",
         labCol=Samples13,labRow=Samples13,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.14 Branch thickness at the center (db)
#--------------------------------------------------------------------------------

## A)
files14 = list.files(path=".", pattern="manwhit_spheres_angles_end_db_size1000_")
names14 <- substr(files14,1,40) #Create list of file names without the ".csv" part 
names14

## B)
list_matrices14 <-
  lapply(files14, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix14 <- add_matrices(list_matrices14)
fullmatrix14 <- fullmatrix14[!rownames(fullmatrix14) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                             !colnames(fullmatrix14) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples14 <- as.factor(rownames(fullmatrix14))
IDcol14 <-as.data.frame(Samples14, colnames="X1")
colnames(IDs)<- c("Samples14", "Morphospecies")
IDcol14 <-merge(IDcol14,IDs, by  = "Samples14") 
Morphospecies14 <- as.factor(IDcol14$Morphospecies)

cols14 <- c("#114260","#4D724D","#D19240")[Morphospecies14]

## E)
pdf(file='Figures/MWtest-pvalbinary_db.pdf')
heatmap3(fullmatrix14, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols14, 
         ColSideLabs = "Morphospecies",
         labCol=Samples14,labRow=Samples14,
         xlab="Samples")
dev.off()

#--------------------------------------------------------------------------------
# 2.15 Branch thickness at the tip (dc)
#--------------------------------------------------------------------------------

## A)
files15 = list.files(path=".", pattern="manwhit_spheres_angles_end_dc_size1000_")
names15 <- substr(files15,1,40) #Create list of file names without the ".csv" part 
names15

## B)
list_matrices15 <-
  lapply(files15, function(x) {
    
    cat("################")
    
    cat(sprintf("Starting to read the file %s\n", x))
    
    single_mat <- as.matrix(read.table(x, sep = ","))
    
    cat("Done!\n")
    
    # Change the rownames
    # first make sure the IDsf is as long as the number of columns
    # and rows
    
    cat("checking that rownames and colnames are of the expected lengths\n")
    
    stopifnot(length(IDsf) == nrow(single_mat))
    stopifnot(length(IDsf) == ncol(single_mat))
    
    cat("OK!\nAssigning the rownames and colnames\n")
    
    # actually change them
    rownames(single_mat) <- IDsf
    colnames(single_mat) <- IDsf
    
    cat("OK!\nTransforming the data\n")
    # Transform data from p-values to dummy variable where alpha=0.05:
    # p-val =< alpha -> Different distributions -> -1
    # p-val > alpha -> Similar distributions -> +1
    
    single_mat[single_mat > 0.05] <- 1
    single_mat[single_mat <= 0.05] <- -1
    
    cat("Done all! returning the result\n")
    # return the result
    return(single_mat)
  })

## C)
fullmatrix15 <- add_matrices(list_matrices15)
fullmatrix15 <- fullmatrix15[!rownames(fullmatrix15) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09"),
                             !colnames(fullmatrix15) %in% c("19Oki01","19Oki02","19Oki05","19Oki07","19Oki09")]

## D)
Samples15 <- as.factor(rownames(fullmatrix15))
IDcol15 <-as.data.frame(Samples15, colnames="X1")
colnames(IDs)<- c("Samples15", "Morphospecies")
IDcol15 <-merge(IDcol15,IDs, by  = "Samples15") 
Morphospecies15 <- as.factor(IDcol15$Morphospecies)

cols15 <- c("#114260","#4D724D","#D19240")[Morphospecies15]

## E)
pdf(file='Figures/MWtest-pvalbinary_dc.pdf')
heatmap3(fullmatrix15, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolors)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols15, 
         ColSideLabs = "Morphospecies",
         labCol=Samples15,labRow=Samples15,
         xlab="Samples")
dev.off()

#------------------------------------------------------------------------------------
# 2.16 Using ALL variables
#------------------------------------------------------------------------------------

## A) Summarize matrices:
add_matrices_vars <- function(...) {
  a <- list(...)
  cols <- sort(unique(unlist(lapply(a, colnames))))
  rows <- sort(unique(unlist(lapply(a, rownames))))
  out <- array(0, dim = c(length(rows), length(cols)), dimnames = list(rows,cols))
  for (m in a) out[rownames(m), colnames(m)] <- out[rownames(m), colnames(m)] + m
  out
}

fullmatrixall <- add_matrices_vars(fullmatrix1, fullmatrix2,
                                   fullmatrix3, fullmatrix4,
                                   fullmatrix5, fullmatrix6,
                                   fullmatrix7, fullmatrix8,
                                   fullmatrix9, fullmatrix10,
                                   fullmatrix11, fullmatrix12,
                                   fullmatrix13, fullmatrix14,
                                   fullmatrix15, fullmatrix15)

# Save output
write.table(fullmatrixall, file="fullmatrixall.txt",
            row.names=TRUE, col.names=TRUE)

## B) Format the matrices
Samples <- as.factor(rownames(fullmatrixall))
IDcol <-as.data.frame(Samples, colnames="X1")
colnames(IDs)<- c("Samples", "Morphospecies")
IDcol <-merge(IDcol,IDs, by  = "Samples") 
Morphospecies <- as.factor(IDcol$Morphospecies)

cols <- c("#114260","#4D724D","#D19240")[Morphospecies]
gradientcolorsall<- c("white","#ecebeb","#a7a7a7","#4d4d4d","#252525")

## C) Plot the heatmap:
heatmap3(fullmatrixall, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolorsall)(1000),
         method="ward.D", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols, 
         ColSideLabs = "Morphospecies",
         labCol=Samples,labRow=Samples,
         xlab="Samples")
dev.off()

#------------------------------------------------------------------------------------
# 2.17 Curvature variables
#------------------------------------------------------------------------------------

## A) Summarize matrices:
fullmatrixcurv <- add_matrices_vars(fullmatrix1, fullmatrix2,
                                   fullmatrix3, fullmatrix4,
                                   fullmatrix5,fullmatrix6,
                                   fullmatrix7, fullmatrix8)

## B) Format the matrices:
Samplescurv <- as.factor(rownames(fullmatrixcurv))
IDcolcurv <-as.data.frame(Samplescurv, colnames="X1")
colnames(IDs)<- c("Samplescurv", "Morphospecies")
IDcolcurv <-merge(IDcolcurv,IDs, by  = "Samplescurv") 
Morphospecies <- as.factor(IDcolcurv$Morphospecies)

cols <- c("#114260","#4D724D","#D19240")[Morphospecies]

## C) Plot the heatmap:
pdf(file='Figures/MWtest-pvalbinary_curvaturevars.pdf')
heatmap3(fullmatrixcurv, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolorsall)(1000),
         method="average", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols, 
         ColSideLabs = "Morphospecies",
         labCol=Samplescurv,labRow=Samplescurv,
         xlab="Samples")
dev.off()

#------------------------------------------------------------------------------------
# 2.18 Branch variables
#------------------------------------------------------------------------------------

## A) Summarize matrices:
fullmatrixbranch <- add_matrices_vars(fullmatrix9, fullmatrix10,
                                      fullmatrix11, fullmatrix12,
                                      fullmatrix13, fullmatrix14,
                                      fullmatrix15, fullmatrix15)

## B) Format the matrices:
Samplesbranch <- as.factor(rownames(fullmatrixbranch))
IDcolbranch <-as.data.frame(Samplesbranch, colnames="X1")
colnames(IDs)<- c("Samplesbranch", "Morphospecies")
IDcolbranch <-merge(IDcolbranch,IDs, by  = "Samplesbranch") 
Morphospecies <- as.factor(IDcolbranch$Morphospecies)

cols <- c("#114260","#4D724D","#D19240")[Morphospecies]

## C) Plot the heatmap:
pdf(file='Figures/MWtest-pvalbinary_branchvars.pdf')
heatmap3(fullmatrixbranch, cexRow=0.4, cexCol=0.4, revC=TRUE, margins=c(5,5),
         col=colorRampPalette(gradientcolorsall)(1000),
         method="average", hclustfun=hclust, symm=TRUE,
         Colv="Rowv", ColSideColors=cols, 
         ColSideLabs = "Morphospecies",
         labCol=Samplesbranch,labRow=Samplesbranch,
         xlab="Samples")
dev.off()

#-------------------------------------------------------------------
# OPTIONAL: Generate a plot for each individual matrix
#-------------------------------------------------------------------

for (i in seq_along(list_matrices)) {
  
  output_filename_i <- paste0(names1[i], ".pdf")
  
  pdf(output_filename_i)
  
  heatmap3(list_matrices[[i]],
           cexRow=0.6, cexCol=0.6, revC=TRUE, 
           col=colorRampPalette(c("blue", "red"))(1000),
           method="ward.D", hclustfun=hclust, symm=TRUE)
  
  dev.off()
  
}
#-------------------------------------------------------------------
