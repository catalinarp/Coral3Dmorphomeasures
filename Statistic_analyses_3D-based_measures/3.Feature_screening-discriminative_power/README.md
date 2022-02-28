# <b>Feature screening and discriminative power of 3D-based measures</b>
+ Variable screening and selection:
  + “Preliminary selected” subset according to interspecific differences (ANOVA and two-sample tests; https://github.com/catalinarp/Coral3Dmorphomeasures/blob/main/Statistic_analyses_3D-based_measures/1.Univariate_values/InterspecificDifferences.R)
  + Multivariate analysis of variance - MANOVA (stats R package v4.1.0) 
  + Box plots (ggplot2 R package v3.3.5) 
  + Assessing correlation:
    + Pearson correlation coefficients (Hmisc R package v4.5-0)
    + Correlation plot (psych R package v2.1.6)
  + Clustering inspection:
    + Hierarchical clustering analysis - HCA (cluster R package v2.1.2)
    + Variable selection for Gaussian model-based clustering (clustvarsel R package v2.3.4) and principal component analysis - PCA (stats R package v4.1.0)
+ Discriminative power of 3D-based measures:
  + Selection of variables with low variance inflation factor (VIF < 10; usdm R package v1.1-18) and linear discriminant analysis - LDA (MASS R package v7.3-54) 
  + Discriminant analysis based on Gaussian finite mixture modelling (mclust R package v5.4.7)
