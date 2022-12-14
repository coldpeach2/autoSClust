# autoSClust

## Description
Deploy the scRNA sequencing data analysis pipeline for downstream cell level and gene level expression analysis. This package is meant to accommodate a broad range of scRNA-seq data types and output the best method according to the nature of the dataset.
## Instalation
To install the latest version of the package:
``` r
require("devtools")
devtools::install_github("<coldpeach2>/<autoSClust>", build_vignettes = TRUE) 
library("<autoSClust>")
```
To run the shinyApp: Under construction

## Overview
``` r
ls("package:<autoSClust>")
      data(package = "<autoSClust>") # optional
      browseVignettes("<autoSClust>")
```
## Contributions
autoSClust currently consists of 4 functions that facilitate scRNA-seq data analysis

The first function `runPreProcess` converts the counts matrix outputs from CellRanger into a Seurat object. Then performs the first pre-processing steps including filtering out genes with insignificant gene expression. Next normalization takes place, where 3 different normalization techniques, SCTransform, SCnorm and LogNorm and are trialed and the best one is returned according to its performance on correctness (Rand index), compactness (Dunn index), and robustness 

The second function `runNormalization` will output 3 different normalization techniques and output the results of the different normalization methods along with their corresponding graphs for visualization purposes. 

Next, `runClust` will take in the three objects and performs the necessary preparations for clustering

The function `computeMetrics` will compute the Rand Index and Dunn index and display the results in a graph.

Finally, 


Next `plotData` is deployed to visualize the output of the corrected genes. It expects a normalized scRNA-seq dataset and returns a plot. However, before the genes are plotted, feature selection is run to further reduce the technical noise of the normalized data and will be using Seurat's FindVariableFeatures function. Further, after feature selection, dimensionality reduction is required to additionally reduce un-informative genes. Dimensionality reduction falls into two main categories; linear and non-linear. Linear dimensionality reduction includes Principle component analysis (PCA) while non-linear methods includes Uniform Approximation and Projection method (UMAP) and t-SNE. UMAP has proven to be the better alternative to transform the data in dimensionally reduced space compared to PCA or t-SNE. This function will output PCA, t-SNE and UMAP to visualize the data after normalization and dimensionality reduction. 



## References

## Acknowledgements
