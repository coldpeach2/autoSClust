## ----setup--------------------------------------------------------------------

library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(tidyverse)

#NOTE!!! because I have configured my functions to run with the Shiny app, this vignette will no longer run as expected !
#TO DO: Please supply the path below to which the package is downloaded to run this function!

#counts_matrix_path <-

runPreProcess(counts_matrix_path, filterMetrics="percent_mitochondrial")


runNormalization(srt.data,  norm="all")

runClust(srt.data, features, res)






