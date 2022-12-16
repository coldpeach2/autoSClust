#' runClust
#'
#' A function that computes Dunn Index and and mean Silhouette score for a clustering output
#'
#' @param srt.data A Seurat object with raw counts in the counts slot.
#' @param features Number of features in the dataset, default is 2000
#' @param res Resolution to find clusters, default is 0.5. The higher the number, the more clusters computed.
#' @return
#'
#' @export
#'
#' @examples
#'
#' @import cluster SingleCellExperiment Seurat clValid scuttle utils
runClust <- function(srt.data, features, res) {

  # Find variable features, identification of highly variable features
  srt.data <- FindVariableFeatures(srt.data, selection.method = "vst", nfeatures = 2000)
  # scale the data
  all.genes <- rownames(srt.data)
  srt.data <- ScaleData(srt.data, features = all.genes)
  # perform linear dimensionality reduction
  srt.data <- RunPCA(srt.data, features = VariableFeatures(object = srt.data))
  # cluster the cells
    # Find neighbours computes a KNN graph based on the euclidean distance in PCA space,
  srt.data <- FindNeighbors(srt.data, dims = 1:10)
    #  Find clusters uses the Louvain algorithm to iteratively group cells together
  srt.data <- FindClusters(srt.data, resolution = 0.5)

  # non linear dimensionality reduction
  srt.data <- RunUMAP(srt.data, dims = 1:30, verbose = FALSE)

  return(srt.data)

}
