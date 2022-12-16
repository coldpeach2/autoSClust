#' runClust
#'
#' A function that computes Dunn Index and and mean Silhouette score for a clustering output
#'
#' @param srt.data A Seurat object with raw counts in the counts slot.
#' @param features
#' @param dims
#' @param res
#' @return A list with the Dunn Index and mean Silhouette score as elements
#'
#' @export
#'
#' @examples
#' # Example:
#' \dontrun{
#' data(embryo)
#' data(embryoClusts)
#' # Compute metrics for the first clustering output on the embryo dataset
#' metrics <- computeUnsupervisedMetrics(embryo, embryoClusts[[1]])
#' }
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
