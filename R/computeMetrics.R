#' computeUnsupervisedMetrics
#'
#' A function that computes Dunn Index and and mean Silhouette score for a clustering output
#'
#' @param srt.data A Seurat object with raw counts in the counts slot.
#'
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
#'

computeMerics <- function(srt.data) {


}
