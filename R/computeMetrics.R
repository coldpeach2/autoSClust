#' computeMetrics
#' @param srt.data A Seurat object with raw counts in the counts slot.
#'
#' @return The Dunn Index score
#' @export
#'
#' @examples
#' @import Seurat clValid

install.packages("clValid")
library(clValid)
computeMerics <- function(srt.data) {

  countsMat <- as.matrix(srt.data@assays[["RNA"]]@data)
  dunnIndex <- clValid::dunn(clusters=srt.data@meta.data[["seurat_clusters"]], Data=countsMat)

}
