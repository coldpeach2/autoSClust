
#' runNormalization
#'
#' A function that takes in user defined normalization methods and outputs a density plot
#'
#' @param srt.data A Seurat object
#' @param norm Character vector of normalization methods.Choices are "Log", "SCT" or "CPM". Default is "Log"
#'
#' @return seurat object of normalized gene counts
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



runNormalizarion <- function(srt.data, norm){
  #SCTransform, SCnorm and LogNorm

  if ((norm == "log_norm") || (norm == "all" )){
    # log normalize data
    data.Log <- NormalizeData(srt.data, normalization.method = "LogNormalize", scale.factor = 10000)

    plot.title3 = "Log Normalized data"
    data.Log_plot <- ggplot2::ggplot(data.Log@meta.data,
                                     ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="#23e809", fill="#189407") +
      ggplot2::ggtitle(plot.title3) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) +
      ggplot2::xlim(0,2000) + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
  }
  if ((norm == "SCTransform") || (norm == "all" )){
    # run sctransform
    data.SCT <- SCTransform(srt.data, vars.to.regress = "percent_mt", verbose = FALSE)

    plot.title2 = "SCTransform Normalized data"
    data.SCT_plot <- ggplot2::ggplot(data.SCT@meta.data,
                                     ggplot2::aes(x=nFeature_SCT, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="#a34fc2", fill="#7b16a1") +
      ggplot2::ggtitle(plot.title2) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())
  }
  if ((norm == "counts_per_million") || (norm == "all" )){
    # CPM normalization
    data.CPM <- NormalizeData(srt.data, normalization.method = "RC", scale.factor = 1e6)

    plot.title1 = "Counts Per Million Normalized data"
    data.CPM_plot <- ggplot2::ggplot(data.CPM@meta.data,
                                     ggplot2::aes(x=nFeature_RNA, color=orig.ident, fill=orig.ident)) +
      ggplot2::geom_density(size=0.5, alpha=0.2, color="#02bdc4", fill="#05797d") +
      ggplot2::ggtitle(plot.title1) + theme(plot.title = element_text(hjust = 0.5, face ="bold")) +
      ggplot2::xlim(0,2000) + NoLegend() +
      theme(axis.title.y = element_blank(),axis.text.y = element_blank(),axis.ticks.y = element_blank())

  }
  if (norm == "all"){
    srt.data <- NormalizeData(srt.data, normalization.method = "LogNormalize", scale.factor = 10000)
    plots <- list(data.Log_plot, data.SCT_plot, data.CPM_plot, srt.data)
    return(plots)

  }
  if (norm == "counts_per_million"){
    lst <- list(data.CPM_plot, data.CPM)
    return(lst)

  }
  if (norm == "log_norm"){
    lst <- list(data.Log_plot, data.Log)
    return(lst)

  }
  if (norm == "SCTransform"){
    lst <- list(data.SCT_plot, data.SCT_plot)
    return(lst)

  }


}
