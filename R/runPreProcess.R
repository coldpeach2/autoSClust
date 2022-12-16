#' runPreProcess
#'
#' A function that runs user-defined pre processing steps of computing new quality control metrics in the dataset and filtering out uninformative genes
#'
#'
#' @param counts_matrix_path Path that points to the location of raw data folder containing a barcodes.tsv, genes.tsv and matrix.mtx
#' @param filterMetrics String indicating which filtration method to apply to dataset according to quality control metrics.
#' Choices are "percent_mt", "percent_ribo", "percent_disso" or "all. Default is "all"
#'
#' @return Seurat object of filtered dataset
#'
#' @export
#'
#' @details

#'
#' @examples

#'
#'
#' @import Seurat Matrix ggplot2 patchwork

#counts_matrix_path <- file.path("~/pbmc/hg19")
#runPreProcess(counts_matrix_path = counts_matrix_path,filterMetrics="all")
# Create a Seurat object from the raw data and visualize QC metrics

library(Seurat)
library(SeuratObject)
library(SeuratWrappers)
library(Matrix)
library(ggplot2)
library(tidyverse)
runPreProcess <- function(counts_matrix_path, filterMetrics) {
  data_matrix <- Read10X(counts_matrix_path)# use seurat function to read 10X count data
  srt.data <- CreateSeuratObject(data_matrix)

  #Violin Plots
  v1 <- VlnPlot(srt.data, features="nFeature_RNA", cols="pink", pt.size=0.02) + NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5))
  v2 <- VlnPlot(srt.data, features="nCount_RNA", cols="blue", pt.size=0.02) + NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5))


  # min and max genes detected in each cell
  max_nFeature <- max(as.numeric(unlist(srt.data[["nFeature_RNA"]])))
  min_nFeature <- min(as.numeric(unlist(srt.data[["nFeature_RNA"]])))
  # min and max number of molecules detected within a cell
  max_nCount <- max(as.numeric(unlist(srt.data[["nCount_RNA"]])))
  min_nCount <- min(as.numeric(unlist(srt.data[["nCount_RNA"]])))


  ### set thresholds vlines and hlines for ggplot
  nFeature.min <- ggplot2::geom_vline(xintercept=min_nFeature, linetype="dashed", color="darkgrey")
  nFeature.max <- ggplot2::geom_vline(xintercept=max_nFeature, linetype="dashed", color="darkgrey")
  #nUMI -> nCount_RNA
  nCount.x.min <- ggplot2::geom_vline(xintercept=min_nCount, linetype="dashed", color="darkgrey")
  nCount.y.max <- ggplot2::geom_hline(yintercept=max_nCount, linetype="dashed", color="darkgrey")

  #umi lines on y axis
  nCount.x.max <- ggplot2::geom_vline(xintercept=max_nCount, linetype="dashed", color="darkgrey")
  nCount.y.min <- ggplot2::geom_hline(yintercept=min_nCount, linetype="dashed", color="darkgrey")
  f1 <- FeatureScatter(srt.data, feature1="nFeature_RNA", feature2="nCount_RNA", pt.size = 0.3) + scale_color_manual(values = c("#03fc4e", "#666666")) + ggtitle("Genes per Cell vs Count reads per Cell") + nFeature.min + nFeature.max + nCount.y.min + nCount.y.max + NoLegend()

  srt.data <- subset(srt.data, subset =
                       nFeature_RNA > min_nFeature &
                       nFeature_RNA < max_nFeature &
                       nCount_RNA   > min_nCount &
                       nCount_RNA   < max_nCount )

  if ((filterMetrics == "percent_mitochondrial") || (filterMetrics == "all" )){
    srt.data[["percent_mt"]] <- PercentageFeatureSet(srt.data, pattern = "^MT-")
    v3 <- VlnPlot(srt.data, features="percent_mt", cols="green", pt.size=0.02) + NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5))

    max_mito <- max(as.numeric(unlist(srt.data[["percent_mt"]])))
    min_mito <- min(as.numeric(unlist(srt.data[["percent_mt"]])))

    #Mt
    mito.x <- ggplot2::geom_vline(xintercept=max_mito, linetype="dashed", color="darkgrey")
    mito.y <- ggplot2::geom_hline(yintercept=max_mito, linetype="dashed", color="darkgrey")
    f2 <- FeatureScatter(srt.data, feature1="nFeature_RNA", feature2="percent_mt", pt.size = 0.3)+ scale_color_manual(values = c("purple", "#666666")) + ggtitle("Percentage of Mitochondrial Genes") +
      nFeature.min + nFeature.max + mito.y + NoLegend()
  }
  if ((filterMetrics == "percent_ribosomal") || (filterMetrics == "all" )) {
    srt.data[["percent_ribo"]] <- PercentageFeatureSet(srt.data, pattern = '^RPL|^RPS|^MRPL|^MRPS')
    v4 <- VlnPlot(srt.data, features="percent_ribo", cols="orange", pt.size=0.02) + NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5))
    # min and max percentage of ribosomal genes
    max_ribo <- max(as.numeric(unlist(srt.data[["percent_ribo"]])))
    min_ribo <- min(as.numeric(unlist(srt.data[["percent_ribo"]])))
    #Ribo
    ribo.x <- ggplot2::geom_vline(xintercept=max_ribo, linetype="dashed", color="darkgrey")
    ribo.y <- ggplot2::geom_hline(yintercept=max_ribo, linetype="dashed", color="darkgrey")
    f4 <- FeatureScatter(srt.data, feature1="nFeature_RNA", feature2="percent_ribo", pt.size = 0.3)+ scale_color_manual(values = c("green", "#666666")) + ggtitle("Percentage of Ribosomal Genes") +
      nFeature.min + nFeature.max + ribo.y  + NoLegend()
  }

    if ((filterMetrics == "percent_dissociation") || (filterMetrics == "all" )){
      # calculate percentage of mitochondria, ribosomal and dissociation genes
      dissociation <- c("ACTG1","ANKRD1","ARID5A","ATF3","ATF4","BAG3","BHLHE40","BRD2","BTG1","BTG2","CCNL1",
                        "CCRN4L","CEBPB","CEBPD","CEBPG","CSRNP1","CXCL1","CYR61","DCN","DDX3X","DDX5","DES",
                        "DNAJA1","DNAJB1","DNAJB4","DUSP1","DUSP8","EGR1","EGR2","EIF1","EIF5","ERF","ERRFI1",
                        "FAM132B","FOS","FOSB","FOSL2","GADD45A","GADD45G","GCC1","GEM","H3F3B","HIPK3",
                        "HSP90AA1","HSP90AB1","HSPA1A","HSPA1B","HSPA5","HSPA8","HSPB1","HSPE1","HSPH1",
                        "ID3","IDI1","IER2","IER3","IER5","IFRD1","IL6","IRF1","IRF8","ITPKC","JUN","JUNB",
                        "JUND","KCNE4","KLF2","KLF4","KLF6","KLF9","LITAF","LMNA","MAFF","MAFK","MCL1","MIDN",
                        "MIR22HG","MT1","MT2","MYADM","MYC","MYD88","NCKAP5L","NCOA7","NFKBIA","NFKBIZ","NOP58",
                        "NPPC","NR4A1","ODC1","OSGIN1","OXNAD1","PCF11","PDE4B","PER1","PHLDA1","PNP","PNRC1",
                        "PPP1CC","PPP1R15A","PXDC1","RAP1B","RASSF1","RHOB","RHOH","RIPK1","SAT1","SBNO2","SDC4",
                        "SERPINE1","SKIL","SLC10A6","SLC38A2","SLC41A1","SOCS3","SQSTM1","SRF","SRSF5","SRSF7",
                        "STAT3","TAGLN2","TIPARP","TNFAIP3","TNFAIP6","TPM3","TPPP3","TRA2A","TRA2B","TRIB1",
                        "TUBB4B","TUBB6","UBC","USP2","WAC","ZC3H12A","ZFAND5","ZFP36","ZFP36L1","ZFP36L2","ZYX")
      srt.data[["percent_disso"]] <- PercentageFeatureSet(srt.data, pattern = dissociation)
      v5 <- VlnPlot(srt.data, features="percent_disso", cols="purple", pt.size=0.02) + NoLegend() + theme(axis.text.x=element_text(angle=0, hjust=0.5))
      # min and max percentage of dissociation genes
      max_dis <- max(as.numeric(unlist(srt.data[["percent_disso"]])))
      min_dis <- min(as.numeric(unlist(srt.data[["percent_disso"]])))

      #Disso
      disso.x <- ggplot2::geom_vline(xintercept=max_dis, linetype="dashed", color="darkgrey")
      disso.y <- ggplot2::geom_hline(yintercept=max_dis, linetype="dashed", color="darkgrey")
      f5 <- FeatureScatter(srt.data, feature1="nFeature_RNA", feature2="percent_disso", pt.size = 0.3)+ scale_color_manual(values = c("#590849", "#666666")) + ggtitle("Percentage of Dissociation Genes") +
        nFeature.min + nFeature.max + disso.y  + NoLegend()
    }

    #plots <- list(v1, f1, v2)

    #return(plots)
    if (exists("v3") && exists("v4") && exists("v5")) {
      #print(patchwork::wrap_plots(v1 + f1 | v2 | v3 + f2 | v4 + f4 | v5 + f5 + plot_layout(guides = 'collect') + NoLegend()) +
             # plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))

      srt.data <- subset(srt.data, subset =
                           percent_mt   < max_mito &
                           percent_ribo < max_ribo &
                           percent_disso < max_dis)
      final_plot <- list(v1, f1, v2, v3, f2, v4, f4, v5, f5, srt.data)
      return(final_plot)
    }
    if (exists("v3")) {
      #print(patchwork::wrap_plots(v1 + f1 | v2 | v3 + f3 + plot_layout(guides = 'collect') + NoLegend()) +
             # plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))
      srt.data <- subset(srt.data, subset = percent_mt   < max_mito)
      final_plot <- list(v1, f1, v2, v3, f2, srt.data)
      return(final_plot)

    }
  if (exists("v4")) {
    #print(patchwork::wrap_plots(v1 + f1 | v2 | v4 + f4 + plot_layout(guides = 'collect') + NoLegend()) +
            #plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))
    srt.data <- subset(srt.data, subset = percent_ribo < max_ribo )
    final_plot <- list(v1, f1, v2, v4, f4, srt.data)
    return(final_plot)

  }
  if (exists("v5")) {
    #print(patchwork::wrap_plots(v1 + f1 | v2 | v5 + f5 + plot_layout(guides = 'collect') + NoLegend()) +
            #plot_annotation(theme=theme(plot.title = element_text(hjust = 0.5, face="bold"))))
    srt.data <- subset(srt.data, subset = percent_disso < max_dis )
    final_plot <- list(v1, f1, v2, v5, f5, srt.data)
    return(final_plot)
  }

}


