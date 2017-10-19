#' @title Retrieve topTable of differentially expressed genes
#'
#' @description Retrieve topTable of differentially expressed genes
#'
#' @details This function fits a linear model to the gene expression data and retrieves a `topTable` output that has the genes sorted by p value
#'
#' @return A dataframe containing information from the linear fit, ranking DE genes on significance
#'
#' @param DGElist An object of class `DGElist` containing information on gene counts and samples, typical output from the `limma` package
#' @param designMartix A model matrix containing information on the metadata, an output from the `getDesignMatrix` function
#'
#' @export

getTopGenes <- function(DGElist, designMatrix) {
  library(edgeR)
  library(limma)
  library(magrittr)
  library(dplyr)
  library(readr)
  library(tibble)
  library(RColorBrewer)
  library(DESeq)
  library(splitstackshape)
  library(pheatmap)
  library(ggdendro)
  library(plotly)

  if (!is.character(variable)) {
    stop("variable should be character string")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }

  fitGeneExpr <- voom(DGElist, designMatrix, plot = FALSE) %>%
    lmFit(design = designMatrix) %>%
    eBayes()
  nGenesTotal <- nrow(fitGeneExpr)
  slopeContrast <- colnames(fitGeneExpr$design)[2]
  topGenes <- topTable(fitGeneExpr,
                       coef = slopeContrast, number = nGenesTotal) %>%
    rownames_to_column("GeneID")
  topGenes
}

