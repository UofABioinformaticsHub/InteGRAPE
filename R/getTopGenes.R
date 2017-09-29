#' @title Retrieve topTable of differentially expressed genes
#'
#' @description Retrieve topTable of differentially expressed genes
#'
#' @details This function fits a linear model to the gene expression data and retrieves a `topTable` output that has the genes sorted by p value
#'
#' @return This tells you what object it returns
#'
#' @param DGElist describe this
#' @param variable the variable of interest
#' @param metadata describe
#'
#' @export

getTopGenes <- function(DGElist, variable, metadata) {
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

  numericVar <- is.numeric(metadata[[variable]])
  int <- if_else(numericVar, 1, 0)


  designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
    model.matrix(data = metadata)
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

