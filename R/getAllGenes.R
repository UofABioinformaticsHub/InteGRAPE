#' @title Retrieve topTable of differentially expressed genes
#'
#' @description Retrieve topTable of differentially expressed genes
#'
#' @details Paragraph that describes what the function is doing
#'
#' @return This tells you what object it returns
#'
#' @param DGElist describe this
#' @param variable the variable of interest
#' @param metadata describe
#'
#' @export
getAllGenes <- function(DGElist, variable, metadata) {


  if (!is.character(variable)) {
    stop("variable should be character string")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }

  numericVar <- is.numeric(metadata[[variable]])
  int <- dplyr::if_else(numericVar, 1, 0)


  designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
    model.matrix(data = metadata)
  fitGeneExpr <- limma::voom(DGElist, designMatrix, plot = FALSE) %>%
    limma::lmFit(design = designMatrix) %>%
    limma::eBayes()
  nGenesTotal <- nrow(fitGeneExpr)
  slopeContrast <- colnames(fitGeneExpr$design)[2]
  topGenes <- limma::topTable(fitGeneExpr,
                       coef = slopeContrast, number = nGenesTotal) %>%
    tibble::rownames_to_column("GeneID")
  topGenes
}

