#' @title Retrieve topTable of differentially expressed genes
#'
#' @description Retrieve topTable of differentially expressed genes
#'
#' @details This function fits a linear model to the gene expression data and retrieves a `topTable` output that has the genes sorted by p value
#'
#' @return This tells you what object it returns
#'
#' @param DGElist describe this
#' @param designMatrix A model matrix that is the output of the `getDesignMatrix` function
#'
#' @export
getAllGenes <- function(DGElist, designMatrix) {


  if (!is.character(variable)) {
    stop("variable should be character string")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }




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

