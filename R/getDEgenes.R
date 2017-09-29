#' @title Retrieve differentially expressed genes
#'
#' @description Retrieve a character vector with a list of differentially expressed genes (should have ENSEMBL IDs)
#'
#' @details This function works in a similar way to that of the `getAllGenes` function in that its fitting a linear model to the gene expression data. Then it runs the `topTable` function and retrieves the genes from that data set, ordered by p value.
#'
#' @return A character vector of differentially expressed genes.
#'
#' @param DGElist A DGElist output from `limma`
#' @param variable The variable of interest
#' @param metadata A dataframe with different variables on the x axis and samples on the y axis.
#'
#' @export

# The input for this

getDEgenes <- function(DGElist, variable, metadata) {

  require(edgeR)
  require(limma)
  require(magrittr)
  require(dplyr)
  require(tibble)

  if (!is.character(variable)) {
    stop("variable should be character string")
  }

  if (!is.data.frame(metadata)) {
    stop("metadata must be a data frame")
  }

  numericVar <- is.numeric(metadata[[variable]])
  int <- if_else(numericVar, 1, 0)


  designMatrix <- as.formula(paste("~", paste(1, variable, sep = "+"))) %>%
    model.matrix(data = metadata)
  fitGeneExpr <- voom(DGElist, designMatrix, plot = FALSE) %>%
    lmFit(design = designMatrix) %>%
    eBayes()
  nGenesTotal <- nrow(fitGeneExpr)
  slopeContrast <- colnames(fitGeneExpr$design)[2]
  topGenes <- topTable(fitGeneExpr,
                       coef = slopeContrast, number = nGenesTotal) %>%
    rownames_to_column("GeneID")

  # Find out how many genes are differentially expressed
  # Go by statistical significance
  summary(0.05 >= topGenes$P.Value)

  # As we can see from the result of our summary, there are 4204 differentially expressed genes
  # Now we need to get those 4204 genes
  # First define how we'll subset it
  de <- 0.05 >= topGenes$P.Value

  # Then subset it, easy peesy
  deGenes <- topGenes[de , ]

  # Now try hypergeometric dist funtion

  de <- deGenes$GeneID
  de
}

