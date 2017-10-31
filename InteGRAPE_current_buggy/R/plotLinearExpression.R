#' @title Retrieve differentially methylated genes
#'
#' @description Retrieve a character vector with a list of differentially methylated genes (Should have ENSEMBL IDs)
#'
#' @details This function fits a linear model to the differential expression data and visualises that model in a scatter plot.
#'
#' @return A plot visualising the linear model of gene expression.
#'
#' @param DGElist A DGElist output from `limma`
#' @param designMatrix A model matrix with information on the variable of interest, the output of `getDesignMatrix`
#' @param allGenes An output with the `topTable` from the `getAllGenes` function
#' @param selectGene This is a subset operator, it takes an integer or a character vector of a gene ID of interest
#'
#' @export

### Visualse linear representation of gene expression levels
#Here is the code that can be placed into the app to select a gene and visualise its linear trend


plotLinearExpression <- function(DGElist, designMatrix, allGenes, selectGene = 1) {

  voomObj <- limma::voom(DGElist, designMatrix, plot = FALSE)
  selectGene <- 1
  Gene <- allGenes[[selectGene, 1]]
  pvalue <- geneExprPValFull[[selectGene, 2]]
  GeneCounts <- voomObj$E[Gene, ]

  LinearRep <- tibble::tibble(factor(names(GeneCounts[sampleOrder]), levels =
                               unique(names(GeneCounts[sampleOrder]))), GeneCounts[sampleOrder])
  LinearRep$Variable <- inputValsTotal[sampleOrderTotal]
  colnames(LinearRep) <- c("Region", "Counts", "Variable")

  #OrderedLinearRep <- LinearRep[sampleOrder, ]
  #levels(OrderedLinearRep$Region) <- unique(names(GeneCounts))

  plotLinearRep <- ggplot2::ggplot(data = LinearRep, aes(x = Variable, y = Counts)) +
    ggplot2::geom_point(stat = "identity")
  plotLinearRep + ggplot2::geom_smooth(method = "lm") + ggtitle(paste(Gene, ": p =", round(pvalue, 11))) + ggplot2::theme_bw()
}
