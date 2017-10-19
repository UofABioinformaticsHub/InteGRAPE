#' @title Retrieve differential methylation information
#'
#' @description Extracts the top differential methylation tags from a `DGElist`
#'
#' @details Differential methylation across a set of samples are fitted with a generalised linear model in the context of a metadata variable and then ranked in order of significance
#'
#' @param DGElist An object of class `DGElist`, a typical output in the `limma` package.
#' @param designMatrix A model matrix that is the output of the `getDesignMatrix` function
#'
#' @export
#'
getTopMethylation <- function(DGElist, designMatrix) {

  # Fit the data to liner model
  fit <- glmFit(DGElist, designMatrix)

  lrt <- glmLRT(fit, coef = colnames(designMatrix)[2])

  # Get the top Tags for the data
  lrt_top <- topTags(lrt, n = nrow(DGElist$counts), adjust.method = "BH", sort.by = "PValue")
}
