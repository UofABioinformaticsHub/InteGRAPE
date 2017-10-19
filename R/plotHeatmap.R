#' @title Plot heatmap
#'
#' @description This function plots a heatmap for direct visualisation of results
#'
#' @details This function will plot a heatmap directly from the count data, an annotation bar at the top of the heatmap will offer information about the plot at a glance. A side bar indicating the pvalue will allow determination of statistical significance at a glance as well.
#'
#' @return A lovely looking heatmap which is interactive
#'
#' @param DGElist A DGElist containing the count and sampling data
#' @param variable The selected variable of interest, should be a character string or vector.
#' @param metadata A dataframe with different variables on the x axis and samples on the y axis.
#' @param nGenes Specify the number of genes you'd like to get data for in the top table at the end
#'
#' @export


plotHeatmap <- function(DGElist, variable, metadata, nGenes = 30) {

  designMatrix <- getDesignMatrix(variable, metadata)
  ScaledCPM <- getScaledCPM(DGElist = DGElist, designMatrix = designMatrix, nGenes = nGenes)
  plottingHeatmap(ScaledCPM = ScaledCPM, variable = variable, metadata = metadata)

}
