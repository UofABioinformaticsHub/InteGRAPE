#' @title Create Design Matrix
#'
#' @description This function is just a quick way to make a design matrix
#'
#' @details `getDesignMatrix` creates a design matrix from the given `metadata` and `variable` of interest.
#' @details This function defined a formula from the supplied information and creates the matrix directly from that information.
#' @detailts The intercept is removed when the input is a categorical variable but retained when the variable is continuous.
#'
#' @return A model matrix that is compatible with `limma` and `InteGRAPE` analyses with information on the current variable of interest
#'
#' @param variable The selected variable of interest, should be a character string or vector.
#' @param metadata A dataframe with different variables on the x axis and samples on the y axis.
#'
#' @export


getDesignMatrix <- function(variable, metadata) {

  numericVar <- is.numeric(metadata[[variable]])
  int <- if_else(numericVar, 1, 0)
  designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>%
    model.matrix(data = metadata)
  designMatrix

}
