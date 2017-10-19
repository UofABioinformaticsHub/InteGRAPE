#' @title Normalisation of reads
#'
#' @description Normalise the data to account for any variation
#'
#' @details Experimental batch effects occur when data has been collected from different experiments, which may lead to variation in the data as a result of those independent experiments. Normalisation of the data will normalise this variation and effectively eliminate any variance from batch effects.
#'
#' @param DGElist An object of class `DGElist`, typically an output of the `limma::DGElist` function
#'
#' @export
normalizeReads <- function(DGElist) {

  # Save an unnormalised version of the data
  DGElist2 <- DGElist

  # Plot unnormalised
  par(mfrow = c(1, 2))
  lcpm <- cpm(DGElist2, log = TRUE)
  boxplot(lcpm, las = 2, col = col, main = "")
  title(main = "Unnormalised data", ylab = "Log-cpm")

  # Calculate Normal Factors, this is the normalisation step
  DGElist <- calcNormFactors(
    DGElist,
    method = "TMM"
  )

  # Plot these results
  lcpm <- cpm(DGElist, log = TRUE)

  boxplot(lcpm, las = 2, col = col, main = "")

  title(main = "Normalised data", ylab = "Log-cpm")

  DGElist
}
