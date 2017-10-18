#' @title Run Fisher's Method to combine p values
#' 
#' @description Fisher's Method combines p values from different statistical models.
FishersMethod <- function(GOde, GOdm) {
  
  # GOdm$P.DM <- GOdm$P.DE
  # GOdm$P.DE <- NULL
  
  GO_Combined <- dplyr::left_join(GOdm, GOde, by = c("Term", "Ont", "N"))
  
  CombinedPVals <- data_frame(GO_Combined$P.DE.x, GO_Combined$P.DE.y)
  
  PvalsCombined <- (1 - pchisq(rowSums(-2 * log(CombinedPVals)), df = 2 * length(CombinedPVals)))

  GO_Combined$Combined_Pvalues <- PvalsCombined
  
  GO_Combined
}


