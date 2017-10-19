#' @title Multi-dimensional scaling of samples
#'
#' @description Plot the clustering of samples to find any similarities or dissimilarities
#'
#' @details A multi-dimensional scaling plots clusters the data in an unsuperivsed manner, meaning the data input isn't trained with an expected outcome. An example of a method using a trained data set is a Principle Component Analysis. Multi-dimensional scaling is an efficient way to try and observe if there is any differences between samples before executing any formal analyses.
#'
#' @param DGElist An object of class `DGElist`, conatining information about counts and samples
#' @param variable A character string of the variable under investigation, should be found in the metadata
#' @param metadata A data frame containing a set of variables and information relating to them
#' @param size Alter the size of the points, 80 is the standard default size
#'
#' @export
#'
plot_MDS <- function(DGElist, variable, metadata, size = 80) {

  designMatrix <- getDesignMatrix(variable, metadata)

  voomResult <- limma::voom(DGElist, designMatrix, plot = FALSE)

  mds <- limma::plotMDS(voomResult$E,
                         labels = colnames(DGElist$counts))

  mds <- mds[[3]]

  mds <- cbind(mds, 1:nrow(mds))

  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    rowmds <- mds[mds[,3] == x$id,]
    rowmeta <- metadata[metadata$id == x$id,]
    paste0("Variable", ": ", rowmeta$Planting_Date, "<br>",
           "Region", ": ", rownames(rowmds), collapse = "<br />")
  }

  Variable <- metadata[[variable]]
  Region <- DGElist$samples$group
  MDScolour <- grDevices::colorRampPalette(c("green", "blue"))(100)

  MDSofVariables <- mds %>%
    set_colnames(c("Dim1", "Dim2", "id")) %>%
    as.data.frame() %>%
    ggvis::ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
    ggvis::layer_points(size := ~size) %>%
    ggvis::add_legend("fill", title = paste("Variable:", variable),
               orient = "left") %>%
    ggvis::add_legend("shape", title = "Regions",
               orient = "right")


  MDSofVariables %>% ggvis::add_tooltip(all_values, "hover")
}
