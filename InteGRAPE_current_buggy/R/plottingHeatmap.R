#' @title Plot heatmap
#'
#' @description This function plots a heatmap from the scaled cpm values
#'
#' @details This function will plot a heatmap directly from the count data, an annotation bar at the top of the heatmap will offer information about the plot at a glance. A side bar indicating the pvalue will allow determination of statistical significance at a glance as well.
#'
#' @return A lovely looking heatmap which is interactive
#'
#' @param ScaledCPM A DGElist containing the scaled count data, the output from the `getScaledCPM` function
#' @param variable The selected variable of interest, should be a character string or vector.
#' @param metadata A dataframe with different variables on the x axis and samples on the y axis.
#'
#' @export

plottingHeatmap <- function(ScaledCPM, variable, metadata) {

  # Plotting

  currentVar <- variable

  ## Set colour information

  nColours <- 101
  varRange <- range(appmetaG[[currentVar]])
  varGradient <- floor(seq(varRange[1], varRange[2], length.out = nColours))
  varPositions <- findInterval(appmetaG[[currentVar]], varGradient)
  varPalette <- grDevices::colorRampPalette(c("brown", "yellow"))(nColours)
  varColours <- varPalette[varPositions]

  saveRDS(varPalette, "~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/varPalette.rds")
  nColoursPVal <- nGenesTotal
  PValrange <- range(geneExprPValFull$P.Value)
  PValgradient <- seq(from = PValrange[1], to = PValrange[2], length.out = nColoursPVal)
  PValpositions <- findInterval(geneExprPVal$P.Value, PValgradient)
  PValpalette <- grDevices::colorRampPalette(c("green", "light green", "yellow", "pink", "red", "dark red"))(nColoursPVal)
  PValcolours <- PValpalette[PValpositions]

  heatmapColours <- grDevices::colorRampPalette(c("dark blue", "blue", "cyan", "white", "pink", "red", "dark red"))(540)

  ## Create colour bars
  ### Create data for colour bar

  i <- nrow(unique(appmetaG[currentVar]))

  if(int == 1) {
    x <- seq(from = min(appmetaG[currentVar]), to = max(appmetaG[currentVar]), length.out = nColours)
    y <- 1
  } else {
    x <- 1:i
    y <- 1
  }

  xy <- tibble::data_frame(x, y)

  # x <- as.factor(appmeta[[currentVar]])
  # plyr::revalue(x, c("Central_Grounds" = 1, "Eden_Valley" = 2))


  ### Create data for FDR bar

  a <- 1:101
  b <- 1
  c <- tibble::data_frame(a, b)


  ### Create data for FDR bar title

  texta <- 1
  textb <- 1
  text <- tibble::data_frame(texta, textb)


  ### Make data for variables ordered

  orderedVariables <- inputVals[sampleOrder]
  yAxis <- 1
  variableAnnotation <- tibble::data_frame(orderedVariables, yAxis)



  ### Organise title

  title <- if(int == 1) {
    paste("Top", nGenes, "Genes - Gene Expression Analysis of", currentVar, ":", round(min(appmetaG[currentVar]), 2), "to",
          round(max(appmetaG[currentVar]), 2))
  } else {
    paste("Analysis by", currentVar)
  }

                       ############################
                       ###   Plot colour bars   ###
                       ############################

  colourBarXaxis <- levels(geneExprCPM2ggplotPVal$Vineyard)
  colourBarXaxis <- as.factor(colourBarXaxis)

  colourBarDF <- tibble::data_frame(colourBarXaxis, orderedVariables)
  colourBarDF$colourBarXaxis <-  factor(colourBarDF$colourBarXaxis, levels = colourBarDF$colourBarXaxis)


                  #######################################
                  #####   Variable Annotation bar   #####
                  #######################################

  varBar <- ggplot2::ggplot(data = colourBarDF, aes(x = colourBarXaxis, y = orderedVariables, fill = orderedVariables)) +
    ggplot2::geom_bar(stat = "identity", width = 1) +
    ggplot2::scale_fill_gradientn(colours = varPalette) +
    ggplot2::coord_cartesian(ylim = c(min(colourBarDF$orderedVariables),
                             (max(colourBarDF$orderedVariables)))) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_void()


                   ####################################################
                   ####   I'm not sure if I use this plot at all   ####
                   ####################################################

  colourBar <- ggplot2::ggplot(data = xy, aes(
    x = seq(from = 0, to = 18.99, length.out = 101),
    y = y, fill = x)) +
    ggplot2::geom_raster(stat = "identity", position = "identity") +
    ggplot2::scale_fill_gradientn(colours = varPalette) +
    ggplot2::ggtitle(paste("Gene Expression Analysis of", currentVar, "from",
                  round(min(appmetaG[[currentVar]]), 2), "to",
                  round(max(appmetaG[[currentVar]]), 2))) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::expand_limits(FALSE) +
    ggplot2::theme_void()
                               ##############################
                               ###   P.value colour bar   ###
                               ##############################

  geneExprPVal <- tibble::as.tibble(geneExprPVal)
  geneExprPVal$GeneID <- factor(geneExprPVal$GeneID)
  geneExprCPM2ggplot <- tibble::as.tibble(geneExprCPM2ggplot)
  levels(geneExprPVal$GeneID) <- levels(geneExprCPM2ggplot$GeneID)

  PValbarWhole <- ggplot2::ggplot(data = geneExprPVal,
                         aes(x = currentVar, y = GeneID, fill = -log10(P.Value))) +
    ggplot2::geom_raster(stat = "identity", position = "identity") +
    ggplot2::scale_fill_gradientn(colours = PValcolours) +
    ggplot2::guides(fill = FALSE) +
    ggplot2::theme_void()

  #ggplot2::ggplot(geneExprPVal)

  PValtitle <- ggplot2::ggplot(data = text, aes(x = texta, y = textb)) +
    ggplot2::geom_text(label = "P Values") +
    ggplot2::theme_void()

  # -
                               ####################################
                               ###  Create interactive heatmap  ###
                               ####################################


  ### Prepare hide object used to remove labels and lines
  # This object will allow us to remove any excess information which will clutter our plot. In doing so, information regarding
  # the samples and genes will be accessible through the hover function. Turning on the sparklines in our heatmap will allow
  # us to keep track of what sample and vineyard it is from
  # (In the shiny app, I plan to add a statement that will allow the user to turn axes on and off)
  # Don't really need to do this for methylation, but I just made the code and it looked pretty

  hide <- list(
  title = "",
  zeroline = FALSE,
  showline = FALSE,
  showticklabels = FALSE,
  showgrid = FALSE
  )


  ### Create the interactive heatmap
  #The values are in log-cpm, which means that the 0 baseline equals 1cpm and any counts with less than 1 cpm are represented as negative.

  GeneExprheatmap <- ggplot2::ggplot(data = geneExprCPM2ggplot,
  aes(x = Vineyard, y = GeneID, fill = scaledCPM)) +
  ggplot2::geom_raster(stat = "identity", position = "identity") +
  ggplot2::xlab("Vineyard") +
  ggplot2::scale_fill_gradientn(colours = heatmapColours) +
  ggplot2::theme(axis.text.x = element_text(angle = 90, hjust = 1))

  GeneExprheatmap <- GeneExprheatmap %>% plotly::ggplotly()


  # This is the original subplot with the two colour bars

  GeneExprheatmapSubplot <- plotly::subplot(varBar, PValtitle, GeneExprheatmap, PValbarWhole,
  nrows = 2, heights = c(0.1, 0.75),
  widths = c(0.8, 0.1),
  shareX = TRUE, shareY = TRUE, titleX = FALSE,
  titleY = FALSE)
  plotly::layout(GeneExprheatmapSubplot, title = title, showlegend = TRUE,
  margin = list(l = 140,
  r = 50,
  b = 60,
  t = 40),
  plot_bgcolor = "white")

}
