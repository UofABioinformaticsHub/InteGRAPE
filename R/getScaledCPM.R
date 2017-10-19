#' @title Retrieve scaled counts per million
#'
#' @description This function retrieves counts per million that have been scaled to allow more effective visualisation
#'
#' @details This function will take the count data from a DGElist and scale it in order to enhance the visualisation appearance, especially when making heatmaps. The scaling method mimics that of the `pheatmap` function.
#'
#' @return A molten dataframe containing scaled cpm values
#'
#' @param DGElist A DGElist containing the count and sampling data
#' @param designMatrix A model matrix that is the output of the `getDesignMatrix` function
#' @param nGenes Specify the number of genes you'd like to get data for in the top table at the end
#'
#' @export


# library(dplyr)

getScaledCPM <- function(DGElist, designMartix, nGenes = 30) {

  require(edgeR)
  require(limma)
  require(magrittr)
  require(dplyr)
  require(readr)
  require(ggplot2)
  require(tibble)
  require(RColorBrewer)
  require(DESeq)
  require(splitstackshape)

  # Fitting models
  ## Voom method

  ### Fit gene expression to linear model using lmFit

  fitGeneExpr <- limma::voom(DGElist, designMatrix, plot = FALSE) %>%
    limma::lmFit(design = designMatrix) %>%
    limma::eBayes()
  nGenesTotal <- nrow(fitGeneExpr)
  slopeContrast <- colnames(fitGeneExpr$design)[2]
  topGenes <- limma::topTable(fitGeneExpr,
                       coef = slopeContrast, number = nGenesTotal) %>%
    tibble::rownames_to_column("GeneID")


  ## Define gene and sample order to subset for in the topGenes object
  inputVals <- metadata[[currentVar]]
  sampleOrder <- order(inputVals)

  GeneOrder <- order(topGenes$GeneID)
  topGenes <- topGenes[GeneOrder, ] %>%
    tibble::as_data_frame()


  ### Select only for FDR and gene names

  geneExprPValFull <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)


  ## Retireve Gene IDs for all samples
  #First we want to get the gene names and FDR for all of the samples, this way we can make an FDR bar that will give a colour indication of the FDRcutoff

  topGenesTotal <- rownames(limma::topTable(fitGeneExpr, coef = slopeContrast, number = nGenesTotal))
  inputValsTotal <- metadata[[currentVar]]
  sampleOrderTotal <- order(inputValsTotal)
  #inputVals <- inputVals[sampleOrder]


  ### Get the cpm for all values
  #This is the step required for the FDR

  geneExprCPM2plot <- cpm(DGElist, log = TRUE)[topGenesTotal, sampleOrder]
  colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)

  colnames(geneExprCPM2plot) <- colnames(geneExprCPM2plot) %>% gsub(pattern = "_Kaesler", replacement = "8") %>%
    gsub(pattern = "_Elder_Nur", replacement = "6") %>%
    gsub(pattern = "_Nuriootpa", replacement = "5") %>%
    gsub(pattern = "_Hens_Keyn", replacement = "22") %>%
    gsub(pattern = "_Hens_Eden", replacement = "21") %>%
    gsub(pattern = "_Heathvale", replacement = "20")


  ## Scale by rows
  #In this function we are calculating the mean and standard deviation for every row of the matrix. Then we subtract the mean from the raw value and divide by the calculated standard deviation to generate our scaled values. Dividing by the standard deviation means that

  scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
  }



  ### Generate CPM dataframe

  nGenes <- nGenes
  topGenes <- limma::topTable(fitGeneExpr, coef = slopeContrast, number = nGenes)
  topGenes <- tibble::rownames_to_column(topGenes, "GeneID")

  #inputVals <- inputVals[sampleOrder]

  topGenes_names <- rownames(limma::topTable(fitGeneExpr, coef = slopeContrast, number = nGenes))
  inputVals <- metadata[[currentVar]]
  sampleOrder <- order(inputVals)


  ### Get the cpm

  geneExprCPM2plot <- cpm(DGElist, log = TRUE)[topGenes_names, sampleOrder]
  colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)

  colnames(geneExprCPM2plot) <- colnames(geneExprCPM2plot) %>% gsub(pattern = "_Kaesler", replacement = "8") %>%
    gsub(pattern = "_Elder_Nur", replacement = "6") %>%
    gsub(pattern = "_Nuriootpa", replacement = "5") %>%
    gsub(pattern = "_Hens_Keyn", replacement = "22") %>%
    gsub(pattern = "_Hens_Eden", replacement = "21") %>%
    gsub(pattern = "_Heathvale", replacement = "20")

  scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)

  dendogramObject <- as.dendrogram(hclust(dist(scaledCPMgeneExpr), method = "ward.D2"))
  clusterOrder <- order.dendrogram(dendogramObject)
  clusteredMatrix <- scaledCPMgeneExpr[clusterOrder, ]
  geneExprCPM2ggplot <- reshape2::melt(clusteredMatrix, value.name = "scaledCPM")
  colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")

  ### Select only for FDR and gene names

  geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)
  GeneOrder <- levels(geneExprCPM2ggplot$GeneID)
  GeneOrder <- data.frame(GeneOrder, 1)
  rownames(GeneOrder) <- GeneOrder$GeneOrder
  colnames(GeneOrder) <- c("GeneID", "boo")
  geneExprPVal <- dplyr::left_join(GeneOrder, geneExprPVal, by = "GeneID")
  geneExprPVal$boo <- NULL
  geneExprPVal$booo <- NULL

  #geneExprCPM2ggplot <- reshape2::melt(scaledCPMgeneExpr, value.name = "scaledCPM")
  #colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")

  geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPVal, by = 'GeneID')
  geneExprCPM2ggplotPVal$GeneID <- as.factor(geneExprCPM2ggplotPVal$GeneID)
  # levels(geneExprCPM2ggplotPVal$GeneID) <- geneExprCPM2ggplotPVal$GeneID

  geneExprCPM2ggplotPVal

}

