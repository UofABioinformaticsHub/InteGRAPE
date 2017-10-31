
# This is the server logic for a Shiny web application.
# You can find out more about building applications with Shiny here:
#
# http://shiny.rstudio.com
#

# library(edgeR)
# library(limma)
# library(ggvis)
# library(magrittr)
# library(dplyr)
# library(tidyr)
# library(readr)
# library(tibble)
# library(RColorBrewer)
# library(DESeq)
# library(splitstackshape)
# library(rtracklayer)
# library(plotly)
# library(SummarizedExperiment)
# library(IRanges)
# library(GenomicRanges)
# library(GenomicFeatures)
# library(Biostrings)
# library(shiny)
# library(shinycssloaders)
# 
# 
# source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getTopGenes.R")
# source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDEgenes.R")
# source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDMgenes.R")
# source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getGOenrichment.R")
# source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotGOterms.R")
# source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotVariableGO.R")

library(InteGRAPE)
library(shiny)
library(magrittr)
library(leaflet)
library(rgdal)
library(raster)

server <- function(input, output, session) {
  

  
  ###########################################################
  ###   Prepare all the data required for the shiny app   ###
  ###########################################################
  # 
  # # Load in the Methylation data
  # # DGE_Methylation <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_Methylation.rds")
  # 
  # # Load in the RNA sequencing data
  # load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_RNAseq.RData")
  # 
  # # Load in the gene expression meta-data
  # load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/appmetaG.RData")
  # 
  # # Load in the methylation meta-data
  # load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/appmetaM.RData")
  # 
  # # Load in the gff file
  # gffFile <- file.path("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/Annotation_Files/Vitis_vinifera_annotation.gff")
  # file.exists(gffFile)
  # genesGR <- import.gff(gffFile) %>%
  #   subset(type == "CDS") %>%
  #   split(f = .$CDS)
  # 
  # # Load in the GRanges
  # GR <- file.path("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/Annotation_Files/Vitis_vinifera.IGGP_12x.36.gff3.gz")
  # file.exists(gffFile)
  # genesGR <- import(gffFile)
  # 
  # # Load in the TSV file
  # TSVfile <- file.path("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/Annotation_Files/mart_export.txt")
  # file.exists(TSVfile)
  # geneloci <- read_tsv(TSVfile)
  # 
  # # Load in the Ordered methylation counts
  # OrderedSubsetMethCounts <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/OrderedSubsetMethCounts.rds")
  # 
  # # Load in the methylation ranges file
  # MethylationInGenes <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/MethylationInGenes.rds")
  # 
  # # Load in the methylation phenotypic data
  # phenoDataM <- read_csv(file = "~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/phenoDataM.csv")
  # 
  # # Load in the subsetted methylation dge list
  # DGE_SubsetMethylation <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_SubsetMethylation.rds")
  # 
  # # Load in the chromosome lengths
  # load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/vitisChrLengths.RData")
  # 
  # ## Load in the map file 
  # BarossaZone <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/mapFiles/BarossaZone.rds")
  # #Soil_AWHC <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/mapFiles/Soil_AWHC.rds")
  # WineRegions <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/mapFiles/WineRegions.rds")
  # 
  # DGE_SubsetMethylation <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_SubsetMethylation.rds")
  # 
  ###############################################
  ###   Rendering of the Barossa leflet map   ###
  ###############################################
  
  output$BarossaMap <- renderLeaflet({
    BarossaMap <- leaflet() %>%
      addProviderTiles(providers$Esri.NatGeoWorldMap) %>%
      addPolygons(data = WineRegions, weight = 10, col = 'black') %>%
      #addPolygons(data = BarossaZone, weight = 5, col = 'red') %>%
      #addPolygons(data = Soil_AWHC, weight = 1) %>%
      addMarkers(lng = 139.033633, lat = -34.416577, popup = "Northern Grounds 1, Barossa Valley")%>%
      addMarkers(lng = 139.001078, lat = -34.413146, popup = "Northern Grounds 2, Barossa Valley") %>%
      addMarkers(lng = 138.976153, lat = -34.414945, popup = "Northern Grounds 3, Barossa Valley") %>%
      addMarkers(lng = 139.041833, lat = -34.427955, popup = "Northern Grounds 4, Barossa Valley") %>%
      
      addMarkers(lng = 139.008384, lat = -34.475305, popup = "Central Grounds 1, Barossa Valley") %>%
      addMarkers(lng = 138.995848, lat = -34.475355, popup = "Central Grounds 2, Barossa Valley") %>%
      addMarkers(lng = 138.946878, lat = -34.482674, popup = "Central Grounds 3, Barossa Valley") %>%
      addMarkers(lng = 138.991688, lat = -34.485512, popup = "Central Grounds 4, Barossa Valley") %>%
      
      addMarkers(lng = 138.951285, lat = -34.556889, popup = "Eastern Edge 1, Barossa Valley") %>%
      addMarkers(lng = 138.993547, lat = -34.521639, popup = "Eastern Edge 2, Barossa Valley") %>%
      addMarkers(lng = 138.991335, lat = -34.515250, popup = "Eastern Edge 3, Barossa Valley") %>%
      addMarkers(lng = 138.982400, lat = -34.546472, popup = "Eastern Edge 4, Barossa Valley") %>%
      
      addMarkers(lng = 138.933188, lat = -34.572833, popup = "Southern Grounds 1, Barossa Valley") %>%
      addMarkers(lng = 138.890905, lat = -34.624119, popup = "Southern Grounds 2, Barossa Valley") %>%
      addMarkers(lng = 138.918583, lat = -34.597124, popup = "Southern Grounds 3, Barossa Valley") %>%
      
      addMarkers(lng = 138.919776, lat = -34.495654, popup = "Western Ridge 1, Barossa Valley") %>%
      addMarkers(lng = 138.939971, lat = -34.451348, popup = "Western Ridge 2, Barossa Valley") %>%
      addMarkers(lng = 138.933749, lat = -34.420239, popup = "Western Ridge 3, Barossa Valley") %>%
      addMarkers(lng = 138.908783, lat = -34.494848, popup = "Western Ridge 4, Barossa Valley") %>%
      
      addMarkers(lng = 139.108200, lat = -34.584654, popup = "Eden Valley 1") %>%
      addMarkers(lng = 139.080000, lat = -34.624119, popup = "Eden Valley 2") %>%
      addMarkers(lng = 139.101978, lat = -34.515250, popup = "Eden Valley 3")
    
    BarossaMap
  })
  
  ############################
  ###   Create MDS plots   ###
  ############################
  

  
  ## Code for plotting the MDS plots ##
  MDStoPlot <- reactive({
  #   currentVar <- input$variable
  # numericVar <- is.numeric(appmetaG[[currentVar]])
  # int <- if_else(numericVar, 1, 0)
  # designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
  #   model.matrix(data = appmetaG)
  # # Run the voom method   
  # voomGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE) 
  # 
  # # Plot the Gene expression MDS  
  # mdsGeneExpr <- plotMDS(voomGeneExpr$E, 
  #                        labels = colnames(DGE_RNAseq$counts))
  # # Write this dodgy little function  
  # all_values <- function(x) {
  #   if(is.null(x)) return(NULL)
  #   paste0(names(x), ": ", format(x), collapse = "<br />")
  # }
  # # Set the variable  
  # Variable <- appmetaG[[currentVar]]
  # # Set the region  
  # Region <- DGE_RNAseq$samples$group
  # # Set the colour palette  
  # MDScolour <- colorRampPalette(c("brown", "yellow"))(100)
  
  # This is the code in which we will make the Methylation mds plot, it is embedded in the function that
  # Makes the gene expression one, I've done this so that the required objects need only load once 
  # to lessen the waiting times
  MethylMDS <- reactive({
    
    plot_MDS(DGElist = DGE_SubsetMethylation,
             variable = input$variable,
             metadata = appmetaM)
    # #First we select our variable
    # numericVar <- is.numeric(appmetaM[[currentVar]])
    # int <- if_else(numericVar, 1, 0)
    # # then develop the design matrix
    # designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
    #   model.matrix(data = appmetaM)
    # 
    # # Next we use the voom method to generate our log counts, note a subsetted dge list is used
    # # This is so we can cut down loading times, using all the methylation sites add a minute onto
    # # Onto the loading time
    # voomMethylation <- voom(DGE_SubsetMethylation, designMatrix, plot = FALSE) 
    # 
    # # Get the mds object ready fo visualisation
    # mdsMethylation <- plotMDS(voomMethylation$E, 
    #                           labels = colnames(DGE_SubsetMethylation$counts))
    # # Now we specify the all_values function which is necessary for our tool_tips
    # all_values <- function(x) {
    #   if(is.null(x)) return(NULL)
    #   paste0(names(x), ": ", format(x), collapse = "<br />")
    # }
    # 
    # # Now we can begin with our visualisation
    # Variable <- appmetaM[[currentVar]]
    # Region <- DGE_SubsetMethylation$samples$group
    # MDScolour <- colorRampPalette(c("brown", "yellow"))(100)
    # 
    # # Visualisation of the methylation mds in ggvis
    # MDSofVariables <- mdsMethylation[[3]] %>% 
    #   set_colnames(c("Dim1", "Dim2")) %>%
    #   as.data.frame() %>%
    #   ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
    #   layer_points() %>%
    #   add_legend("fill", title = paste("Variable:", currentVar),
    #              orient = "left") %>%
    #   add_legend("shape", title = "Regions",
    #              orient = "right")
    # 
    # 
    # MDSofVariables %>% add_tooltip(all_values, "hover")
  })
  
  MethylMDS %>% bind_shiny("mMDS", "mp_ui")
  
  # MDSofVariables <- mdsGeneExpr[[3]] %>% 
  #   set_colnames(c("Dim1", "Dim2")) %>%
  #   as.data.frame() %>%
  #   ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
  #   layer_points() %>%
  #   add_legend("fill", title = paste("Variable:", currentVar),
  #              orient = "left") %>%
  #   add_legend("shape", title = "Regions",
  #              orient = "right")
  # 
  # 
  # MDSofVariables %>% add_tooltip(all_values, "hover") %>% bind_shiny("MDS", "p_ui")
  
  # Plot the MDS plot for gene expression using ggvis  
  # mdsGeneExpr[[3]] %>%
  #   set_colnames(c("Dim1", "Dim2")) %>%
  #   as.data.frame() %>%
  #   ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
  #   layer_points() %>%
  #   add_legend("fill", title = paste("Variable:", currentVar),
  #              orient = "left") %>%
  #   add_legend("shape", title = "Regions",
  #              orient = "right") %>%
  #   add_tooltip(all_values, "hover")
  
  plot_MDS(DGElist = DGE_RNAseq,
           variable = input$variable,
           metadata = appmetaG)
  
  })
  
  # Bind it to shiny  
  MDStoPlot %>% bind_shiny("MDS", "p_ui")
  
  
  #MDSofVariables <- ggplotly(MDSofVariables)
  #MDSofVariables
  
  ########################
  ###   Heatmap code   ###
  ########################
  
  ## Here is where we start the heatmap plot ##
  output$Heatmap <- renderPlotly({
    
    # # Reactive function to select the variable in the app  
    # currentVar <- input$variable1
    # # If else statement to prepare whether or not it is discreet or categorical  
    # numericVar <- is.numeric(appmetaG[[currentVar]])
    # int <- if_else(numericVar, 1, 0)
    # 
    # # DesignMatrix of the metadata variable  
    # designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
    #   model.matrix(data = appmetaG)
    # # fitting the gene expression data to a linear fit  
    # fitGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE) %>%
    #   lmFit(design = designMatrix) %>%
    #   eBayes()
    # # number of total genes  
    # nGenesTotal <- nrow(fitGeneExpr)
    # # Define the column of the design matrix so we don't plot the slope  
    # slopeContrast <- colnames(fitGeneExpr$design)[2]
    # # Select the top genes based on p-values  
    # topGenes <- topTable(fitGeneExpr, coef = slopeContrast, number = nGenesTotal) %>%
    #   rownames_to_column("GeneID") %>%
    #   as_data_frame()
    # # Selecting the topGenes for the Pvalue bar  
    # geneExprPValFull <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)
    # # Determining the sample order  
    # topGenesTotal <- rownames(topTable(fitGeneExpr, coef = slopeContrast, number = nGenesTotal))
    # inputValsTotal <- appmetaG[[currentVar]]
    # sampleOrder <- order(inputValsTotal)
    # # Create the gene Expression dataset to plot  
    # # geneExprCPM2plot <- cpm(DGE_RNAseq, log = TRUE)[topGenesTotal, sampleOrder]
    # # colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)
    # # # Edit the column names  
    # # colnames(geneExprCPM2plot) <- colnames(geneExprCPM2plot) %>% gsub(pattern = "_Kaesler", replacement = "8") %>%
    # #   gsub(pattern = "_Elder_Nur", replacement = "6") %>%
    # #   gsub(pattern = "_Nuriootpa", replacement = "5") %>%
    # #   gsub(pattern = "_Hens_Keyn", replacement = "22") %>%
    # #   gsub(pattern = "_Hens_Eden", replacement = "21") %>%
    # #   gsub(pattern = "_Heathvale", replacement = "20")
    # # Scale rows function 
    # scale_rows = function(x){
    #   m = apply(x, 1, mean, na.rm = T)
    #   s = apply(x, 1, sd, na.rm = T)
    #   return((x - m) / s)
    # }
    # # Scale rows with that function 
    # #scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)
    # 
    # ## It seems like I have some repetition of code in this section, this may be why my heatmap is rendering slower ##
    # ##  and the clustering function may not be executing in some cases, will need to double check that there is no  ##
    # ##                                     unnecessary iteration of code.                                           ##
    # 
    # # Clustering of the heat map samples  
    # # dendogramObject <- as.dendrogram(hclust(dist(scaledCPMgeneExpr), method = "ward.D2"))
    # # clusterOrder <- order.dendrogram(dendogramObject)
    # # clusteredMatrix <- scaledCPMgeneExpr[clusterOrder, ]
    # # geneExprCPM2ggplot <- reshape2::melt(clusteredMatrix, value.name = "scaledCPM")
    # # colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")
    # # # Selecting the genes and the p-value from the data frame 
    # # geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)
    # # GeneOrder <- levels(geneExprCPM2ggplot$GeneID)
    # # GeneOrder <- data.frame(GeneOrder, 1)
    # # rownames(GeneOrder) <- GeneOrder$GeneOrder
    # # colnames(GeneOrder) <- c("GeneID", "boo")
    # # geneExprPVal <- left_join(GeneOrder, geneExprPVal, by = "GeneID")
    # # geneExprPVal$boo <- NULL
    # # # Left joining the data of the pvalue onto the gene names 
    # # geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPVal, by = 'GeneID')
    # # geneExprCPM2ggplotPVal$GeneID <- as.factor(geneExprCPM2ggplotPVal$GeneID)
    # 
    # # Generating top table of the genes and setting rownames to a column 
    # nGenes <- input$Slider
    # topGenes <- topTable(fitGeneExpr, coef = slopeContrast, number = nGenes)
    # topGenes <- rownames_to_column(topGenes, "GeneID")
    # inputVals <- appmetaG[[currentVar]]
    # sampleOrder <- order(inputVals)
    # # Select the genes and p-value  
    # geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)
    # 
    # # this code here may not be necessary, don't know why I have it repeated  
    # # nGenes <- input$Slider
    # topGenesNames <- rownames(topTable(fitGeneExpr, coef = slopeContrast, number = nGenes))
    # # inputVals <- appmetaG[[currentVar]]
    # # sampleOrder <- order(inputVals)
    # 
    # # Once again this code appears to be iterated  
    # geneExprCPM2plot <- cpm(DGE_RNAseq, log = TRUE)[topGenesNames, sampleOrder]
    # colnames(geneExprCPM2plot) <- make.names(colnames(geneExprCPM2plot), unique = TRUE)
    # 
    # colnames(geneExprCPM2plot) <- colnames(geneExprCPM2plot) %>% gsub(pattern = "_Kaesler", replacement = "8") %>%
    #   gsub(pattern = "_Elder_Nur", replacement = "6") %>%
    #   gsub(pattern = "_Nuriootpa", replacement = "5") %>%
    #   gsub(pattern = "_Hens_Keyn", replacement = "22") %>%
    #   gsub(pattern = "_Hens_Eden", replacement = "21") %>%
    #   gsub(pattern = "_Heathvale", replacement = "20")
    # # Seems repeated  
    # scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)
    # 
    # # Clustering of the heat map samples  
    # dendogramObject <- as.dendrogram(hclust(dist(scaledCPMgeneExpr), method = "ward.D2"))
    # clusterOrder <- order.dendrogram(dendogramObject)
    # clusteredMatrix <- scaledCPMgeneExpr[clusterOrder, ]
    # geneExprCPM2ggplot <- reshape2::melt(clusteredMatrix, value.name = "scaledCPM")
    # colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")
    # # Selecting the genes and the p-value from the data frame 
    # geneExprPVal <- dplyr::select(topGenes, .data$GeneID, .data$P.Value)
    # GeneOrder <- levels(geneExprCPM2ggplot$GeneID)
    # GeneOrder <- data.frame(GeneOrder, 1)
    # rownames(GeneOrder) <- GeneOrder$GeneOrder
    # colnames(GeneOrder) <- c("GeneID", "boo")
    # geneExprPVal <- left_join(GeneOrder, geneExprPVal, by = "GeneID")
    # geneExprPVal$boo <- NULL
    # # Left joining the data of the pvalue onto the gene names 
    # geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPVal, by = 'GeneID')
    # geneExprCPM2ggplotPVal$GeneID <- as.factor(geneExprCPM2ggplotPVal$GeneID)
    # 
    # #geneExprCPM2ggplot <- reshape2::melt(scaledCPMgeneExpr, value.name = "scaledCPM")
    # 
    # #colnames(geneExprCPM2ggplot) <- c("GeneID", "Vineyard", "scaledCPM")
    # 
    # # Yep need to delete repeated code   
    # geneExprCPM2ggplotPVal <- left_join(geneExprCPM2ggplot, geneExprPVal, by = 'GeneID')
    # 
    # # Date in preparation for the visualisation of colour bars and annotation plots  
    # nColours <- 101
    # varRange <- range(appmetaG[[currentVar]])
    # varGradient <- floor(seq(varRange[1], varRange[2], length.out = nColours))
    # varPositions <- findInterval(appmetaG[[currentVar]], varGradient)
    # varPalette <- colorRampPalette(c("brown", "yellow"))(nColours)
    # varColours <- varPalette[varPositions]
    # 
    # # This data is for the pvalue bar  
    # nColoursPVal <- nGenesTotal
    # PValrange <- range(geneExprPValFull$P.Value)
    # PValgradient <- seq(from = PValrange[1], to = PValrange[2], length.out = nColoursPVal)
    # PValpositions <- findInterval(geneExprPVal$P.Value, PValgradient)
    # PValpalette <- colorRampPalette(c("green", "light green", "yellow", "pink", "red", "dark red"))(nColoursPVal)
    # PValcolours <- PValpalette[PValpositions]
    # 
    # # Colour palette for the heatmap  
    # heatmapColours <- colorRampPalette(c("dark blue", "blue", "cyan", "white", "pink", "red", "dark red"))(540)
    # 
    # # Define object for function  
    # i <- nrow(unique(appmetaG[currentVar]))
    # # This is a function which is an expedient for the generation of a colour bar for either categorical or continuous data  
    # if(int == 1) {
    #   x <- seq(from = min(appmetaG[currentVar]), to = max(appmetaG[currentVar]), length.out = nColours)
    #   y <- 1
    # } else {
    #   x <- 1:i
    #   y <- 1
    # }
    # # Data created as an expedient for the generation of annotation  
    # xy <- data_frame(x, y)
    # 
    # a <- 1:101
    # b <- 1
    # c <- data_frame(a, b)
    # 
    # texta <- 1
    # textb <- 1
    # text <- data_frame(texta, textb)
    # # Define the variable annotation data to plot  
    # orderedVariables <- inputVals[sampleOrder]
    # yAxis <- 1
    # variableAnnotation <- data_frame(orderedVariables, yAxis)
    # # Different titles to give meaningful explanation of continuous and categorical variables  
    # title <- if(int == 1) {
    #   paste("Top", nGenes,"Gene Expression Analysis of", currentVar, ":", round(min(appmetaG[currentVar]), 2), "to",
    #         round(max(appmetaG[currentVar]), 2))
    # } else {
    #   paste("Analysis by", currentVar)
    # }
    # # Define the colour bar X axis  
    # colourBarXaxis <- levels(geneExprCPM2ggplotPVal$Vineyard)
    # colourBarXaxis <- as.factor(colourBarXaxis)
    # # Define the colour bar data frame  
    # colourBarDF <- data_frame(colourBarXaxis, orderedVariables)
    # colourBarDF$colourBarXaxis <-  factor(colourBarDF$colourBarXaxis, levels = colourBarDF$colourBarXaxis)
    # 
    # # Creation of the variable information annotation bar    
    # varBar <- ggplot(data = colourBarDF, aes(x = colourBarXaxis, y = orderedVariables, fill = orderedVariables)) + 
    #   geom_bar(stat = "identity", width = 1) + 
    #   scale_fill_gradientn(colours = varPalette) + 
    #   coord_cartesian(ylim = c(min(colourBarDF$orderedVariables),
    #                            (max(colourBarDF$orderedVariables)))) +
    #   guides(fill = FALSE) + 
    #   theme_void()
    # # This hasnt been defined as an object? Is it needed in the subplot?
    # ggplot(data = colourBarDF, aes(x = colourBarXaxis, y = currentVar, fill = orderedVariables)) +
    #   geom_bar(stat = "identity", width = 1) +
    #   guides(fill = FALSE) +
    #   theme_void()
    # # Defined the colour bar on top of the heatmap  
    # colourBar <- ggplot(data = xy, aes(
    #   x = seq(from = 0, to = 18.99, length.out = 101),
    #   y = y, fill = x)) +
    #   geom_raster(stat = "identity", position = "identity") +
    #   scale_fill_gradientn(colours = varPalette) +
    #   ggtitle(paste("Gene Expression Analysis of", currentVar, "from",
    #                 round(min(appmetaG[[currentVar]]), 2), "to",
    #                 round(max(appmetaG[[currentVar]]), 2))) +
    #   guides(fill = FALSE) +
    #   expand_limits(FALSE) +
    #   theme_void()
    # # Defined the P value bar   
    # PValbarWhole <- ggplot(data = geneExprPVal,
    #                        aes(x = currentVar, y = GeneID, fill = P.Value)) +
    #   geom_raster(stat = "identity", position = "identity") + 
    #   #scale_fill_continuous(colours = colorRampPalette(c("green", "dark green"))(30)) +
    #   scale_fill_gradientn(colours = PValcolours) +
    #   guides(fill = FALSE) +
    #   theme_void()
    # 
    # ggplot(geneExprPVal)
    # # defines the title that is placed on top of the pvalue bar 
    # PValtitle <- ggplot(data = text, aes(x = texta, y = textb)) +
    #   geom_text(label = "P Values") +
    #   theme_void()
    # 
    # # Hides unnecessary graphics  
    # hide <- list(
    #   title = "",
    #   zeroline = FALSE,
    #   showline = FALSE,
    #   showticklabels = FALSE,
    #   showgrid = FALSE
    # )
    # # Creation of the heatmap  
    # GeneExprheatmap <- ggplot(data = geneExprCPM2ggplot,
    #                           aes(x = Vineyard, y = GeneID, fill = scaledCPM)) + 
    #   geom_raster(stat = "identity", position = "identity") +
    #   xlab("Vineyard") + 
    #   scale_fill_gradientn(colours = heatmapColours) +
    #   theme(axis.text.x = element_text(angle = 90, hjust = 1)) +
    #   theme(axis.text.y = element_blank(),
    #         axis.ticks.y = element_blank())
    # 
    # GeneExprheatmap <- GeneExprheatmap %>% ggplotly()
    # 
    # # Visualisation within a subplot  
    # GeneExprheatmapSubplot <- subplot(varBar, PValtitle, GeneExprheatmap, PValbarWhole, 
    #                                   nrows = 2, heights = c(0.1, 0.9),
    #                                   widths = c(0.8, 0.1),
    #                                   shareX = TRUE, shareY = TRUE, titleX = FALSE,
    #                                   titleY = FALSE)
    # 
    # layout(GeneExprheatmapSubplot, title = title, showlegend = TRUE,
    #        width = 1200,
    #        height = 800,
    #        margin = list(l = 80,
    #                      r = 50,
    #                      b = 100,
    #                      t = 50),
    #        plot_bgcolor = "white")
    
    plotHeatmap(DGElist = DGE_RNAseq,
                variable = input$variable1,
                metadata = appmetaG,
                nGenes = input$Slider)
    
    
  })
  ## This is some practice code for the once contemplated gene info download ##
  # observeEvent(input$saveGenes, {
  #   session$sendCustomMessage(type = 'testmessage',
  #                             message = 'Here are your genes')
  # })
  # 
  # output$downloadData <- downloadHandler(
  #   filename = function() { paste0(currentVar, "_top_genes.txt") },
  #   content = function(file) {
  #     write()
  #   }
  # )
  
  ##################################################################################### 
  ###   This section describes the code used to generate the manhattan plots used   ###
  ###          to visualise the p-values of gene expression and methylation         ###
  #####################################################################################
  
  GeneExprPlot <- reactive({ 
    
  # currentVar <- input$variable2
  #
  # # Same as before we set our variable parameters
  # numericVar <- is.numeric(appmetaG[[currentVar]])
  # int <- if_else(numericVar, 1, 0)
  # 
  # # Develop the design matrix and fit the statistical model
  # designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
  #   model.matrix(data = appmetaG)
  # fitGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE) %>%
  #   lmFit(design = designMatrix) %>%
  #   eBayes()
  # nGenes <- nrow(fitGeneExpr)
  # slopeContrast <- colnames(fitGeneExpr$design)[2]
  # 
  # # Run the toptable function
  # topGenes <- topTable(fitGeneExpr, coef = slopeContrast, adjust.method = "BH", number = nGenes) %>%
  #   rownames_to_column("GeneID") %>%
  #   as_data_frame()
  # 
  # # Filter genes based on the top Genes
  # GenePositions <- filter(geneloci, `Gene stable ID` %in% topGenes$GeneID)
  # 
  # # set new vrownames
  # GPcol <- c("chr", "start", "end", "GeneID")
  # colnames(GenePositions) <- GPcol
  # 
  # # Get rid of bad names
  # GenePositions <- GenePositions[!(GenePositions$chr == "Un"),]
  # bad <- 1:19
  # bad <- sub("$", "_random", bad)
  # GenePositions <- GenePositions[!(GenePositions$chr %in% bad),]
  # GenePositions$chr <- sub('*', 'chr', GenePositions$chr) %>%
  #   as.factor()
  # 
  # ## Create seqinfo object for the GRanges ##
  # # Load in FASTA file
  # fastafile <- file.path("~/R/R projects/Grapevine_data/fasta_files/Vitis_vinifera.IGGP_12x.dna.toplevel.fa")
  # file.exists(fastafile)
  # 
  # # Get the fasta seqlengths and modify columns
  # VviniferaSeqLengths <- fasta.seqlengths(fastafile)
  # VviniferaSeq <- data.frame(VviniferaSeqLengths) %>%
  #   rownames_to_column("Genome")
  # 
  # # Separate the info
  # VviniferaSeqSeparated <- separate(VviniferaSeq, Genome, c("seqnames", "Strand", "Genome", "REF"), " ")
  # 
  # # remove unnecessary information
  # remove <- c("Un", "18_random", "13_random", "12_random", "7_random", "3_random", "17_random", "10_random",
  #             "16_random", "1_random", "9_random", "5_random", "11_random", "4_random", "2_random", "6_random",
  #             "19_random", "14_random", "15_random", "8_random")
  # VviniferaSeqSeparated <- VviniferaSeqSeparated[-which(VviniferaSeqSeparated$seqnames %in% remove), ]
  # VviniferaSeqSeparated$Genome <- sub('x.*', 'x', VviniferaSeqSeparated$Genome)
  # VviniferaSeqSeparated$Genome <- gsub("chromosome", "VitisVinifera", VviniferaSeqSeparated$Genome)
  # VviniferaSeqSeparated$Strand <- NULL
  # VviniferaSeqSeparated$REF <- NULL
  # VviniferaSeqSeparated$isCircular <- FALSE
  # VviniferaSeqSeparated$seqnames <- sub('*', 'chr', VviniferaSeqSeparated$seqnames)
  # 
  # # Generate the seqinfo object
  # VvSeqInfo <- Seqinfo(seqnames = VviniferaSeqSeparated$seqnames,
  #                      seqlengths = VviniferaSeqSeparated$VviniferaSeqLengths, 
  #                      isCircular = VviniferaSeqSeparated$isCircular,
  #                      genome = VviniferaSeqSeparated$Genome)
  # 
  # # Create GRanges from Gene Expression data
  # VitisGR <- makeGRangesFromDataFrame(GenePositions,
  #                                     ignore.strand = TRUE,
  #                                     seqnames.field = "chr",
  #                                     start.field = "start",
  #                                     end.field = "end",
  #                                     seqinfo = VvSeqInfo,
  #                                     keep.extra.columns = TRUE,
  #                                     starts.in.df.are.0based = FALSE)
  # 
  # # Specify promter regions
  # PromotersVitisGR <- promoters(VitisGR, upstream = 1000, downstream = 5)
  # #ManhattanToPlot %>% bind_shiny("GeneExpr", "g_ui")
  
  
  ## Make the methylation manhattan plot ###
  MethylationPlot <- reactive({
    # MethylSE <- SummarizedExperiment(assays = list(counts = OrderedSubsetMethCounts),
    #                                  rowData = MethylationInGenes, colData = phenoDataM)
    # 
    # # Extract cut sites
    # cutSites <- rowRanges(MethylSE)
    # 
    # # Adjust the cut sites to overlap recognition site on each strand
    # start(cutSites) <- ifelse(test = strand(cutSites) == '+',
    #                           yes = start(cutSites) - 1, no = start(cutSites) - 2)
    # end(cutSites) <- ifelse(test = strand(cutSites) == '+',
    #                         yes = end(cutSites) + 2, no = end(cutSites) + 1)
    # MethylSE <- subsetByOverlaps(MethylSE, cutSites)
    
    # correctCuts <- checkCuts(cutSites = cutSites,
    # genome = fastafile,
    # fasta = TRUE,
    # seq = "CCGG")
    
    # top <- diffMeth(se = MethylSE, cateogory = "Region",
    #                 condition1 = "Eden_Valley", condition2 = "Central_Grounds",
    #                 cpmThreshold = 1, thresholdSamples = 1)
    
    # Make a DGE_list from subsetted data
    # numericVar <- is.numeric(appmetaM[[currentVar]])
    # int <- if_else(numericVar, 1, 0)
    # designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>% model.matrix(data = appmetaM)
    # 
    # # Fit the data to liner model
    # fit <- glmFit(DGE_SubsetMethylation, designMatrix)
    # 
    # lrt <- glmLRT(fit, coef = colnames(designMatrix)[2])
    # 
    # # Get the top Tags for the data
    # lrt_top <- topTags(lrt, n = nrow(DGE_SubsetMethylation$counts), adjust.method = "BH", sort.by = "PValue")
    # 
    # MethylPval <- lrt_top$table
    # 
    # MethylPval$logFC <- NULL
    # MethylPval$logCPM <- NULL
    # MethylPval$LR <- NULL
    # 
    # MethylPval <- rownames_to_column(MethylPval)
    # 
    # # Separate the column into separate information
    # MethylPval <- separate(MethylPval,
    #                        rowname,
    #                        c("Chromosome", "Position", "Strand"), "\\:")
    # MethylPvalsep <- separate(MethylPval,
    #                           Position,
    #                           c("start", "end"), "\\-")
    # MethylPval <- as.tibble(MethylPval)
    # 
    # MethylPvalsep$Chromosome <- MethylPval$Chromosome %>% as.factor()
    # 
    # # Create methylation GRanges
    # MethylationGR <- makeGRangesFromDataFrame(MethylPvalsep,
    #                                           ignore.strand = TRUE,
    #                                           seqnames.field = "Chromosome",
    #                                           start.field = "start",
    #                                           end.field = "end",
    #                                           seqinfo = VvSeqInfo,
    #                                           keep.extra.columns = TRUE,
    #                                           starts.in.df.are.0based = FALSE)
    # 
    # # Reactive component to select chromosome
    # Chr <- input$Chromosome
    # 
    # MethylChrGR <- MethylationGR[seqnames(MethylationGR) == Chr]
    # MethylChrGR_df <- as_data_frame(MethylChrGR)
    # 
    # # Create axis information
    # methyax <- round(abs(seq(from = max(log10(MethylChrGR_df$PValue)), to =
    #                            min(log10(MethylChrGR_df$PValue)), length.out = 8)), 1)
    # 
    # #########################################################
    # ######    Find out way to use reactive function    ######
    # #########################################################
    # 
    # # observeEvent(input$button, {
    # #   cat("Displaying range between", input$range[1], "and", input$range[2], "/n")
    # # })
    # # 
    # # MethylChrGR_df <- eventReactive(input$button, {
    # #                                 subset(MethylChrGR_df, start > input$range[1] & start < input$range[2])
    # #                                 })
    # 
    # 
    # # Create methylation manhattan plot
    # MethylationPlot <- MethylChrGR_df %>%
    #   ggvis(x = ~start, y = ~-log10(PValue), fill := "darkred") %>%
    #   layer_points() %>%
    #   set_options(height = 380, width = 800) %>%
    #   add_axis("x",
    #            title = paste("Position on", Chr),
    #            orient = "bottom",
    #            properties = axis_props(
    #              labels = list(angle = 45,
    #                            align = "left",
    #                            fontSize = 10)
    designMatrix <- getDesignMatrix(variable = input$variable2,
                                    metadata = appmetaM)
    topMethylation <- getTopMethylation(DGElist = DGE_SubsetMethylation,
                                        designMatrix = designMatrix)
    plotManhattanMethylation(topMethylation = topMethylation,
                             SeqInfo = VvSeqInfo,
                             Chromosome = input$Chromosome)
               #))
    #    %>%
    #   add_axis("y",
    #            title = paste("-log10 of p value Methylation"),
    #            values = methyax,
    #            orient = "left")
    # 
    # all_values <- function(x) {
    #   if(is.null(x)) return(NULL)
    #   row <- MethylChrGR_df[MethylChrGR_df$start == x$start]
    #   paste0("Position: ", format(row), collapse = "<br />")}
    # 
    # MethylationPlot %>% layer_lines(y = -log10(1.5*10e-6), stroke := "red", strokeWidth := 4) %>% add_tooltip(all_values, "hover")
  })
  
  # Bind the plot to the shiny framework
  MethylationPlot %>%
    #add_tooltip(all_values, "hover") %>%
    bind_shiny("Methylation", "m_ui")
  
  
  # topGenes <- getTopGenes(DGE_RNAseq, "Planting_Date", appmetaG)
  # # Set the topGenes 
  # GeneExprPval <- topGenes
  # 
  # # Get rid of the useless information
  # GeneExprPval$logFC <- NULL
  # GeneExprPval$AveExpr <- NULL
  # GeneExprPval$t <- NULL
  # GeneExprPval$B <- NULL
  # 
  # # Manipulate and alter the gene information
  # colnames(GeneExprPval) <- c("GeneID", "PValue", "FDR")
  # 
  # GeneExprPositions <- filter(geneloci, `Gene stable ID` %in% GeneExprPval$GeneID)
  # 
  # colnames(GeneExprPositions) <- c("Chromosome", "start", "end", "GeneID")
  # 
  # GeneExprPositions$Chromosome <- sub("*", "chr", GeneExprPositions$Chromosome)
  # 
  # # Get ris of the useless information
  # bad <- 1:19
  # bad <- sub("$", "_random", bad)
  # bad <- sub("*", "chr", bad) %>% c("chrUn")
  # 
  # GeneExprPositions <- GeneExprPositions[!(GeneExprPositions$Chromosome %in% bad),]
  # 
  # GeneExprPval <- left_join(GeneExprPositions, GeneExprPval, by = "GeneID")
  # 
  # # Set The chromsome information as factors
  # GeneExprPval$Chromosome <- GeneExprPval$Chromosome %>% as.factor()
  # 
  # # Create a GRanges using the gene expression data
  # GeneExprGR <- makeGRangesFromDataFrame(GeneExprPval,
  #                                        ignore.strand = TRUE,
  #                                        seqnames.field = "Chromosome",
  #                                        start.field = "start",
  #                                        end.field = "end",
  #                                        seqinfo = VvSeqInfo,
  #                                        keep.extra.columns = TRUE,
  #                                        starts.in.df.are.0based = FALSE)
  # 
  # 
  # #Select Chromosome for Gene Expression
  # Chr <- input$Chromosome
  # 
  # # Set for the chromosome information
  # GeneExprChrGR <- GeneExprGR[seqnames(GeneExprGR) == Chr]
  # GeneExprChrGR_df <- as_data_frame(GeneExprChrGR)
  # 
  # # Set the y axis information
  # genyax <- round(seq(from = max(-log10(GeneExprChrGR_df$PValue)), to =
  #                       min(-log10(GeneExprChrGR_df$PValue)), length.out = 8), 2)
  # 
  # #################################################
  # ###   Find out way to use reactive function   ###
  # #################################################
  # 
  # goodboisG <- (GeneExprChrGR_df$FDR < 0.05)
  # goodboisG <- GeneExprChrGR_df[goodboisG, ]
  # G <- order(goodboisG$FDR)
  # goodboisG <- goodboisG[G, ]
  # maxG <- length(goodboisG$PValue)
  # cut_offG <- -log10(goodboisG$PValue[maxG])
  # 
  # # observeEvent(input$button, {
  # #   cat("Displaying ranges between", input$range[1], "and", input$range[2], "/n")
  # # })
  # # 
  # # GeneExprChrGR_df <- eventReactive(input$button, {
  # #                                   subset(GeneExprChrGR_df, start > input$range[1] & start < input$range[2])
  # #                                   })
  # # Create the Gene Expression manhattan plot
  # GeneExprPlot <- GeneExprChrGR_df %>%
  #   ggvis(x = ~start, y = ~-log10(PValue), fill := "darkblue") %>%
  #   layer_points() %>%
  #   set_options(height = 380, width = 800) %>%
  #   add_axis("x",
  #            title = "") %>%
  #   # add_axis("x",
  #   #          title = paste("Position on", Chr),
  #   #          orient = "bottom",
  #   #          properties = axis_props(
  #   #            labels = list(angle = 45,
  #   #                          align = "left",
  #   #                          fontSize = 10))) %>%
  #   add_axis("y",
  #            title = paste("-log10 of p value Gene Expression"),
  #            values = genyax,
  #            orient = "left")
  # 
  # all_values <- function(x) {
  #   if(is.null(x)) return(NULL)
  #   paste0(names(x), ": ", format(x), collapse = "<br />")
  # }
  # 
  # GeneExprPlot %>% layer_lines(y = cut_offG, stroke := "red", strokeWidth := 4) %>% add_tooltip(all_values, "hover")
  
  designMatrix <- getDesignMatrix(variable = input$variable2,
                                  metadata = appmetaG)
  topGenes <- getTopGenes(DGElist = DGE_RNAseq,
                          designMatrix = designMatrix)
  plotManhattanGeneExpr(topGenes = topGenes,
                        SeqInfo = VvSeqInfo,
                        Chromosome = input$Chromosome)
  
  })
  
  
  
  # Bind the plot to shiny
  GeneExprPlot %>%
    #add_tooltip(all_values, "hover") %>% 
    bind_shiny("GeneExpr", "g_ui")
  
  
  
  GOplot <- reactive({
    
    designMatrixG <- getDesignMatrix(variable = input$variable3,
                                     metadata = appmetaG)
    
    topGenes <- getTopGenes(DGElist = DGE_RNAseq,
                            designMatrix = designMatrixG)
    
    de <- getDEgenes(DGElist = DGE_RNAseq,
                     designMatrix = designMatrixG)
    
    designMatrixM <- getDesignMatrix(variable = input$variable3,
                                     metadata = appmetaM)
    
    dm <- getDMgenes(DGElist = DGE_SubsetMethylation,
                     designMatrix = designMatrixM,
                     allGenes = topGenes,
                     SeqInfo = VvSeqInfo)
    
    plotVariableGO(de, dm, ALLGO2ENS)
  })
  
  GOplot %>%
    bind_shiny("GO", "o_ui")
  
}
