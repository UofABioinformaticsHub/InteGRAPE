#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(edgeR)
library(limma)
library(ggvis)
library(magrittr)
library(dplyr)
library(tidyr)
library(ggthemes)
library(readr)
library(tibble)
library(RColorBrewer)
library(DESeq)
library(splitstackshape)
library(rtracklayer)
library(plotly)
library(SummarizedExperiment)
library(IRanges)
library(GenomicRanges)
library(GenomicFeatures)
library(leaflet)
library(rgdal)
library(raster)
library(Biostrings)
library(shiny)
library(shinycssloaders)


source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getTopGenes.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDEgenes.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotHeatmap.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDMgenes.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getGOenrichment.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotGOterms.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotVariableGO.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getIntegration.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDesignMatrix.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotManhattamGeneExpr.R")
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotManhattanMethylation.R")


# Define UI for application that draws a histogram
ui <- navbarPage(title = "Inte-GRAPE-tion",
                 theme = "bootstrap.css",
                 tabPanel(title = "Geography",
                          tags$h1("Locations of each of the different samples"),
                          tags$h2("This map shows the locations of all the different samples within the Barossa wine growing regions."),
                          withSpinner(leafletOutput("BarossaMap",
                                                    width = "80%",
                                                    height = 800),
                                      type = getOption("spinner.type", default = 3),
                                      color = getOption("spinner.color", default = "#0275D8"),
                                      color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
                 tabPanel(title = "MDS Plots",
                          tags$h1("Multi-dimensional scaling plots for each sample"),
                          tags$h2("These plots visualise the separation of each of the variables across the samples"),
                          sidebarPanel(
                            selectInput("variable",
                                        label = "Variable of Interest",
                                        choices = list("Planting_Date" = c("Planting_Date", "Site", "Clone", "Pruning"),
                                                       "Organisation" = c("Space_rows", "Space_vines", "Bud_number", "Orientation"),
                                                       "Midrow_Management" = c("Herbicide_Mid", "Mechanical_Mid", "VGrowth_Mid",
                                                                               "Mulch_Mid", "Mid_Management"),
                                                       "Undervine_Management" = c("Herbicide_Under", "Mechanical_Under",
                                                                                  "VGrowth_Under", "Mulch_Under", "Under_Management"),
                                                       "Canopy_Management" = c("Canopy_Manage"),
                                                       "Irrigation" = c("Irrigation"),
                                                       "Subregion" = c("Subregion"),
                                                       "Soil" = c("Soil", "pH", "P", "Salinity", "Moisture", "Water_hold_cap"),
                                                       "Elevation" = c("Elevation"),
                                                       "Climate" = c("Mean_rain", "Seas_rain", "Jan_temp",
                                                                     "Web_record", "Grow_temp", "Grow_days"),
                                                       "HarvestDate" = c("Harvest_Date"),
                                                       "Berry" = c("50Berry_weight", "TSS_berry", "pH_berry",
                                                                   "TA_berry", "E_conc", "Col_berry", "Col_berry_weight",
                                                                   "Phenolics_Ber", "Phenolics_Ber_weight"),
                                                       "Wine_Chemistry" = c("Citric", "Gluconic", "Malic", "Pyruvic",
                                                                            "Tartaric", "Alcohol_wine", "pH_Wine", "TA_Wine",
                                                                            "Tannins", "Phenolics_wine", "Age_chem",
                                                                            "antho_ionisation", "total_antho"),
                                                       "Wine_Sensory" = c("Colour_dens", "Colour_dens_SO2", "Hue",
                                                                          "SO2_resis", "Colour", "Fruit_sweet", "Acidity",
                                                                          "Red_fruit", "Dark_fruit", "Spice", "Choc_Mocha",
                                                                          "Bitter", "Astringency", "Tannin_texture", "Mouthfeel")
                                        ))
                            ),
                          ggvisOutput("MDS"),
                          uiOutput("p_ui"),
                          fluidRow(column(4),
                                   column(8,
                          ggvisOutput("mMDS"),
                          uiOutput("mp_ui")
                          ))),
                 tabPanel(title = "Gene Expression Analysis",
                          tags$h1("Heat-map visualisation of gene expression"),
                          tags$h2("The gene expression levels within the grapevine were fitted to a linear model and visualised within a heat-map"),
                          sidebarLayout(position = "left",
                          sidebarPanel(
                            selectInput("variable1",
                                        label = "Variable of Interest",
                                        choices = list("Planting_Date" = c("Planting_Date", "Site", "Clone", "Pruning"),
                                                       "Organisation" = c("Space_rows", "Space_vines", "Bud_number", "Orientation"),
                                                       "Midrow_Management" = c("Herbicide_Mid", "Mechanical_Mid", "VGrowth_Mid",
                                                                               "Mulch_Mid", "Mid_Management"),
                                                       "Undervine_Management" = c("Herbicide_Under", "Mechanical_Under",
                                                                                  "VGrowth_Under", "Mulch_Under", "Under_Management"),
                                                       "Canopy_Management" = c("Canopy_Manage"),
                                                       "Irrigation" = c("Irrigation"),
                                                       "Subregion" = c("Subregion"),
                                                       "Soil" = c("Soil", "pH", "P", "Salinity", "Moisture", "Water_hold_cap"),
                                                       "Elevation" = c("Elevation"),
                                                       "Climate" = c("Mean_rain", "Seas_rain", "Jan_temp",
                                                                     "Web_record", "Grow_temp", "Grow_days"),
                                                       "HarvestDate" = c("Harvest_Date"),
                                                       "Berry" = c("50Berry_weight", "TSS_berry", "pH_berry",
                                                                   "TA_berry", "E_conc", "Col_berry", "Col_berry_weight",
                                                                   "Phenolics_Ber", "Phenolics_Ber_weight"),
                                                       "Wine_Chemistry" = c("Citric", "Gluconic", "Malic", "Pyruvic",
                                                                            "Tartaric", "Alcohol_wine", "pH_Wine", "TA_Wine",
                                                                            "Tannins", "Phenolics_wine", "Age_chem",
                                                                            "antho_ionisation", "total_antho"),
                                                       "Wine_Sensory" = c("Colour_dens", "Colour_dens_SO2", "Hue",
                                                                          "SO2_resis", "Colour", "Fruit_sweet", "Acidity",
                                                                          "Red_fruit", "Dark_fruit", "Spice", "Choc_Mocha",
                                                                          "Bitter", "Astringency", "Tannin_texture", "Mouthfeel")
                                                       )),
                            sliderInput("Slider",
                                        label = "Number of Genes",
                                        min = 2,
                                        max = 100,
                                        value = 30)),
                           mainPanel(
                             fluidRow(
                               column(8, 
                          withSpinner(plotlyOutput("Heatmap", height = "800px"),
                                      type = getOption("spinner.type", default = 3),
                                      color = getOption("spinner.color", default = "#0275D8"),
                                      color.background = getOption("spinner.color.background", default = "#FFFFFF"))),
                          column(4,
                                 plotOutput(outputId = "LinearTrend")
                          ))))),
                          
                 tabPanel(title = "Gene Expression vs Methylation",
                          tags$h1("Comparison of p-values from gene expression and methylation"),
                          tags$h2("Through the comparison of p-values, genes which are both differentially methylated and differentially expressed can be identified"),
                          sidebarPanel(
                            selectInput("variable2",
                                        label = "Variable of Interest",
                                        choices = list("Planting_Date" = c("Planting_Date", "Site", "Clone", "Pruning"),
                                                       "Organisation" = c("Space_rows", "Space_vines", "Bud_number", "Orientation"),
                                                       "Midrow_Management" = c("Herbicide_Mid", "Mechanical_Mid", "VGrowth_Mid",
                                                                               "Mulch_Mid", "Mid_Management"),
                                                       "Undervine_Management" = c("Herbicide_Under", "Mechanical_Under",
                                                                                  "VGrowth_Under", "Mulch_Under", "Under_Management"),
                                                       "Canopy_Management" = c("Canopy_Manage"),
                                                       "Irrigation" = c("Irrigation"),
                                                       "Subregion" = c("Subregion"),
                                                       "Soil" = c("Soil", "pH", "P", "Salinity", "Moisture", "Water_hold_cap"),
                                                       "Elevation" = c("Elevation"),
                                                       "Climate" = c("Mean_rain", "Seas_rain", "Jan_temp",
                                                                     "Web_record", "Grow_temp", "Grow_days"),
                                                       "HarvestDate" = c("Harvest_Date"),
                                                       "Berry" = c("50Berry_weight", "TSS_berry", "pH_berry",
                                                                   "TA_berry", "E_conc", "Col_berry", "Col_berry_weight",
                                                                   "Phenolics_Ber", "Phenolics_Ber_weight"),
                                                       "Wine_Chemistry" = c("Citric", "Gluconic", "Malic", "Pyruvic",
                                                                            "Tartaric", "Alcohol_wine", "pH_Wine", "TA_Wine",
                                                                            "Tannins", "Phenolics_wine", "Age_chem",
                                                                            "antho_ionisation", "total_antho"),
                                                       "Wine_Sensory" = c("Colour_dens", "Colour_dens_SO2", "Hue",
                                                                          "SO2_resis", "Colour", "Fruit_sweet", "Acidity",
                                                                          "Red_fruit", "Dark_fruit", "Spice", "Choc_Mocha",
                                                                          "Bitter", "Astringency", "Tannin_texture", "Mouthfeel")
                                        )),
                            
                            # Call seqinfo object here
                            selectInput("Chromosome",
                                        label = "Select Chromosome",
                                        choices = paste0("chr", 1:19)),
                            uiOutput("range"),
                            sliderInput("range",
                                         label = "Select range of bases",
                                         min = 0,
                                         max = 34000000,
                                         value = c(1, 20000000),
                                         animate = TRUE
                                         )
                          ), 
                          ggvisOutput("GeneExpr"),
                          uiOutput("g_ui"),
                          fluidRow(column(4),
                                   column(8,
                          ggvisOutput("Methylation"),
                          uiOutput("m_ui")))),
                 
                 tabPanel(title = "GO term correlations",
                          tags$h1("Identifying GO terms from genes that are both differentially expressed and methylated"),
                          tags$h2(""),
                          sidebarPanel(
                            selectInput("variable3",
                                        label = "Variable of interest",
                                        choices = list("Planting_Date" = c("Planting_Date", "Site", "Clone", "Pruning"),
                                                       "Organisation" = c("Space_rows", "Space_vines", "Bud_number", "Orientation"),
                                                       "Midrow_Management" = c("Herbicide_Mid", "Mechanical_Mid", "VGrowth_Mid",
                                                                               "Mulch_Mid", "Mid_Management"),
                                                       "Undervine_Management" = c("Herbicide_Under", "Mechanical_Under",
                                                                                  "VGrowth_Under", "Mulch_Under", "Under_Management"),
                                                       "Canopy_Management" = c("Canopy_Manage"),
                                                       "Irrigation" = c("Irrigation"),
                                                       "Subregion" = c("Subregion"),
                                                       "Soil" = c("Soil", "pH", "P", "Salinity", "Moisture", "Water_hold_cap"),
                                                       "Elevation" = c("Elevation"),
                                                       "Climate" = c("Mean_rain", "Seas_rain", "Jan_temp",
                                                                     "Web_record", "Grow_temp", "Grow_days"),
                                                       "HarvestDate" = c("Harvest_Date"),
                                                       "Berry" = c("50Berry_weight", "TSS_berry", "pH_berry",
                                                                   "TA_berry", "E_conc", "Col_berry", "Col_berry_weight",
                                                                   "Phenolics_Ber", "Phenolics_Ber_weight"),
                                                       "Wine_Chemistry" = c("Citric", "Gluconic", "Malic", "Pyruvic",
                                                                            "Tartaric", "Alcohol_wine", "pH_Wine", "TA_Wine",
                                                                            "Tannins", "Phenolics_wine", "Age_chem",
                                                                            "antho_ionisation", "total_antho"),
                                                       "Wine_Sensory" = c("Colour_dens", "Colour_dens_SO2", "Hue",
                                                                          "SO2_resis", "Colour", "Fruit_sweet", "Acidity",
                                                                          "Red_fruit", "Dark_fruit", "Spice", "Choc_Mocha",
                                                                          "Bitter", "Astringency", "Tannin_texture", "Mouthfeel")
                                        ))),
                            ggvisOutput("GO"),
                            uiOutput("o_ui")
                            )
                 
# 
#                  tabPanel(title = "GO Enrichment Analysis",
#                           tags$h1("Table of enriched GO terms"),
#                           tags$h2("Table showing terms that are enriched for both differential gene expression and methylation, each GO term is linked to an information page"),
#                           sidebarPanel(
#                             selectInput("variable4",
#                                        label = "Variable of interest",
#                                        choices = list("Planting_Date" = c("Planting_Date", "Site", "Clone", "Pruning"),
#                                                                 "Organisation" = c("Space_rows", "Space_vines", "Bud_number", "Orientation"),
#                                                                 "Midrow_Management" = c("Herbicide_Mid", "Mechanical_Mid", "VGrowth_Mid",
#                                                                                         "Mulch_Mid", "Mid_Management"),
#                                                                 "Undervine_Management" = c("Herbicide_Under", "Mechanical_Under",
#                                                                                            "VGrowth_Under", "Mulch_Under", "Under_Management"),
#                                                                 "Canopy_Management" = c("Canopy_Manage"),
#                                                                 "Irrigation" = c("Irrigation"),
#                                                                 "Subregion" = c("Subregion"),
#                                                                 "Soil" = c("Soil", "pH", "P", "Salinity", "Moisture", "Water_hold_cap"),
#                                                                 "Elevation" = c("Elevation"),
#                                                                 "Climate" = c("Mean_rain", "Seas_rain", "Jan_temp",
#                                                                               "Web_record", "Grow_temp", "Grow_days"),
#                                                                 "HarvestDate" = c("Harvest_Date"),
#                                                                 "Berry" = c("50Berry_weight", "TSS_berry", "pH_berry",
#                                                                             "TA_berry", "E_conc", "Col_berry", "Col_berry_weight",
#                                                                             "Phenolics_Ber", "Phenolics_Ber_weight"),
#                                                                 "Wine_Chemistry" = c("Citric", "Gluconic", "Malic", "Pyruvic",
#                                                                                      "Tartaric", "Alcohol_wine", "pH_Wine", "TA_Wine",
#                                                                                      "Tannins", "Phenolics_wine", "Age_chem",
#                                                                                      "antho_ionisation", "total_antho"),
#                                                                 "Wine_Sensory" = c("Colour_dens", "Colour_dens_SO2", "Hue",
#                                                                                    "SO2_resis", "Colour", "Fruit_sweet", "Acidity",
#                                                                                    "Red_fruit", "Dark_fruit", "Spice", "Choc_Mocha",
#                                                                                    "Bitter", "Astringency", "Tannin_texture", "Mouthfeel")
#                                        )),
#                             sliderInput("number",
#                                         label = "Number of terms",
#                                         min = 1,
#                                         max = 50,
#                                         value = 30)),
#                           tableOutput("GOenrich"))
                 )

                 

# Define server logic required to draw a histogram
server <- function(input, output, session) {
  
  # output$range <- renderUI({
  #   
  #   load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/vitisChrLengths.RData")
  #   
  #   sliderInput("range",
  #               label = "Select range of bases",
  #               min = 0,
  #               max = vitisChrLengths[input$Chromosome, ],
  #               value = c(1, 1000000)
  #               )
  # })
  
  # all_values <- function(x) {
  #   if(is.null(x)) return(NULL)
  #   paste0(names(x), ": ", format(x), collapse = "<br />")
  # }
  
source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/GrapevineApp/Global.R")
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
  MDStoPlot <- reactive({currentVar <- input$variable
  numericVar <- is.numeric(appmetaG[[currentVar]])
  int <- if_else(numericVar, 1, 0)
  designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
    model.matrix(data = appmetaG)
  # Run the voom method   
  voomGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE) 
  
  # Plot the Gene expression MDS  
  mdsGeneExpr <- plotMDS(voomGeneExpr$E, 
                         labels = colnames(DGE_RNAseq$counts))
  # Write this dodgy little function  
  all_values <- function(x) {
    if(is.null(x)) return(NULL)
    paste0(names(x), ": ", format(x), collapse = "<br />")
  }
  # Set the variable  
  Variable <- appmetaG[[currentVar]]
  # Set the region  
  Region <- DGE_RNAseq$samples$group
  # Set the colour palette  
  MDScolour <- colorRampPalette(c("brown", "yellow"))(100)
  
  # This is the code in which we will make the Methylation mds plot, it is embedded in the function that
  # Makes the gene expression one, I've done this so that the required objects need only load once 
  # to lessen the waiting times
  MethylMDS <- reactive({
    #First we select our variable
    numericVar <- is.numeric(appmetaM[[currentVar]])
    int <- if_else(numericVar, 1, 0)
    # then develop the design matrix
    designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
      model.matrix(data = appmetaM)
    
    # Next we use the voom method to generate our log counts, note a subsetted dge list is used
    # This is so we can cut down loading times, using all the methylation sites add a minute onto
    # Onto the loading time
    voomMethylation <- voom(DGE_SubsetMethylation, designMatrix, plot = FALSE) 
    
    # Get the mds object ready fo visualisation
    mdsMethylation <- plotMDS(voomMethylation$E, 
                              labels = colnames(DGE_SubsetMethylation$counts))
    # Now we specify the all_values function which is necessary for our tool_tips
    all_values <- function(x) {
      if(is.null(x)) return(NULL)
      paste0(names(x), ": ", format(x), collapse = "<br />")
    }
    
    # Now we can begin with our visualisation
    Variable <- appmetaM[[currentVar]]
    Region <- DGE_SubsetMethylation$samples$group
    MDScolour <- colorRampPalette(c("brown", "yellow"))(100)
    
    # Visualisation of the methylation mds in ggvis
    MDSofVariables <- mdsMethylation[[3]] %>% 
      set_colnames(c("Dim1", "Dim2")) %>%
      as.data.frame() %>%
      ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
      layer_points() %>%
      add_legend("fill", title = paste("Variable:", currentVar),
                 orient = "left") %>%
      add_legend("shape", title = "Regions",
                 orient = "right")
    
    
    MDSofVariables %>% add_tooltip(all_values, "hover")
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
  mdsGeneExpr[[3]] %>%
    set_colnames(c("Dim1", "Dim2")) %>%
    as.data.frame() %>%
    ggvis(x = ~Dim1, y = ~Dim2, fill = ~Variable, shape = ~Region) %>%
    layer_points() %>%
    add_legend("fill", title = paste("Variable:", currentVar),
               orient = "left") %>%
    add_legend("shape", title = "Regions",
               orient = "right") %>%
    add_tooltip(all_values, "hover")
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
    plotHeatmap(DGElist = DGE_RNAseq, variable = input$variable1, metadata = appmetaG, nGenes = input$Slider)
  })
  
  output$LinearTrend <- renderPlot({
    eventdata <- event_data("plotly_click")
    
    #validate(need(!is.null(eventdata), "Hover over the time series chart to populate this heatmap"))
    designMatrix <- getDesignMatrix(input$variable3, appmetaG)

     datapoint <- eventdata$pointNumber
     datapoint <- unlist(datapoint)[[1]]
    Gene <- rownames(clusteredMatrix)[datapoint]
    
     #window <- as.numeric(input$window)
     
     #datapoint
    
     # rng <- (datapoint - window):(datapoint + window)
     # cormat <- round(cor(stockdata[rng, 1:5]), 2)
    
      voomObj <- limma::voom(DGE_RNAseq, designMatrix, plot = FALSE)
      #selectGene <- geneExprCPM2ggplot[geneExprCPM2ggplot$scaledCPM == datapoint, ]
     
      #Gene <- selectGene$GeneID
     
      #GeneP <- topGenes[topGenes$GeneID == Gene, ]
      #pvalue <- GeneP$P.Value
      GeneCounts <- voomObj$E[Gene, ]
     
      LinearRep <- tibble::tibble(factor(names(GeneCounts[sampleOrder]), levels =
                                           unique(names(GeneCounts[sampleOrder]))), GeneCounts[sampleOrder])
     LinearRep$Variable <- inputValsTotal[sampleOrderTotal]
      colnames(LinearRep) <- c("Region", "Counts", "Variable")
     
     
      plotLinearRep <- ggplot2::ggplot(data = LinearRep, aes(x = Variable, y = Counts)) +
        ggplot2::geom_point(stat = "identity")
      plotLinearRep + ggplot2::geom_smooth(method = "lm") #+ ggtitle(paste(Gene, ": p =", round(pvalue, 11))) + ggplot2::theme_bw()

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
  
  GeneExprPlot <- reactive({ currentVar <- input$variable2
  
  # Same as before we set our variable parameters
  numericVar <- is.numeric(appmetaG[[currentVar]])
  int <- if_else(numericVar, 1, 0)
  
  # Develop the design matrix and fit the statistical model
  designMatrix <- as.formula(paste("~", paste(int, currentVar, sep = "+"))) %>%
    model.matrix(data = appmetaG)
  fitGeneExpr <- voom(DGE_RNAseq, designMatrix, plot = FALSE) %>%
    lmFit(design = designMatrix) %>%
    eBayes()
  nGenes <- nrow(fitGeneExpr)
  slopeContrast <- colnames(fitGeneExpr$design)[2]
  
  # Run the toptable function
  topGenes <- topTable(fitGeneExpr, coef = slopeContrast, adjust.method = "BH", number = nGenes) %>%
    rownames_to_column("GeneID") %>%
    as_data_frame()
  
  # Filter genes based on the top Genes
  GenePositions <- filter(geneloci, `Gene stable ID` %in% topGenes$GeneID)
  
  # set new vrownames
  GPcol <- c("chr", "start", "end", "GeneID")
  colnames(GenePositions) <- GPcol
  
  # Get rid of bad names
  GenePositions <- GenePositions[!(GenePositions$chr == "Un"),]
  bad <- 1:19
  bad <- sub("$", "_random", bad)
  GenePositions <- GenePositions[!(GenePositions$chr %in% bad),]
  GenePositions$chr <- sub('*', 'chr', GenePositions$chr) %>%
    as.factor()
  
  ## Create seqinfo object for the GRanges ##
  # Load in FASTA file
  fastafile <- file.path("~/R/R projects/Grapevine_data/fasta_files/Vitis_vinifera.IGGP_12x.dna.toplevel.fa")
  file.exists(fastafile)
  
  # Get the fasta seqlengths and modify columns
  VviniferaSeqLengths <- fasta.seqlengths(fastafile)
  VviniferaSeq <- data.frame(VviniferaSeqLengths) %>%
    rownames_to_column("Genome")
  
  # Separate the info
  VviniferaSeqSeparated <- separate(VviniferaSeq, Genome, c("seqnames", "Strand", "Genome", "REF"), " ")
  
  # remove unnecessary information
  remove <- c("Un", "18_random", "13_random", "12_random", "7_random", "3_random", "17_random", "10_random",
              "16_random", "1_random", "9_random", "5_random", "11_random", "4_random", "2_random", "6_random",
              "19_random", "14_random", "15_random", "8_random")
  VviniferaSeqSeparated <- VviniferaSeqSeparated[-which(VviniferaSeqSeparated$seqnames %in% remove), ]
  VviniferaSeqSeparated$Genome <- sub('x.*', 'x', VviniferaSeqSeparated$Genome)
  VviniferaSeqSeparated$Genome <- gsub("chromosome", "VitisVinifera", VviniferaSeqSeparated$Genome)
  VviniferaSeqSeparated$Strand <- NULL
  VviniferaSeqSeparated$REF <- NULL
  VviniferaSeqSeparated$isCircular <- FALSE
  VviniferaSeqSeparated$seqnames <- sub('*', 'chr', VviniferaSeqSeparated$seqnames)
  
  # Generate the seqinfo object
  VvSeqInfo <- Seqinfo(seqnames = VviniferaSeqSeparated$seqnames,
                       seqlengths = VviniferaSeqSeparated$VviniferaSeqLengths, 
                       isCircular = VviniferaSeqSeparated$isCircular,
                       genome = VviniferaSeqSeparated$Genome)
  
  # Create GRanges from Gene Expression data
  VitisGR <- makeGRangesFromDataFrame(GenePositions,
                                      ignore.strand = TRUE,
                                      seqnames.field = "chr",
                                      start.field = "start",
                                      end.field = "end",
                                      seqinfo = VvSeqInfo,
                                      keep.extra.columns = TRUE,
                                      starts.in.df.are.0based = FALSE)
  
  # Specify promter regions
  PromotersVitisGR <- promoters(VitisGR, upstream = 1000, downstream = 5)
  #ManhattanToPlot %>% bind_shiny("GeneExpr", "g_ui")
  
  
  ## Make the methylation manhattan plot ###
  MethylationPlot <- reactive({
    plotManhattanMethylation(DGElist = DGE_SubsetMethylation, variable = input$variable2, metadata = appmetaM,
                             range_start = input$range[1], range_end = input$range[2],
                             SeqInfo = VvSeqInfo, Chromosome = input$Chromosome)
  })
  
  # Bind the plot to the shiny framework
  MethylationPlot %>%
    #add_tooltip(all_values, "hover") %>%
    bind_shiny("Methylation", "m_ui")
  
  GeneExprPlot <- plotManhattanGeneExpr(topGenes, VvSeqInfo, Chromosome = input$Chromosome,
                                        range_start = input$range[1], range_end = input$range[2])
  
  })
  

  # Bind the plot to shiny
  GeneExprPlot %>%
    #add_tooltip(all_values, "hover") %>% 
    bind_shiny("GeneExpr", "g_ui")
  


GOplot <- reactive({
  plotVariableGO(DGE_RNAseq, DGE_SubsetMethylation, appmetaG, appmetaM, input$variable3, VvSeqInfo, EG.GO)
})

GOplot %>%
  bind_shiny("GO", "o_ui")

# enrichtable <- reactive({
#   
#   currentVar <- input$variable3
#   
#   topGenes <- getTopGenes(DGElist = DGE_RNAseq, variable = currentVar, metadata = appmetaG)
#   
#   de <- getDEgenes(DGElist = DGE_RNAseq, variable = currentVar, metadata = appmetaG)
#   
#   dm <- getDMgenes(DGElist = DGE_SubsetMethylation, variable = currentVar, metadata = appmetaM,
#              allGenes = topGenes, SeqInfo = VvSeqInfo)
#   
#   enrichment <- getIntegration(de = de, dm = dm, EG.GO = ALLGO2ENS)
#   
#   enrichment %>% head(input$number)
# })
# 
# output$GOenrich <- renderDataTable(enrichtable)

}

shinyApp(ui = ui, server = server)
