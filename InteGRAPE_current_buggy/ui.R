
# This is the user-interface definition of a Shiny web application.
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
# library(leaflet)
# library(rgdal)
# library(raster)
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
                                          withSpinner(plotlyOutput("Heatmap"),
                                                      type = getOption("spinner.type", default = 3),
                                                      color = getOption("spinner.color", default = "#0275D8"),
                                                      color.background = getOption("spinner.color.background", default = "#FFFFFF"))
                                        ))),
                 
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
                            uiOutput("range")
                            # sliderInput("range",
                            #             label = "Select range of bases",
                            #             min = 0,
                            #             max = vitisChrLengths[input$Chromosome, ],
                            #             value = c(1, 20000000),
                            #             animate = TRUE
                            #             ),
                            #actionButton("button", "Adjust Range")
                          ), 
                          ggvisOutput("GeneExpr"),
                          uiOutput("g_ui"),
                          fluidRow(column(4),
                                   column(8,
                                          ggvisOutput("Methylation"),
                                          uiOutput("m_ui")))),
                 
                 tabPanel(title = "GO analysis",
                          tags$h1("Identifying GO terms from genes that are both differentially expressed and methylated"),
                          tags$h2(""),
                          sidebarPanel(
                            selectInput("variable3",
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
                                        ))),
                          ggvisOutput("GO"),
                          uiOutput("o_ui")
                 )
)
