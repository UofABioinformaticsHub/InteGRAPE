library(edgeR)
library(limma)
library(ggvis)
library(magrittr)
library(dplyr)
library(tidyr)
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



###########################################################
###   Prepare all the data required for the shiny app   ###
###########################################################

# Load in the Methylation data
# DGE_Methylation <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_Methylation.rds")

# Load in the Methylation data
# DGE_Methylation <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_Methylation.rds")

# This is where you should keep a file path to your DGElist to be used in the app, this one here
# for me is where I've stored my RNA-seq DGE-list, so any transcriptomic data should be put here.
# Typically it should be RNA-seq count data (integers) as microarray data yields continuous values.
# The limma methods are designed to allow RNA-seq data to be analysed using methods designed for
# micro-arrays, but I haven't tested the workflow with microarray data yet.
load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_RNAseq.RData")


# This is where you should load in your meta-data. Because of the difference in samples between the 
# RNA-seq and msgbs data, I had to make two separate meta-data objects that link up to each of the 
# different data sets. I'm not sure if you'll have this problem, but it's work adressing here.
# There is a way to set keys for the table but thats a future step for me I guess.
load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/appmetaG.RData")


# This is the matrix of gene counts clustered with the `hclust` algorithm. Not sure why it's here but it
# may be important. 
clusteredMatrix <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/clusteredMatrix.rds")


# This is where you can load in your DGE-list for the methylation data from the msgbs, the one I had was
# subsetted to include only those methylation sites which were located within genes.
DGE_SubsetMethylation <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/DGE_SubsetMethylation.rds")


# Once again, the meta-data object is specific for each DGElist so as to match up with the dimensions
# of the DGElist count data. 
load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/appmetaM.RData")


# This is where you can load in the counts for the msgbs data, I need to make sure I revise all of the methods
# in the app to see if this step is needed, but for the sake of running the app, I'll have this
# included.
OrderedSubsetMethCounts <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/OrderedSubsetMethCounts.rds")


# This file will contain those methyltion sites that occur within genes. I'm still not sure
# if this file is strictly needed but it's good to be on the safe side. You can easily generate
# this kind of file by simply creating a gene ranges object for the RNA-seq and msgbs data and
# then subsetting by the methylation sites that overlap with genes. The function is called `subsetByOverlaps`
# within the IRanges package.
MethylationInGenes <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/MethylationInGenes.rds")


# Load in the methylation phenotypic data, this is just the annotations for the methylation (msgbs) count 
# data. Essentially this just contains the informatio you'd need to create the DGElist object.
phenoDataM <- read_csv(file = "~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/phenoDataM.csv")



               ####################################
               ######                        ######
               #####     Annotation files     #####
               ######                        ######
               ####################################



# This file here is used to create the GRanges of all genes in vitis vinifera which can then be
# subsetted by the gene name in order to generate a GRanges of all of your genes, and then subsequently
# find all methylation sites within genes
gffFile <- file.path("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/Annotation_Files/Vitis_vinifera_annotation.gff")
file.exists(gffFile)
genesGR <- import.gff(gffFile) %>%
  subset(type == "CDS") %>%
  split(f = .$CDS)

# This does the exact same thing as above. If one of the annotation files doesn't work, it's always 
# good to have a spare on hand so you don't need to go back to ENSEMBL plants to download another file
# just to find out that ENSEMBL is down and then you have to go home and drink wine until it's back
# up.
GR <- file.path("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/Annotation_Files/Vitis_vinifera.IGGP_12x.36.gff3.gz")
file.exists(gffFile)
genesGR <- import(gffFile)

# This file contains information on the location of all of the genes, I'm pretty sure this is used for 
# something really important.
TSVfile <- file.path("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/Annotation_Files/mart_export.txt")
file.exists(TSVfile)
geneloci <- read_tsv(TSVfile)


# I think having the chromosome legths for vitis vinifera was important when generating the SeqInfo object
# but it looks like this info might not be needed anymore...
load("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/vitisChrLengths.RData")

## This is just used to draw out the map
BarossaZone <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/mapFiles/BarossaZone.rds")
#Soil_AWHC <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/mapFiles/Soil_AWHC.rds")
WineRegions <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/mapFiles/WineRegions.rds")


# You'll need this universal set of genes mapped back to GO terms for your GO term over-representation analysis
# in order for the function to work. This is the EG.GO argument
ALLGO2ENS <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/ALLGO2ENS.rds")


# SeqInfo object required to annotate the GRanges object
VvSeqInfo <- readRDS("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/dataFiles/VvSeqInfo.rds")

