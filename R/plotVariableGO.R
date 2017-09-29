#' @title Plot GO terms from count data
#'
#' @description This function allows the plotting of GO terms directly from the two DGElists
#'
#' @details This function has the same output as the `plotGOterms` function, but it runs the entire analysis directly from the DGElists. The body will contain the functions used, and it may be better to run the analysis from each step as this function requires several inputs.
#'
#' @param DGE_RNAseq A DGElist output from `limma` including data from differential gene expression
#' @param DGE_SubsetMethylation A DGElist output from `limma` including data from differential gene methylation
#' @param variable The variable of interest, a character string or vector
#' @param appmetaG A dataframe with different variables on the x axis and samples on the y axis. Specifically for the gene expression data
#' @param appmetaM A dataframe with different variables on the x axis and samples on the y axis. Specifically for the gene methylation data
#' @param VvSeqInfo A SeqInfo object containing genomic information on your species of interest. Can either find one on ENSEMBL or create one from an object using `seqinfo` https://www.rdocumentation.org/packages/GenomeInfoDb/versions/1.8.3/topics/seqinfo
#' @param EG.GO A dataframe containing the information relating to ENSEMBL Gene IDs mapped to their respective GO terms


plotVariableGO <- function(DGE_RNAseq, DGE_SubsetMethylation, appmetaG, appmetaM, variable, VvSeqInfo, EG.GO) {

  source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDEgenes.R")
  source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getDMgenes.R")
  source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getGOenrichment.R")
  source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/plotGOterms.R")


  de <- getDEgenes(DGElist = DGE_RNAseq, variable = variable, metadata = appmetaG)
  topGenes <- getTopGenes(DGElist = DGE_RNAseq, variable = variable, metadata = appmetaG)
  dm <- getDMgenes(DGElist = DGE_SubsetMethylation, variable = variable, metadata = appmetaM,
                   allGenes = topGenes, SeqInfo = VvSeqInfo)
  GOde <- getGOenrichment(de = de, EG.GO = EG.GO)
  GOdm <- getGOenrichment(de = dm, EG.GO = EG.GO)
  GOdm$P.DM <- GOdm$P.DE
  GOdm$P.DE <- NULL

  plotGOterms(de = de, dm = dm, EG.GO = EG.GO)


}
