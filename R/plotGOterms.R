#' @title Plot significant GO terms
#'
#' @description This is a good check to see if we can expect any functionally enriched ontologies before running a formal analysis
#'
#' @details Plotting the GO terms associated with differentially expressed genes and differentially methylated genes will give an overview on what to expect when running further GO analyses such as a Fishers Exact Test or Hypergeometric Distribution Model.
#'
#' @return A plot showing which GO terms are likely to be enriched
#'
#' @param de A character vector of differentially expressed genes, an output from the `getDEgenes` function
#' @param dm A character vector of differentially methylated genes, an output from the `getDMgenes` function
#' @param EG.GO A dataframe containing the information relating to ENSEMBL Gene IDs mapped to their respective GO terms
#'
#' @export
plotGOterms <- function(de, dm, EG.GO) {

  library(biomaRt)
  library(magrittr)
  library(readr)
  library(tibble)
  library(dplyr)
  library(limma)
  library(edgeR)
  library(R.utils)
  library(GO.db)
  library(tidyr)
  library(splitstackshape)

  source("~/R/R projects/Grapevine_data/Grapevine_Data_Analysis/R_files/R_code/Scripts/getGOenrichment.R")


  GOde <- getGOenrichment(de = de, EG.GO = EG.GO)
  GOdm <- getGOenrichment(de = dm, EG.GO = EG.GO)

  GOdm$P.DM <- GOdm$P.DE
  GOdm$P.DE <- NULL

GO_Combined <- dplyr::left_join(GOdm, GOde, by = c("Term", "Ont", "N"))

GO_Combined$id <- 1:nrow(GO_Combined)

all_valuesGO <- function(x) {
  if(is.null(x)) return(NULL)
  GO <- GO_Combined[GO_Combined$id == x$id,]
  paste0("GO Term", ": ", GO$Term, collapse = "<br/ >")
}

GOplot <- GO_Combined %>%
  ggvis::ggvis(x = ~-log10(P.DM), y = ~-log10(P.DE), key := ~id) %>%
  ggvis::layer_points(fill := "blue")

GOplot %>%
  ggvis::layer_lines(x = ~-log10(0.05), stroke := "red", strokeWidth := 2) %>%
  ggvis::layer_lines(y = ~-log10(0.05), stroke := "red", strokeWidth := 2) %>%
  ggvis::add_tooltip(all_valuesGO, "click")
}

#plotGOterms(GOdm = GOdm, GOde = GOde)
