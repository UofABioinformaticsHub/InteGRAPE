#' @title Integration of GO terms
#'
#' @description Retrieve a dataframe of GO IDs and terms that are enriched for both gene expression and methylation
#'
#' @details This function runs an intersect between the gene expression GO terms vs the differential methylation GO terms. Thus finding those enriched for both and integrating the information.
#'
#' @param de A character vector of differentially expressed genes, retrieved from the `getDEgenes` function.
#' @param dm A character vector of differentially methylated genes, retrieved from the `getDMgenes` function.
#' @param EG.GO A data frame containing information of ENSEMBL gene IDs mapped to GO IDs and terms
#'
#' @export

getIntegration <- function(de, dm, EG.GO) {

  # I was going to include this in an earlier idea for this function but decided against it

  # topGenes <- getAllGenes(DGElist = DGE_RNAseq, variable = variable, metadata = appmetaG)

  # de <- getDEgenes(DGElist = DGE_RNAseq, variable = variable, metadata = appmetaG)
  # dm <- getDMgenes(DGElist = DGE_SubsetMethylation, variable = variable, metadata = appmetaM,
  #                  allGenes = topGenes, SeqInfo = SeqInfo)

  GOde <- getGOenrichment(de, EG.GO = EG.GO)
  GOdm <- getGOenrichment(dm, EG.GO = EG.GO)

  sigGOde <- GOde[GOde$P.DE < 0.05, ]
  sigGOdm <- GOdm[GOdm$P.DE < 0.05, ]

  enrichment <- intersect(sigGOdm$Term, sigGOde$Term)

  enrichment <- sigGOdm[which(sigGOdm$Term %in% enrichment), ]

  enrichment <- data.frame(GO_ID = rownames(enrichment), GO_Term = enrichment$Term)

  enrichment
}
