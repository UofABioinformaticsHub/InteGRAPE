#' @title Plot gene positions with differential expression
#'
#' @description This function plots the results of a differential gene expression analysis within a Manhattan Plot.
#'
#' @details Plotting the information regarding postions of differentially expressed genes allows a whole genome view of gene expression and also makes possible a comparison by visualisation of sites with differential methylation.
#'
#' @param topGenes A dataframe output from the `getTopGenes` function
#' @param SeqInfo An object of class `Seqinfo` containing information on the chromosomes within an organism
#' @param Chromosome A character string containing the chromosome to plot, ie "chr1", "chr4", "chr13", etc.
#'
#' @export

plotManhattanGeneExpr <- function(topGenes, SeqInfo, Chromosome = "chr1") {
GeneExprPval <- topGenes

## Sort out unnecessary columns
#Now we will want to sort out any columns that are not needed in subsequent analyses
GeneExprPval$logFC <- NULL
GeneExprPval$AveExpr <- NULL
GeneExprPval$t <- NULL
GeneExprPval$B <- NULL

colnames(GeneExprPval) <- c("GeneID", "PValue", "FDR")

GeneExprPositions <- filter(geneloci, `Gene stable ID` %in% GeneExprPval$GeneID)

colnames(GeneExprPositions) <- c("Chromosome", "start", "end", "GeneID")

GeneExprPositions$Chromosome <- sub("*", "chr", GeneExprPositions$Chromosome)

bad <- 1:19
bad <- sub("$", "_random", bad)
bad <- sub("*", "chr", bad) %>% c("chrUn")


GeneExprPositions <- GeneExprPositions[!(GeneExprPositions$Chromosome %in% bad),]

GeneExprPval <- left_join(GeneExprPositions, GeneExprPval, by = "GeneID")

GeneExprPval$Chromosome <- GeneExprPval$Chromosome %>% as.factor()



GeneExprGR <- makeGRangesFromDataFrame(GeneExprPval,
                                       ignore.strand = TRUE,
                                       seqnames.field = "Chromosome",
                                       start.field = "start",
                                       end.field = "end",
                                       seqinfo = SeqInfo,
                                       keep.extra.columns = TRUE,
                                       starts.in.df.are.0based = FALSE)


## Single Chromosomes
### Select Chromosome

#Select Chromosome for Gene Expression
Chr <- Chromosome

GeneExprChrGR <- GeneExprGR[seqnames(GeneExprGR) == Chr]
GeneExprChrGR_df <- as_data_frame(GeneExprChrGR)

yaxis <- round(abs(seq(from = max(-log10(GeneExprChrGR_df$PValue)), to = min(log10(MethylChrGR_df$PValue)), length.out = 10)), 0)

genyax <- round(abs(seq(from = max(log10(GeneExprChrGR_df$PValue)), to =
                          min(log10(GeneExprChrGR_df$PValue)), length.out = 8)), 1)

## Find FDR cutof

goodboisG <- (GeneExprChrGR_df$FDR < 0.05)
goodboisG <- GeneExprChrGR_df[goodboisG, ]
G <- order(goodboisG$FDR)
goodboisG <- goodboisG[G, ]
maxG <- length(goodboisG$PValue)
cut_offG <- -log10(goodboisG$PValue[maxG])

GeneExprChrGR_df$id <- 1:nrow(GeneExprChrGR_df)

all_valuesG <- function(x) {
  if(is.null(x)) return(NULL)
  rowG <- GeneExprChrGR_df[GeneExprChrGR_df$id == x$id,]
  paste0("Gene ID", ": ", rowG$GeneID, "<br>",
         "P.value", ": ", rowG$PValue, collapse = "<br />")
}

manhattanPlotDGE <- GeneExprChrGR_df %>%
  ggvis(x = ~start, y = ~-log10(PValue), fill := "darkblue", key := ~id) %>%
  layer_points() %>%
  add_axis("x",
           title = paste("Position on", Chr),
           orient = "bottom") %>%
  add_axis("y",
           title = paste("-log10 of p value"),
           values = genyax,
           orient = "left")
manhattanPlotDGE %>%
  layer_lines(y = cut_offG, stroke := "red") %>%
  add_tooltip(all_valuesG, "hover")

manhattanPlotDGE

}
