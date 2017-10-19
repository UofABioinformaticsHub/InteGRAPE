#' @title Plot Methylation positions
#'
#' @description This function plots the results of a differential methylation analysis within a Manhattan Plot.
#'
#' @details Plotting the information regarding methylation positions allows a whole genome view of methylation and also makes possible a comparison by visualisation of sites for differential gene expression.
#'
#' @param topMethylation A dataframe contianing information on the significance of differential methylation sites, an output from the `getTopMethylation` function
#' @param SeqInfo An object of class `Seqinfo` containing information on the chromosomes within an organism
#' @param Chromosome A character string containing the chromosome to plot, ie "chr1", "chr4", "chr13", etc.
#'
#' @export

plotManhattanMethylation <- function(topMethylation, SeqInfo, Chromosome = "chr1") {

  # -
  # P-value comparison
  # We now have the top genes and methylatin sites along with their p-values as a result of our previous code
  # First we should put these two into meaningful objects so we don't forget them

  MethylPval <- topMethylation$table


  ## Sort out unnecessary columns
  # Now we will want to sort out any columns that are not needed in subsequent analyses

  MethylPval$logFC <- NULL
  MethylPval$logCPM <- NULL
  MethylPval$LR <- NULL


  ## Adjust column names
  ## Separate the methylation position info

  MethylPval <- tibble::rownames_to_column(MethylPval)

  MethylPval <- tidyr::separate(MethylPval, rowname, c("Chromosome", "Position", "Strand"), "\\:")

  MethylPvalsep <- tidyr::separate(MethylPval, Position, c("start", "end"), "\\-")

  MethylPval <- tibble::as.tibble(MethylPvalsep)



  MethylPvalsep$Chromosome <- MethylPval$Chromosome %>% as.factor()


  MethylationGR <- GenomicRanges::makeGRangesFromDataFrame(MethylPvalsep,
                                            ignore.strand = TRUE,
                                            seqnames.field = "Chromosome",
                                            start.field = "start",
                                            end.field = "end",
                                            seqinfo = SeqInfo,
                                            keep.extra.columns = TRUE,
                                            starts.in.df.are.0based = FALSE)


  ## Single Chromosomes
  ### Select Chromosome

  # Select Chromosome for Methylation
  Chr <- Chromosome

  MethylChrGR <- MethylationGR[seqnames(MethylationGR) == Chr]
  MethylChrGR_df <- tibble::as_data_frame(MethylChrGR)


  ## Get y axis ready

  yaxis <- round(abs(seq(from = max(-log10(GeneExprChrGR_df$PValue)), to = min(log10(MethylChrGR_df$PValue)), length.out = 10)), 0)

  methyax <- round(abs(seq(from = max(log10(MethylChrGR_df$PValue)),
                           to = min(log10(MethylChrGR_df$PValue)),
                           length.out = 8)), 1)


  MethylChrGR_df$id <- 1:nrow(MethylChrGR_df)
  MethylChrGR_df$start <- as.character(MethylChrGR_df$start)
  MethylChrGR_df$PValue <- as.numeric(MethylChrGR_df$PValue)


  all_valuesM <- function(x) {
    if(is.null(x)) return(NULL)
    rowM <- MethylChrGR_df[MethylChrGR_df$id == x$id,]
    paste0("Position", ": ", rowM$start, "<br>",
           "P.value", ": ", rowM$PValue, collapse = "<br />")
    }


  manhattanPlotDM <- MethylChrGR_df %>%
    ggvis::ggvis(x = ~start, y = ~-log10(PValue), fill := "darkred") %>%
    ggvis::layer_points() %>%
    ggvis::add_axis("x",
             title = paste("Position on", Chr),
             orient = "bottom") %>%
    ggvis::add_axis("y",
             title = paste("-log10 of p value"),
             values = methyax,
             orient = "left")
  manhattanPlotDM %>%
    ggvis::layer_lines(y = -log10(1.5*10e-6), stroke := "red") %>%
    ggvis::add_tooltip(all_valuesM, "hover")

  manhattanPlotDM
}

