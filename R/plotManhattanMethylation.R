
plotManhattanMethylation <- function(DGElist, variable, metadata, SeqInfo, Chromosome = "chr1") {

  numericVar <- is.numeric(metadata[[variable]])
  int <- if_else(numericVar, 1, 0)
  designMatrix <- as.formula(paste("~", paste(int, variable, sep = "+"))) %>% model.matrix(data = metadata)

  # Fit the data to liner model
  fit <- glmFit(DGElist, designMatrix)

  lrt <- glmLRT(fit, coef = colnames(designMatrix)[2])

  # Get the top Tags for the data
  lrt_top <- topTags(lrt, n = nrow(DGElist$counts), adjust.method = "BH", sort.by = "PValue")

  # colnames(lrt_top$table)[1] <- "site"
  # rownames(lrt_top$table) <- NULL



  # - 
  # P-value comparison
  # We now have the top genes and methylatin sites along with their p-values as a result of our previous code
  # First we should put these two into meaningful objects so we don't forget them

  MethylPval <- lrt_top$table


  ## Sort out unnecessary columns
  # Now we will want to sort out any columns that are not needed in subsequent analyses

  MethylPval$logFC <- NULL
  MethylPval$logCPM <- NULL
  MethylPval$LR <- NULL


  ## Adjust column names
  ## Separate the methylation position info

  MethylPval <- rownames_to_column(MethylPval)

  MethylPval <- separate(MethylPval, rowname, c("Chromosome", "Position", "Strand"), "\\:")

  MethylPvalsep <- separate(MethylPval, Position, c("start", "end"), "\\-")

  MethylPval <- as.tibble(MethylPvalsep)



  MethylPvalsep$Chromosome <- MethylPval$Chromosome %>% as.factor()


  MethylationGR <- makeGRangesFromDataFrame(MethylPvalsep,
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
  MethylChrGR_df <- as_data_frame(MethylChrGR)


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
    ggvis(x = ~start, y = ~-log10(PValue), fill := "darkred") %>%
    layer_points() %>%
    add_axis("x",
             title = paste("Position on", Chr),
             orient = "bottom") %>%
    add_axis("y",
             title = paste("-log10 of p value"),
             values = methyax,
             orient = "left")
  manhattanPlotDM %>%
    layer_lines(y = -log10(1.5*10e-6), stroke := "red") %>%
    add_tooltip(all_valuesM, "hover")

  manhattanPlotDM 
}

