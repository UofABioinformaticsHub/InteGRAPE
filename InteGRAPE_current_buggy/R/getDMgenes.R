#' @title Retrieve differentially methylated genes
#'
#' @description Retrieve a character vector with a list of differentially methylated genes (Should have ENSEMBL IDs)
#'
#' @details This function works in a similar way to that of the `getAllGenes` function in that its fitting a linear model to the differential methylation data.
#' @details This variation of the workflow uses a Generalised linear fit model and runs a `topTags` function on the genes.
#'
#' @return A character vector of differentially expressed genes.
#'
#' @param DGElist A DGElist output from `limma`
#' @param designMatrix A model matrix that is the output of the `getDesignMatrix` function
#' @param allGenes An output with the `topTable` from the `getAllGenes` function
#' @param SeqInfo A SeqInfo object containing genomic information on your species of interest. Can either find one on ENSEMBL or create one from an object using `seqinfo` https://www.rdocumentation.org/packages/GenomeInfoDb/versions/1.8.3/topics/seqinfo
#'
#' @export
getDMgenes <- function(DGElist, designMatrix, allGenes, SeqInfo) {

  require(magrittr)
  require(readr)
  require(tibble)
  require(dplyr)
  require(limma)
  require(edgeR)
  require(R.utils)
  require(tidyr)
  require(GenomicRanges)

  if (is.null(DGElist)) {
    stop("DGElist is missing, this must contain the methylation data. Please see the package edgeR for more information")
  }
  if (is.null(variable)) {
    stop("Please provide a variable to assess. Make sure it is the same variable as the one used for your allGenes")
  }
  if (is.null(metadata)) {
    stop("Please provide the metadata. This should contain information on your variable and be in a data frame format")
  }
  if (is.null(allGenes)) {
    stop("In order to find differentially methylated genes, we need a GRanges object in which to overlap our methylation positions with
          which contains the gene ranges")}
  if (is.null(SeqInfo)) {
    stop("It is very important that an SeqInfo object characterising your organism is loaded")
  }


  # First get the DE genes GRanges which we'll need later
  GeneExprPval <- allGenes

  GenePositions <- filter(geneloci, `Gene stable ID` %in% allGenes$GeneID)
  GPcol <- c("chr", "start", "end", "GeneID")
  colnames(GenePositions) <- GPcol

  GeneExprPositions <- filter(geneloci, `Gene stable ID` %in% GeneExprPval$GeneID)
  colnames(GeneExprPositions) <- c("Chromosome", "start", "end", "GeneID")
  GeneExprPositions$Chromosome <- sub("*", "chr", GeneExprPositions$Chromosome)

  bad <- 1:19
  bad <- sub("$", "_random", bad)
  bad <- sub("*", "chr", bad) %>% c("chrUn")

  GeneExprPositions <- GeneExprPositions[!(GeneExprPositions$Chromosome %in% bad),]
  GeneExprPval <- dplyr::left_join(GeneExprPositions, GeneExprPval, by = "GeneID")
  GeneExprPval$Chromosome <- GeneExprPval$Chromosome %>% as.factor()

  GeneExprGR <- GenomicRanges::makeGRangesFromDataFrame(GeneExprPval,
                                         ignore.strand = TRUE,
                                         seqnames.field = "Chromosome",
                                         start.field = "start",
                                         end.field = "end",
                                         seqinfo = SeqInfo,
                                         keep.extra.columns = TRUE,
                                         starts.in.df.are.0based = FALSE)



  # Fit the data to liner model
  fit <- edgeR::glmFit(DGElist, designMatrix)

  lrt <- edgeR::glmLRT(fit, coef = colnames(designMatrix)[2])

  # Get the top Tags for the data
  lrt_top <- edgeR::topTags(lrt, n = nrow(DGElist$counts), adjust.method = "BH", sort.by = "PValue")

  # Check how many statistically significant sites there are
  summary(lrt_top$table$PValue <= 0.05)

  # First define how we'll subset it
  dm <- 0.05 >= lrt_top$table$PValue

  # Then subset it, easy peesy

  getDiffMeth <- function(dm, top) {
    dmGenes <- lrt_top[dm , ]

    dmGenes <- dmGenes$table

    dmGenes <- tibble::rownames_to_column(dmGenes)

    dmGenes <- tidyr::separate(dmGenes,
                        rowname,
                        c("Chromosome", "Position", "Strand"), "\\:")
    dmGenes <- tidyr::separate(dmGenes,
                        Position,
                        c("start", "end"), "\\-")
    dmGenes <- tibble::as.tibble(dmGenes)

    dmGenes$Chromosome <- dmGenes$Chromosome %>% as.factor()

    dmGenes
  }

  dmGenes <- getDiffMeth(dm = dm, top = lrt_top)

  MethylationGR <- GenomicRanges::makeGRangesFromDataFrame(dmGenes,
                                            ignore.strand = TRUE,
                                            seqnames.field = "Chromosome",
                                            start.field = "start",
                                            end.field = "end",
                                            seqinfo = SeqInfo,
                                            keep.extra.columns = TRUE,
                                            starts.in.df.are.0based = FALSE)
  PromotersVitisGR <- GenomicRanges::promoters(VitisGR, upstream = 5000, downstream = 5000)

  dmGenesGR <- IRanges::subsetByOverlaps(PromotersVitisGR, MethylationGR)
  dmGenes <- tibble::as_data_frame(dmGenesGR)

  dm <- dmGenes$GeneID
  dm
}
