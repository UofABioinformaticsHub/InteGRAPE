#' @title Filter reads based by set threshold
#' 
#' @description This function subsets the count data from a DGElist by a threshold specified by the user
#' 
#' @details The threshold set within this function will subset reads that do not meet the threshold of x amount of reads in x amount of samples. If the user sets `reads` to `1` and `samples` to `3`, then this function will remove those genes that have less than 3 samples with any reads in them. Hence, if only two samples had reads in them, they will be removed.
#' 
#' @param DGElist An object of class `DGElist`, typically an output from the `limma` package
#' @param reads The set minimum number of reads that needs to be within a sample
#' @param samples The set minimum number of samples that need to have the specified minimum number of reads within them.
#' 
#' @export
#' 
filterReads <- function(DGElist, reads, samples, plot = FALSE) {
  # take the counts per million of the counts in the DGElist
  cpm <- edgeR::cpm(DGElist)
  lcpm <- edgeR::cpm(DGElist, log = TRUE)
  
  # Set the threshold, this code will specify the threshold of the minimum number of reads in a minimum amount of samples
  keep.exprs <- rowSums(cpm > reads) >= samples
  
  # Now execute the subsetting function
  DGElist <- DGElist[
    keep.exprs,
    keep.lib.sizes = FALSE
    ]
  return(DGElist)
  
  if(plot == TRUE) {
    nsamples <- ncol(DGElist)
    
    col <- brewer.pal(nsamples, "Paired")
    
    par(mfrow = c(1,2))
    
    plot(
      density(lcpm[,1]),
      col = col[1],
      lwd = 2,
      ylim = c(0, 0.21),
      las = 2,
      main = "",
      xlab = "")
    
    title(main = "A. Raw data", 
          xlab = "Log-cpm")
    
    abline(v = 0,
           lty = 3)
    
    for (i in 2:nsamples) {
      den <- density(
        lcpm[,i])
      
      lines(den$x,
            den$y,
            col = col[i],
            lwd = 2)
    }
    
    legend(
      "topright",
      samplenames,
      text.col = col,
      bty = "n")
    
    lcpm <- edgeR::cpm(
      DGElist,
      log = TRUE)
    
    plot(
      density(lcpm[,1]),
      col = col[1],
      lwd = 2,
      ylim = c(0, 0.21),
      las = 2,
      main = "",
      xlab = "")
    
    title(main = "B. Filtered data",
          xlab = "Log-cpm")
    
    abline(v = 0,
           lty = 3)
    
    for (i in 2:nsamples) {
      den <- density(
        lcpm[,1])
      
      lines(den$x, den$y,
            col = col[i],
            lwd = 2)}
    
    legend("topright",
           samplenames,
           text.col = col,
           bty = "n")
  } else {break}
}

