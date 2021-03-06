% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filterReads.R
\name{filterReads}
\alias{filterReads}
\title{Filter reads based by set threshold}
\usage{
filterReads(DGElist, reads, samples, plot = FALSE)
}
\arguments{
\item{DGElist}{An object of class `DGElist`, typically an output from the `limma` package}

\item{reads}{The set minimum number of reads that needs to be within a sample}

\item{samples}{The set minimum number of samples that need to have the specified minimum number of reads within them.}
}
\description{
This function subsets the count data from a DGElist by a threshold specified by the user
}
\details{
The threshold set within this function will subset reads that do not meet the threshold of x amount of reads in x amount of samples. If the user sets `reads` to `1` and `samples` to `3`, then this function will remove those genes that have less than 3 samples with any reads in them. Hence, if only two samples had reads in them, they will be removed.
}
