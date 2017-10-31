#' @title Run a GO enrichment analysis
#'
#' @description This function retrieves functionally enriched GO terms
#'
#' @details This function uses a hypergeometric distribution model to perform the functional enrichment analysis. Essentially
#' @details this model finds the probability of retrieving two genes from the entire `universe` with the same GO term. It finds
#' @details the probabilities of finding two of the same, one with and one without, or both without the GO term. From this it
#' @details assigns a p value to tells whats the probability of viewing this functional enrichment if the null hypothesis is true.
#'
#' @return A dataframe containing information on the Gene ID, the GO term and ontology, and the assigned p value from the hypergeometric distribution model (`phyper`).
#'
#' @param de list of differentially expressed or differentially methylated genes
#' @param EG.GO A dataframe containing the information relating to ENSEMBL Gene IDs mapped to their respective GO terms
#'
#' @export
getGOenrichment <- function(de, EG.GO) {

  require(biomaRt)
  require(magrittr)
  require(readr)
  require(tibble)
  require(dplyr)
  require(limma)
  require(edgeR)
  require(R.utils)
  require(GO.db)
  require(tidyr)
  require(splitstackshape)

  colnames(EG.GO) <- c("gene_id", "go_id", "Ontology")

  # Find all duplicated rows
  d <- duplicated(EG.GO[, c("gene_id", "go_id", "Ontology")])

  # Use subsetting to remove all of the duplicated rows
  EG.GO <- EG.GO[!d, ]

  # Create the universe of genes by finding the unique genes
  universe <- unique(EG.GO$gene_id)

  # Then change them to characters
  universe <- as.character(universe)

  # This if statement is from goana, it checks to see whether the universe genes are actually there
  Total <- length(unique(EG.GO$gene_id))
  if (Total <1L) {stop("No genes found in universe")}

  # Checks to see if the de object is a list, as it needs to be this class for downstream analyses
  if (!is.list(de)) {
    de <- list(DE = de)
  }

  # This is a check to make sure all of the components of the de are vectors
  if (!all(vapply(de, is.vector, TRUE))) {
    stop("components of de should be vectors")
  }

  # Now the de object gets changed to a character vector
  de <- lapply(de, as.character)

  # Here we specify how many sets of genes we have here, we should have 1
  nsets <- length(de)

  # Now we do some groovey stuff with the names of the de object
  names(de) <- limma::trimWhiteSpace(names(de))
  NAME <- names(de)
  i <- which(NAME == "" | is.na(NAME))
  NAME[i] <- paste0("DE", i)
  names(de) <- limma::makeUnique(NAME)
  #universe <- universe[!dup]

  # This finds the total number of the genes
  Total <- length(unique(EG.GO$gene_id))

  # We've done this check twice now, do we need it here?
  if (Total < 1L) {
    stop("No genes found in universe")
  }

  # Now we want to find the differentially expressed genes within the universe style object
  isDE <- lapply(de, function(x) EG.GO$gene_id %in% x)

  # Then take the length of the subset
  TotalDE <- lapply(isDE, function(x) length(unique(EG.GO$gene_id[x])))

  # This simply finds the number of differentially expressed bois
  nDE <- length(isDE)

  # not sure why we use `do.call()` here but this function created a dataframe of the
  # de genes as well as another column `N`. The `N` column will serve the purpose of recording the
  # total number of genes associated with a GO term.
  X <- do.call(cbind, c(N = 1, isDE))

  # This step is in preparation for the `rowsum` function. This will add together all rows which have the
  # same GO term. The beauty of this is that as a result, the `N` column will then tell us exactly how
  # many genes we have within that set
  group <- paste(EG.GO$go_id, EG.GO$Ontology, sep = ".")

  # here is where we execute this command, so S is the total number of
  # different GO terms within the EG.GO data. But the pvalue will tell us
  # the probability of the hypergeometric model which indicates the enrichment
  # of specific GO terms
  S <- rowsum(X, group = group, reorder = FALSE)

  # We prepare this matrix to run the for loop which will retrieve all of the p-values for us
  P <- matrix(0, nrow = nrow(S), ncol = nsets)

  # Now we run the hypergeometric distribution model. Essentially this step is measuring the probability
  # in which two genes that are pulled out of the data at random will possess the GO term of interest.
  # If both have it, then the results is considered TRUE TRUE. So this model will assess the possibilities
  # of getting TRUE TRUE, TRUE FALSE, and FALSE FALSE. From there it determines the enrichment of the GO terms
  # and assigns a pvalue telling us the probability of the statistical model result was found by chance
  # or if the models interpretations are statistically significant. I could be completely wrong about this
  # as well, but a great source is http://blog.nextgenetics.net/?e=94
  # I think this method is likened to bootstrapping

  for (j in 1:nsets) P[, j] <- phyper(q = S[, 1 + j] - 0.5,
                                      m = TotalDE[[j]],
                                      n = Total - TotalDE[[j]],
                                      k = S[, "N"],
                                      lower.tail = FALSE)

  # Now we wish to put all our data back together to see which terms are enriched where
  # first separate the names into matrix columns
  g <- strsplit2(rownames(S), split = "\\.")

  #Then we want to link up the GO terms
  TERM <- suppressMessages(AnnotationDbi::select(GO.db::GO.db,
                                                 keys = g[, 1], columns = "TERM"))

  # Finally save it all into a dataframe
  Results <- data.frame(Term = TERM[[2]], Ont = g[, 2], S,
                        P, stringsAsFactors = FALSE)

  # Change the rownames to the gene ids
  rownames(Results) <- g[, 1]

  # Sort out the column names
  colnames(Results)[3 + nsets + (1L:nsets)] <- paste0("P.", names(de))

  # Now view the fruits of your efforts
  Results }
