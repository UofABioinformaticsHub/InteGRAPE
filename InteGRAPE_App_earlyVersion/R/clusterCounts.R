#' @title Cluster gene counts
#' 
#' @description An easy way to cluster gene counts for enhanced visualisation within a heatmap
#' 
#' @details This function uses the `hclust` algorithm, specifically the `ward.D2` method by default, other methods include
#' @details `ward.D`, `single`, `complete`, `average` which is UPGMA, `mcquitty` which is WPGMC, and finally `centroid` which is UPGMC
#' 
#' @param counts a matrix of discrete gene counts
#' @param method the clustering method to be used, see details for the selection
#' 
#' @return a clustered matrix of original counts matrix
#' 
#' @export
#' 

clusterCounts <- function(counts, method = "ward.D2") {
  
  scale_rows = function(x){
    m = apply(x, 1, mean, na.rm = T)
    s = apply(x, 1, sd, na.rm = T)
    return((x - m) / s)
  }
  
  geneExprCPM2plot <- edgeR::cpm(counts)
  scaledCPMgeneExpr <- scale_rows(geneExprCPM2plot)
  
  dendogramObject <- as.dendrogram(hclust(dist(scaledCPMgeneExpr), method = "ward.D2"))
  clusterOrder <- order.dendrogram(dendogramObject)
  clusteredMatrix <- scaledCPMgeneExpr[clusterOrder, ]
  
  clusteredMatrix
}