filterPathwayGenes <- function(
  x
  ,g
  ,y
  ,numClusters_1
  ,numClusters_2
){

  #---------------------------------------------------------------------------------
  #Step 1: Filter X and G Correlations with Y, individually-------------------------
  #---------------------------------------------------------------------------------

  #----------Calculate absolute correlations----------
  corXY <- apply(x, 2, function(z) abs(cor(z, y)))
  corGY <- apply(g, 2, function(z) abs(cor(z, y)))

  #----------Perform k-means on correlations----------
  clustersXY <- kmeans(corXY, centers=numClusters_1)
  clustersGY <- kmeans(corGY, centers=numClusters_1)

  #----------Remove Smallest Cluster------------------
  #Filter the X items
  clusterMeansX_1 <- clustersXY$centers
  smallestX <- which(clusterMeansX_1 == min(clusterMeansX_1))
  removeX_1 <- which(clustersXY$cluster == smallestX)
  step1X <- x

  #Filter the G items
  clusterMeansG_1 <- clustersGY$centers
  smallestG <- which(clusterMeansG_1 == min(clusterMeansG_1))
  removeG_1 <- which(clustersGY$cluster == smallestG)
  step1G <- g

  #---------------------------------------------------------------------------------
  #Step 2: Filter X and G Correlations with eachother-------------------------------
  #---------------------------------------------------------------------------------

  #------------Correlation between X and G-------------
  corXG <- abs(cor(step1X, step1G))
  xSums <- rowSums(corXG)
  gSums <- colSums(corXG)

  #Kmeans on the correlation sums
  clustersXSums <- kmeans(xSums, centers=numClusters_2)
  clustersGSums <- kmeans(gSums, centers=numClusters_2)

  #----------Remove Smallest Cluster------------------
  #Filter the X items
  clusterMeansX_2 <- clustersXSums$centers
  smallestX <- which(clusterMeansX_2 == min(clusterMeansX_2))
  removeX_2 <- which(clustersXSums$cluster == smallestX)
  removeX_both <- removeX_1[which(removeX_1 %in% removeX_2)]
  step2X <- x
  if(length(removeX_both) > 0){
    step2X <- step2X[,-removeX_both]
  }

  #Filter the G items
  clusterMeansG_2 <- clustersGSums$centers
  smallestG <- which(clusterMeansG_2 == min(clusterMeansG_2))
  removeG_2 <- which(clustersGSums$cluster == smallestG)
  removeG_both <- removeG_1[which(removeG_1 %in% removeG_2)]
  step2G <- g
  if(length(removeG_both) > 0){
    step2G <- step2G[,-removeG_both]
  }

  #------------Return results-------------------------
  returnVal <- list(corXY, corGY, xSums, gSums, removeX_1, removeG_1
                    , removeX_2, removeG_2, removeX_both, removeG_both)
  names(returnVal) <- c('corXY', 'corGY', 'xSums', 'gSums'
                        , 'removeX_1', 'removeG_1', 'removeX_2', 'removeG_2'
                        , 'removeX_both', 'removeG_both')

  return(returnVal)

}
