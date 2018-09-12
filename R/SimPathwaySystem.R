#Function to simulate a pathway system

simPathwaySystem <- function(
  alpha
  ,Rsq
  ,nsample
  ,rind
  ,nx
  ,nxPrime
  ,ns
  ,ntranscriptNoise
  ,ng
  ,ngPrime
  ,nh
  ,ngeneNoise
  ,noiseSD
){

  sim <- gen_eff_obs_v2 (alpha            = alpha,
                         Rsq              = Rsq,
                         nsample          = nsample,
                         rind             = rind,
                         nx               = nx,
                         nxPrime          = nxPrime,
                         ns               = ns,
                         ntranscriptNoise = ntranscriptNoise,
                         ng               = ng,
                         ngPrime          = ngPrime,
                         nh               = nh,
                         ngeneNoise       = ngeneNoise,
                         noiseSD          = noiseSD
  )

  system <- sim$Obs

  x <- data.frame(system$transcripts)
  colnames(x) <- paste0('x',1:ncol(x))
  g <- data.frame(system$genes)
  colnames(g) <- paste0('g',1:ncol(g))
  y <- data.frame(system$y)
  indicatorX <- system$truetranscripts
  indicatorX[indicatorX==1] <- 'pathway'
  indicatorX[indicatorX==0] <- 'not in pathway'
  indicatorG <- system$truegenes
  indicatorG[indicatorG==1] <- 'pathway'
  indicatorG[indicatorG==0] <- 'not in pathway'


  #------------Return results-------------------------
  returnVal <- list(x, g, y, indicatorX, indicatorG)

  names(returnVal) <- c('x', 'g', 'y', 'indicatorX', 'indicatorG')

  return(returnVal)

}
