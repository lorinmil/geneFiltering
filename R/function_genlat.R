
genLatent <- function(sigma2prime){

  latent <- MASS::mvrnorm(1, rep(0, 7), sigma2prime)

	return(latent)
}
