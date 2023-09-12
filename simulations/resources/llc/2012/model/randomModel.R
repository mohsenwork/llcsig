randomModel <- function(n,Bdensity=0.2,Cdensity=0.15) {
  #Creates a random model in canonical form (zero mean, no self loops).
  #INPUT:
  # n - number of observed variables
  # Bdensity - how many direct effects should be non-zero
  # Cdensity - how many covariances should be non-zero
  #OUTPUT:
  #  B  - Direct effects matrix.
  #  Ce - Covariance matrix of the disturbances.  
  #  G  - Directed graph corresponding to causal structure in B.
  #  Ge - Undirected graph corresponding to latent confounder structure in Ce.
  #WARNING: Models produced by this code for n > 20 might be unrealistic and/or unstable-
  #        Adjusting the densities helps.


  # Number of variables
  #n <- 4 # testing!

  # Random structure (density 0 < density < 1)

  B <- 1*(matrix(runif(n*n),n,n) < Bdensity)

  # Make sure diagonal is zero (we specify all models in this way
  # to take out unnecessary degrees of underdetermination)
  diag(B) <- 0

  # Random coefficients
  #old coefficient generation by Patrik
  #B <- B*matrix(rnorm(n*n),n,n)

  #saving the graph
  G<-B
  #new with a lower bound
  B <- B*matrix(runif(n*n,0.1,0.8),n,n)*matrix(sample(c(-1,1),n*n,replace=TRUE),n,n)

  #drawing the structure for the covariance matrix
  Adist <- 1*(matrix(runif(n*n),n,n) < Cdensity/2)
  diag(Adist) <- 1
  Ge<-Adist%*%t(Adist)

  Adist <- Adist*matrix(rnorm(n*n),n,n)
  Ce <- Adist %*% t(Adist)        

  # Return both
  res <- list()
  res$B <- B
  res$Ce <- Ce
  res$G <- G
  res$Ge <- Ge
  res
}

