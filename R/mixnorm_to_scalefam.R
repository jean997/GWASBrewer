#'@export
mixnorm_to_scale_fam <- function(sigma, pi){
  K <- length(sigma)
  sigma <- check_scalar_or_numeric(sigma, "sigma", K)
  pi <- check_scalar_or_numeric(pi, "pi", K)
  check_01(pi, "pi")

  if(!all(sigma >=0)) stop("sigma must all be >= 0.")

  tot <- sum(pi)
  if(abs(tot -1) > 1e-8) stop("pi should sum to 1.\n")

  if(K == 1){
    f <- function(n, sd, ...){rnorm(n = n, mean = 0, sd = sd)}
    return(f)
  }

  Vbase <- sum(pi*sigma^2)
  f <- function(n, sd, ...){
    Vtarget <- sd^2
    a <- sqrt(Vtarget/Vbase)
    rnormalmix(n = n, sd = a*sigma, pi = pi, mu = 0, return.Z = FALSE)
  }
  return(f)
}


#'@title Simulate from a normal mixture distribution
#'@param n Number of points to simulate
#'@param sigma Standard deviations
#'@param mu Means
#'@param pi Mixture proportions
#'@param return.Z if TRUE, also return a vector of indicators indicating which of the K classes each sample belongs to
#'@return If return.Z=TRUE, returns a list with elements beta (samples) and Z (indicators). Otherwise returns a length n vector of samples.
#'@export
rnormalmix <- function(n, sd, pi, mu =0, return.Z=FALSE){
  K <- length(sd)
  sd <- check_scalar_or_numeric(sd, "sd", K)
  pi <- check_scalar_or_numeric(pi, "pi", K)
  mu <- check_scalar_or_numeric(mu, "mu", K )
  tot <- sum(pi)
  if(abs(tot -1) > 1e-8) stop("pi should sum to 1.\n")
  pi <- pi/tot
  if(!all(sd >=0)) stop("sd must all be >= 0.")


  Z <- sample(1:K, size=n, replace = TRUE, prob=pi)
  beta <- rep(NA, n)
  for(k in 1:K){
    nk <- sum(Z==k)
    if(nk==0) next
    beta[Z==k] <- rnorm(n=nk , mean=mu[k], sd = sd[k])
  }
  if(return.Z) return(list("beta"=beta, "Z"=Z))
  return(beta)
}
