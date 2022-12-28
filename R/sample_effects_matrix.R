sample_effects_matrix <- function(J, pi, sigma, f,
                                  sporadic_pleiotropy,
                                  pi_exact,
                                  h2_exact,
                                  R_LD = NULL, af = NULL){

  M <- length(pi)
  #skip checks for internal function
  #pi <- check_scalar_or_numeric(pi, "pi", M)
  #sigma <- check_scalar_or_numeric(sigma, "sigma", M)
  #pi <- check_01(pi)

  if(all(sigma == 0)){
    eff <- matrix(0, nrow = J, ncol = M)
    return(eff)
  }

  if(pi_exact){
    sig_j <- sqrt( (1/round(pi*J)) * sigma^2)
  }else{
    sig_j <- sqrt( (1/(pi*J)) * sigma^2)
  }

  if(pi_exact & sporadic_pleiotropy){
    eff <- purrr::map(seq(M), function(i){
      n <- round(pi[i]*J)
      val <- rep(0, J)
      if(n > 0){
        t <- sample(seq(J), size = n, replace = FALSE)
        val[t] <- f(n=n, sd = sig_j[i], af = af)
      }
      return(val)
    }) %>% do.call(cbind, .)
  }else if(sporadic_pleiotropy){
    eff <- purrr::map(seq(M), function(i){
      t <- rbinom(n=J, size=1, prob = pi[i])
      n <- sum(t==1)
      if(n > 0) t[t==1] <- f(n=n, sd = sig_j[i], af = af)
      return(t)
    }) %>% do.call(cbind, .)
  }else{
    ix <- which(sigma > 0)
    if(sum(pi[ix]) > 1){
      stop("You have requested too many traits and too many causal variants to use sporadic_pleiotropy = FALSE.\n")
    }
    if(pi_exact){
      nz_ix <- rep(0, J)
      for(i in seq(M)){
        n <- round(pi[i]*J)
        ii <- sample(which(nz_ix == 0), size = n, replace = FALSE)
        if(n > 0) nz_ix[ii] <- i
      }
    }else{
      p <- sum(pi_theta[ix])
      nz_ix <- sample(c(0, ix), size = J, replace = TRUE, prob = c(1-p, pi[ix]) )
    }
    eff <- purrr::map(seq(M), function(i){
      t <- which(nz_ix == i)
      n <- length(t)
      val <- rep(0, J)
      if(n > 0) val[t] <- f(n=n, sd = sig_j[i], af = af)
      return(val)
    }) %>% do.call(cbind, .)
  }

  if(h2_exact){
    h2_eff <- compute_h2(b_joint_std = eff, R_LD = R_LD, af = af)
    scale <- sqrt((sigma^2)/h2_eff)
    eff <- t(t(eff)*scale)
  }
  return(eff)
}



