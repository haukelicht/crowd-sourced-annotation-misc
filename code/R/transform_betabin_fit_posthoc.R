#' Transform MCMC estimates of beta-binomial by annotator model to opposite assignment
#' 
#' @param mcmc.fit an \code[rjags]{mcmc} \code[rjags]{mcmc.list} object 
#'      containing parameter estimates of beta-binomial by annotator model
#'      
#' @return \code{mcmc.fit} with transformed parameter estimates
transform_betabin_fit_posthoc <- function(mcmc.fit){
  
  trans_params <- function(fit) {
    params <- colnames(fit)
    
    out <- matrix(nrow = nrow(fit), ncol = ncol(fit))
    colnames(out) <- rep(NA_character_, ncol(fit))
    
    if (length(dic <- grep("^deviance$", params)) > 0){
      out[, dic] <- fit[, dic]
      colnames(out)[dic] <- params[dic]
    }
    
    if (length(pi <- grep("^pi$", params)) > 0){
      out[, pi] <- 1-fit[, pi]
      colnames(out)[pi] <- params[pi]
    }
    
    if (length(cs <- grep("^c\\[\\d+\\]$", params)) > 0){
      out[, cs] <- 1-fit[, cs]
      colnames(out)[cs] <- params[cs]
    }
    
    if (length(theta0s <- grep("^theta0\\[\\d+\\]$", params)) > 0 & length(theta1s <- grep("^theta1\\[\\d+\\]$", params)) > 0){
      out[, theta1s] <- 1-fit[, theta0s]
      colnames(out)[theta1s] <- params[theta1s]
      
      out[, theta0s] <- 1-fit[, theta1s]
      colnames(out)[theta0s] <- params[theta0s]
    }
    
    if (length(par0 <- grep("^(alpha|beta)0$", params)) > 0 & length(par1 <- grep("^(alpha|beta)1$", params)) > 0){
      out[, rev(par1)] <- fit[, par0]
      colnames(out)[par1] <- params[par1]
      
      out[, rev(par0)] <- fit[, par1]
      colnames(out)[par0] <- params[par0]
    }
    
    return(out)
  }
  
  if (inherits(mcmc.fit, "mcmc.list")){
    out <- map(mcmc.fit, trans_params)
    class(out) <- class(mcmc.fit)
  } else if (inherits(mcmc.fit, "mcmc")){
    out <- trans_params(mcmc.fit)
    class(out) <- class(mcmc.fit)
  } else{
    stop(sprintf("Cannot apply `transform_betabin_fit_posthoc` to object of class '%a'", class(mcmc.fit)))
  }
  
  return(out)
}
