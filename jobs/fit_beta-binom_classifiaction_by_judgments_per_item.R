# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title:  Fit Bayesian beta-binomial by annotator models to subsamples of 
#          simulated judgements where subsamples vary in the number of 
#          judgements per item
# Author: Hauke Licht
# File:   fit_beta-binom_classification_by_judgments_per_item.R
# Date:   2019-03-05
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #


# setup ----

# load required namespaces
library(rjags)
library(purrr)
library(dplyr)
library(tidyr)

# set seed
set.seed(1234)

# load DIC modeul
load.module("dic")

# define helper function

#' Function creating JAGS model data input from codings data
#' 
#' @description Given a codings dataframe with columngs 'item', 'coder', and 'judgement', 
#'     function returns a list object that can be passed to \code{data} argument 
#'     of \[rjags]{jags.model} function fitting a beta-binomial by annotator model
#'     
#' @param Given codings dataframe with columngs 'item', 'coder', and 'judgement'
#' 
#' @return A list object with elements 
#'     'N' (the number of items, scalar integer),
#'     'M' (the number of coders, scalar integer),
#'     'J' (the number of judgements, scalar integer), and
#'     'y' (dataframe with columns 'item', 'coder' and 'judgment').
get_codings_data <- function(codings.df) {
  
  out <- list() 
  
  out$N <- length(unique(codings.df$item))
  out$M <- length(unique(codings.df$coder))
  out$J <- length(codings.df$judgment)
  out$y <- codings.df %>% 
    arrange(item, coder) %>% 
    select(item, coder, judgment)
  
  return(out)
}

# simulate codings data ----

n_items <- 1000
n_coders <- 20
pi <- .2
alpha0 <- 40
beta0 <- 8
alpha1 <- 20
beta1 <- 8
missing_rate <- .5

source("./jobs/sim_codings_data.R")

sim_codings


# General model setup ----

n_chains <- 3
n_burnin <- 500
n_iter <- 1000
n_thinning <- 1
model_file_path <- file.path(".", "models", "beta-binomial_by_annotator.jags")
fit_file_path <- file.path(".", "fits", "betabinom_by_annotator_vary_n_i.RData")

# inspect
system(paste("less", model_file_path))

# creat list of lists with initialization values for thetas and pi 
init_vals <- lapply(1:n_chains, function(chain) {
  
  out <- list()
  out[["theta0"]] <- rbeta(n_coders, 8, 2)
  out[["theta1"]] <- rbeta(n_coders, 8, 2)
  out[["pi"]] <- rbeta(1, 2, 8)
  out[[".RNG.name"]] <- "base::Wichmann-Hill"
  out[[".RNG.seed"]] <- 1234
  
  return(out)
})

# inspect
map(init_vals, "theta0")
map(init_vals, "theta1")
map(init_vals, "pi")


# do only run once
if (!file.exists(fit_file_path)) {
  s_fitted <- map(
    3:ceiling(n_coders*missing_rate)
    , function(
        s
        , max.s = ceiling(n_coders*missing_rate)
        , n.trials = 3L
        , model.file = model_file_path
        , init.vals = init_vals
        , n.chains = n_chains
        , n.burnin = n_burnin
        , n.iter = n_iter
        , n.thinning = n_thinning
    ){
    
    out <- list()
    
    for (t in 1:n.trials) {
      # sample s judgements per item
      dat <- if (s < max.s) {
        sim_codings %>% 
          group_by(item) %>% 
          sample_n(s, replace = FALSE)
      } else {
        sim_codings
      }
      
      # write sampled indices to output list
      out$samples <- dat %>% select(item, coder)
      
      # try to initialize model on sampled judgments
      out$model <- tryCatch(
        jags.model(
          file = model.file
          , data = get_codings_data(dat)
          , inits = init.vals
          , n.chains = n.chains
        )
        , error = function(err) err)
        
      # proceed if initialization succeeded
      if (inherits(out$model, "jags"))
          break
      
      # else return without fit 
      if (t == n.trials){
        out$model <- NA
        out$fit <- NA
        warning("Could not initialize model on subsample of ", s, " judgments after ", t, " trials.")
        return(out)
      }
    }
  
    # try to update model
    updated <- tryCatch(update(out$model, n.burnin), error = function(err) err)
      
    # if update fails, return without fit 
    if (inherits(updated, "error")){
      out$fit <- NA
      warning("Could not update model on subsample of ", s, " judgments.")
      return(out)
    }
    
    # try to fit model
    out$fit <- tryCatch(
      coda.samples(
        out$model
        , variable.names = c(
          "deviance"
          , "pi"
          , "c"
          , "theta0", "theta1"
          # , "alpha0", "beta0"
          # , "alpha1", "beta1"
        )
        , n.iter = n.iter
        , thin = n.thinning
      )
      , error = function(err) err
    )
    
    # if fitting fails, return without fit 
    if (inherits(out$fit, "error")){
      out$fit <- NA
      warning("Could not fit model on subsample of ", s, " judgments.")
      return(out)
    }
    
    # else, return without output list including fit 
    return(out)
  })
} else {
  s_fitted <- readRDS(fit_file_path)
}

