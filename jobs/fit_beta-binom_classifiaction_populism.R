# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title:  Fit beta-binomial by coder model to crowd-sourced validation 
#          dataset on populist instances among political statements on 
#          Twitter and Facebook.
# Author: Hauke Licht
# File:   fit_beta-binom_classification_populism.R
# Date:   2019-03-06
# Update: 2019-03-14
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #


# set the file path
file_path <- file.path("~", "switchdrive", "Documents", "work", "phd", "methods", "crowd-sourced_annotation")

# load required namespaces

library(dplyr)
library(rjags)

# set seed
set.seed(1234)

# load internal helpers

helpers <- c("get_mcmc_estimates.R", "get_codings_data.R", "transform_betabin_fit_posthoc")
sapply(file.path(file_path, "code", "R", helpers), source)


# load and prep data ----
dat <- data.table::fread(file.path(file_path, "exdata", "f1248095.csv"), stringsAsFactors = FALSE)

# inspect
table(dat$filter)

codings <- dat %>% 
  # keep only judgements that are 'ok'
  filter("ok" == filter) %>%
  # replace "" with NA in string vectors
  mutate_if(is.character, Vectorize(function(x) if (x == "") NA_character_ else x)) %>%
  mutate(
    # judgement index
    index = row_number(),
    # item index
    item = group_indices(., `_unit_id`) ,
    # coder index
    coder = group_indices(., `_worker_id`),
    # populism indiactors
    elites = elites == "yes",
    exclusionary = exclusionary == "yes",
    people = people == "yes",
    populist = people & elites,
    right_populist = people & elites & exclusionary
  ) %>% 
  tbl_df()

(n_judgments <- nrow(codings))
(n_coders <- length(unique(codings$coder)))
(n_items <- length(unique(codings$item)))

# fit Bayesian beta-binomial by annotator models ----

# load DIC modeul
load.module("dic")

# global model parameters
n_chains <- 3
model_file_path <- file.path(file_path, "models", "beta-binomial_by_annotator.jags")
fit_file_path <- file.path(file_path, "fits", "betabinom_by_annotator_populism.RData")

init_vals <- lapply(1:n_chains, function(chain) {
  
  out <- list()
  out[["pi"]] <- .2 + rnorm(1, 0, .05)
  out[[".RNG.name"]] <- "base::Wichmann-Hill"
  out[[".RNG.seed"]] <- 1234
  
  return(out)
})

model_parameters <- c(
  "deviance"
  , "pi"
  , "c"
  , "theta0", "theta1"
  , "alpha0", "beta0"
  , "alpha1", "beta1"
)

# inspect
system(paste("less", model_file_path))


# fit to people-centrism ---

# subset codings
people_codings <- codings %>%
  mutate(judgment = as.integer(people)) %>% 
  filter(!is.na(judgment)) %>% 
  select(index, item, coder, judgment) 

# contruct model-compatible MCMC data object
people_mcmc_data <- get_codings_data(people_codings)

# initialize model
people_mcmc_model <- jags.model(
  file = model_file_path
  , data = people_mcmc_data
  , inits = init_vals
  , n.chains = n_chains
)

update(people_mcmc_model, 10000)

people_mcmc_fit <- coda.samples(
  people_mcmc_model
  , variable.names = model_parameters
  , n.iter = 40000
  , thin = 20
)


# fit to anti-elitism ----

# subset codings
elite_codings <- codings %>%
  mutate(judgment = as.integer(elites)) %>% 
  filter(!is.na(judgment)) %>% 
  select(index, item, coder, judgment) 

# contruct model-compatible MCMC data object
elite_mcmc_data <- get_codings_data(elite_codings)

# initialize model
elite_mcmc_model <- jags.model(
  file = model_file_path
  , data = elite_mcmc_data
  , inits = init_vals
  , n.chains = n_chains
)

update(elite_mcmc_model, 5000)

elite_mcmc_fit <- coda.samples(
  elite_mcmc_model
  , variable.names = model_parameters
  , n.iter = 15000
  , thin = 15
)



# fit to exclusionism ---

# subset codings
exclusio_codings <- codings %>%
  mutate(judgment = as.integer(exclusionary)) %>% 
  filter(!is.na(judgment)) %>% 
  select(index, item, coder, judgment) 

# contruct model-compatible MCMC data object
exclusio_mcmc_data <- get_codings_data(exclusio_codings)

exclusio_mcmc_model <- jags.model(
  file = model_file_path
  , data = exclusio_mcmc_data
  , inits = init_vals
  , n.chains = n_chains
)

update(exclusio_mcmc_model, 10000)

exclusio_mcmc_fit <- coda.samples(
  exclusio_mcmc_model
  , variable.names = model_parameters
  , n.iter = 100000
  , thin = 50
)

# save fitted objects ----

list(
  people_mcmc_fit = people_mcmc_fit
  , elite_mcmc_fit  = elite_mcmc_fit
  , exclusio_mcmc_fit = exclusio_mcmc_fit
) %>% 
  saveRDS(fit_file_path)
