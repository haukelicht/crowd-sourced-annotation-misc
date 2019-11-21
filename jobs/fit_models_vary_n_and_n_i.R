# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title:  Fit BBA models to simulated people-centrism, anti-elitism, and 
#           exclusionism codings datasets that vary in n (No. items) and
#           n[i] (No. judgments per item)
# Author: Hauke Licht
# Date:   2019-04-16
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #


# setup ----

# set the file path
file_path <- file.path("~", "switchdrive", "Documents", "work", "phd", "methods", "crowd-sourced_annotation")

# load required packages

library(dplyr)
library(purrr)
library(tidyr)
library(rjags)
library(ggplot2)

# MCMCM misc 

# set seed
this_seed <- 1234
set.seed(this_seed)

# define JAGS model file path
model_file_path <- file.path(file_path, "models", "beta-binomial_by_annotator_truncated.jags")

# load DIC module
load.module("dic")

# define parameters to estimate
model_parameters <- c(
    "deviance" # Decision Information Criterion (DIC)
    , "pi" # Prevalence of positive class
    , "c" # item classes
    , "theta0", "theta1" # coder specificities and sensitivities, respectively
    , "alpha0", "beta0" # shape parameters of coders' specificity hyperdistribution
    , "alpha1", "beta1" # shape parameters of coders' sensitivity hyperdistribution
  )

# fix number of chains
n_chains <- 3L
# read grid with burnin parameters
n_burnin <- read.csv(file.path(file_path, "config", "burnin_grid.csv"), row.names = 1)
# read grid with thinning parameters
thin <- read.csv(file.path(file_path, "config", "thinning_grid.csv"), row.names = 1)


# load internal helpers

helpers <- c(
  "simulate_binary_codings.R"
  , "get_mcmc_estimates.R"
  , "get_codings_data.R"
  , "transform_betabin_fit_posthoc.R"
)

{sapply(file.path(file_path, "code", "R", helpers), source); NULL}


# determine sample sizes of individual codings data frames ----

# generate group indeces
group_idx <- 0:500 %/% 50
group_idx <- ifelse(group_idx > 3, group_idx-3, 0)[-length(group_idx)]

# inspect group sizes
map_int(unique(group_idx), function(g) sum(group_idx <= g))


# pop_fitted <- readRDS(file.path(file_path, "fits", "betabinom_by_annotator_populism.RData"))
pop_fitted <- readRDS(file.path(file_path, "fits", "betabinom_by_annotator_allf8.RData"))
str(pop_fitted)

pop_params <- pop_fitted %>% 
  # post-hoc transform parameter assignment if reversed
  map(function(fit) {
    if (mean(unlist(fit[, "pi"])) <= .5) 
      return(fit)
    
    transform_betabin_fit_posthoc(fit)
  }) %>% 
  # obtain tidy data frame with relevant parameter estimates
  map(
    get_mcmc_estimates
    , params = c("pi", "alpha0", "alpha1", "beta0", "beta1")
    , use.regex = FALSE
  ) %>% 
  # rename list elements
  # set_names(c("people", "elite", "exclusio"))
  set_names("elite")


# # fit BBA models to simulated people-centrism classificationd ----
# 
# # compute mean posterior estimates of shape parameters of 
# # coder ability hyperdistributions obtained from BBA model fitted to
# # real codings in Hua et al.'s validation data
# people_param_means <- pop_params$people %>% 
#   group_by(parameter) %>% 
#   summarise_at(3, mean) %>% 
#   spread(parameter, est)
# 
# # inspect
# ggplot(data.frame(x = seq(0,1,.1)), aes(x)) + 
#   stat_function(fun=function(x) dbeta(x, people_param_means$alpha0, people_param_means$beta0), aes(color = "theta[0]")) +
#   stat_function(fun=function(x) dbeta(x, people_param_means$alpha1, people_param_means$beta1), aes(color = "theta[1]")) +
#   labs(
#     title = expression(paste("Prior densities ", theta[.*0], " and ", theta[.*1]))
#     , y = "Density"
#     , x = expression(theta[..])
#     , color = ""
#   ) + 
#   theme_bw() +
#   theme(legend.position = "top") 
# 
# # simulate codings using these ability distribution parameters
# sim_people_codings <- simulate_binary_codings(
#   n.items = 500L
#   , n.coders = 40L
#   , missing.rate = .75
#   , pi = people_param_means$pi
#   , alpha0 = people_param_means$alpha0
#   , beta0 = people_param_means$beta0
#   , alpha1 = people_param_means$alpha1
#   , beta1 = people_param_means$beta1
# )
# 
# # inspect
# sim_people_codings$codings %>% 
#   group_by(item) %>% 
#   summarise(n_judgments = n_distinct(coder)) %>% 
#   group_by(n_judgments) %>% 
#   summarise(n = n_distinct(item))
# 
# sim_people_codings$codings %>% 
#   group_by(item) %>% 
#   mutate(position = row_number()) %>% 
#   head(20)
# 
# # generate list containing data frames people-centrism codings that 
# # vary by n (No. items) and n[i] (No. judments per item)
# if (
#   !exists("people_dat")
#   ||
#   !(is.list(people_dat) & all(map_lgl(people_dat, function(x) all(lengths(x) == 2))))
# ) {
#   # create a list of lists with subsetted dataframes
#   people_dat <- map(seq(200, 500, 50), function(g) {
#     
#     dat <- sim_people_codings$codings %>% 
#       filter(item <= g) %>% 
#       group_by(item) %>% 
#       mutate(
#         sample_size = g
#         , position = row_number()
#       ) %>% 
#       ungroup()
#     
#     map(3:10, function(n, d = dat) { 
#       out <- d %>% 
#         filter(position <= n) %>% 
#         select(-position) %>% 
#         mutate(n_i = n)
#       
#       list(data = out)
#     })
#     
#   })
# }
# 
# # inspect 
# str(people_dat, 2)
# 
# # fit BBA models to each codings data frame
# if (
#   # any only one element
#   any(map_lgl(people_dat, function(x) any(lengths(x) == 1)))
#   &&
#   # any has no fitted object
#   any(
#     lengths(
#       idxs <- people_dat %>%
#         map(map_lgl, map_lgl, function(x) any(inherits(x, "mcmc.list"))) %>%
#         map(which)
#     ) == 0
#   )
# ){
# 
#   # for each sample size ...
#   for(g in seq_along(people_dat)) {
#     # ... and number of judgments per item
#     for (n in 1:8) {
#       
#       # skip fitting if mcmc.list already in output object
#       if (n %in% idxs[[g]])
#         next
#       
#       # otherwise initialize ...
#       model <- jags.model(
#         file = model_file_path
#         , data = get_codings_data(people_dat[[g]][[n]]$data)
#         , n.chains = n_chains
#       )
#       
#       # ... udpate ...
#       update(model, n_burnin[n,g])
#       
#       # ... and fit model
#       people_dat[[g]][[n]]$fit <- coda.samples(
#         model
#         , variable.names = model_parameters
#         , n.iter = 1000L*thin[n,g]
#         , thin = thin[n,g]
#       )
#       
#     } 
#   }
# }


# fit BBA models to simulated anti-elitism classifications ----

# compute mean posterior estimates of shape parameters of
# coder ability hyperdistributions obtained from BBA model fitted to
# real codings in Hua et al.'s validation data
elite_param_means <- pop_params$elite %>%
  group_by(parameter) %>%
  summarise_at(3, mean) %>%
  spread(parameter, est)

# inspect
ggplot(data.frame(x = seq(0,1,.1)), aes(x)) +
  stat_function(fun=function(x) dbeta(x, elite_param_means$alpha0, elite_param_means$beta0), aes(color = "theta[0]")) +
  stat_function(fun=function(x) dbeta(x, elite_param_means$alpha1, elite_param_means$beta1), aes(color = "theta[1]")) +
  labs(
    title = expression(paste("Prior densities ", theta[.*0], " and ", theta[.*1]))
    , y = "Density"
    , x = expression(theta[..])
    , color = ""
  ) +
  theme_bw() +
  theme(legend.position = "top")


# simulate codings using these ability distribution parameters
sim_elite_codings <- simulate_binary_codings(
  n.items = 500L
  , n.coders = 40L
  , missing.rate = .75
  , pi = elite_param_means$pi
  , alpha0 = elite_param_means$alpha0
  , beta0 = elite_param_means$beta0
  , alpha1 = elite_param_means$alpha1
  , beta1 = elite_param_means$beta1
)

# inspect
sim_elite_codings$codings %>%
  group_by(item) %>%
  summarise(n_judgments = n_distinct(coder)) %>%
  group_by(n_judgments) %>%
  summarise(n = n_distinct(item))


# generate list containing data frames anti-elitism codings that
# vary by n (No. items) and n[i] (No. judments per item)
if (
  TRUE
  # !exists("elite_dat")
  # ||
  # !(is.list(elite_dat) & all(map_lgl(elite_dat, function(x) all(lengths(x) == 2))))
) {
  # create a list of lists with subsetted dataframes
  elite_dat <- map(seq(200, 500, 50), function(g) {

    dat <- sim_elite_codings$codings %>%
      filter(item <= g) %>%
      group_by(item) %>%
      mutate(
        sample_size = g
        , position = row_number()
      ) %>%
      ungroup()

    map(3:10, function(n, d = dat) {
      out <- d %>%
        filter(position <= n) %>%
        select(-position) %>%
        mutate(n_i = n)

      list(data = out)
    })

  })
}

# inspect
str(elite_dat, 2)


# fit BBA models to each codings data frame
if (
  TRUE
  # # any only one element
  # any(map_lgl(elite_dat, function(x) any(lengths(x) == 1)))
  # &&
  # # any has no fitted object
  # any(
  #   lengths(
  #     idxs <- elite_dat %>%
  #     map(map_lgl, map_lgl, function(x) any(inherits(x, "mcmc.list"))) %>%
  #     map(which)
  #   ) == 0
  # )
) {
  for(g in seq_along(elite_dat)) {
    for (n in 1:8) {

      # skip fitting if mcmc.list already in output object
      if (n %in% idxs[[g]])
        next

      # otherwise initialize ...
      model <- jags.model(
        file = model_file_path
        , data = get_codings_data(elite_dat[[g]][[n]]$data)
        , n.chains = n_chains
      )

      # ... update ...
      update(model, n_burnin[n,g])

      # ..- and fit model
      elite_dat[[g]][[n]]$fit <- coda.samples(
        model
        , variable.names = model_parameters
        , n.iter = 1000L*thin[n,g]
        , thin = thin[n,g]
      )

    }
  }
}




list(
  # people = list(
  #   fitted = people_dat
  #   , sim_params = append(
  #     as.list(people_param_means)
  #     , sim_people_codings$abilities
  #   )
  # )
  # , 
  elite = list(
    fitted = elite_dat
    , sim_params = append(
      as.list(elite_param_means)
      , sim_elite_codings$abilities
    )
  )
  # , exclusion = list(
  #   fitted = exclusio_dat
  #   , sim_params = append(
  #     as.list(exclusio_param_means)
  #     , sim_exclusio_codings$abilities
  #   )
  # )
  , seed = this_seed
) %>% 
  saveRDS(file.path(file_path, "fits", "sim_antielite_allf8_vary_n_and_n_i.RData"))

# # fit BBA models to simulated people-centrism classificationd ----
# 
# # compute mean posterior estimates of shape parameters of 
# # coder ability hyperdistributions obtained from BBA model fitted to
# # real codings in Hua et al.'s validation data
# exclusio_param_means <- pop_params$exclusio %>% 
#   group_by(parameter) %>% 
#   summarise_at(3, mean) %>% 
#   spread(parameter, est)
# 
# # inspect
# ggplot(data.frame(x = seq(0,1,.1)), aes(x)) + 
#   stat_function(fun=function(x) dbeta(x, exclusio_param_means$alpha0, exclusio_param_means$beta0), aes(color = "theta[0]")) +
#   stat_function(fun=function(x) dbeta(x, exclusio_param_means$alpha1, exclusio_param_means$beta1), aes(color = "theta[1]")) +
#   labs(
#     title = expression(paste("Prior densities ", theta[.*0], " and ", theta[.*1]))
#     , y = "Density"
#     , x = expression(theta[..])
#     , color = ""
#   ) + 
#   theme_bw() +
#   theme(legend.position = "top") 
# 
# 
# # simulate codings using these ability distribution parameters
# sim_exclusio_codings <- simulate_binary_codings(
#   n.items = 500L
#   , n.coders = 40L
#   , missing.rate = .75
#   , pi = exclusio_param_means$pi
#   , alpha0 = exclusio_param_means$alpha0
#   , beta0 = exclusio_param_means$beta0
#   , alpha1 = exclusio_param_means$alpha1
#   , beta1 = exclusio_param_means$beta1
# )
# 
# # inspect
# sim_exclusio_codings$codings %>% 
#   group_by(item) %>% 
#   summarise(n_judgments = n_distinct(coder)) %>% 
#   group_by(n_judgments) %>% 
#   summarise(n = n_distinct(item))
# 
# # generate list containing data frames people-centrism codings that 
# # vary by n (No. items) and n[i] (No. judments per item)
# if (TRUE) {
#   exclusio_dat <- map(seq(200, 500, 50), function(g) {
#     
#     dat <- sim_exclusio_codings$codings %>% 
#       filter(item <= g) %>% 
#       group_by(item) %>% 
#       mutate(
#         sample_size = g
#         , position = row_number()
#       ) %>% 
#       ungroup()
#     
#     map(3:10, function(n, d = dat) { 
#       out <- d %>% 
#         filter(position <= n) %>% 
#         select(-position) %>% 
#         mutate(n_i = n)
#       
#       list(data = out)
#     })
#     
#   })
# }
# 
# # fit BBA models to these codings data frames
# if (
#   # any only one element
#   any(map_lgl(exclusio_dat, function(x) any(lengths(x) == 1)))
#   && 
#   # any has no fitted object
#   any(
#     lengths(
#       idxs <- exclusio_dat %>% 
#         map(map_lgl, map_lgl, function(x) any(inherits(x, "mcmc.list"))) %>% 
#         map(which)
#     ) == 0
#   )
# ){
#   for(g in seq_along(exclusio_dat)) {
#     for (n in 1:8) {
#       
#       # skip fitting if mcmc.list already in output object
#       if (n %in% idxs[[g]])
#         next
#       
#       # otherwise initialize
#       model <- jags.model(
#         file = model_file_path
#         , data = get_codings_data(exclusio_dat[[g]][[n]]$data)
#         , n.chains = n_chains
#       )
#       
#       # ... update ...
#       update(model, n_burnin[n,g])
#       
#       # ... and fit model
#       exclusio_dat[[g]][[n]]$fit <- coda.samples(
#         model
#         , variable.names = model_parameters
#         , n.iter = 1000L*thin[n,g]
#         , thin = thin[n,g]
#       )
#       
#     } 
#   }
# }


# save simulated data and fitted objects ----

# list(
#   people = list(
#     fitted = people_dat
#     , sim_params = append(
#       as.list(people_param_means)
#       , sim_people_codings$abilities
#     )
#   )
#   , elite = list(
#     fitted = elite_dat
#     , sim_params = append(
#       as.list(elite_param_means)
#       , sim_elite_codings$abilities
#     )
#   )
#   , exclusion = list(
#     fitted = exclusio_dat
#     , sim_params = append(
#       as.list(exclusio_param_means)
#       , sim_exclusio_codings$abilities
#     )
#   )
#   , seed = this_seed
# ) %>% 
#   saveRDS(file.path(file_path, "fits", "sim_populism_codings_vary_n_and_n_i.RData"))





