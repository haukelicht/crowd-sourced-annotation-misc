

# set the file path
file_path <- file.path("~", "switchdrive", "Documents", "work", "phd", "methods", "crowd-sourced_annotation")

# load required namespaces

library(dplyr)
library(purrr)
source(file.path(file_path, "code", "R", "get_mcmc_estimates.R"))


sim_fitted <- readRDS(file.path(file_path, "fits", "sim_antielite_allf8_vary_n_and_n_i.RData"))

cs <- map_df(
  sim_fitted$elite$fitted
  , function(g) {
    map_df(g, function(l) {
      get_mcmc_estimates(fit.obj = l$fit, params = "c\\[\\d+\\]") %>%
        mutate(
          n_i = unique(l$data$n_i)
          , sample_size = unique(l$data$sample_size)
        )
    })
  }
)

list(
  codings = sim_fitted$elite$fitted[[7]][[8]]$data
  , pi = map_df(sim_fitted$elite$fitted, function(g) {
    map_df(g, function(l) {
      get_mcmc_estimates(fit.obj = l$fit, params = "pi", use.regex = FALSE) %>%
        mutate(
          n_i = unique(l$data$n_i)
          , sample_size = unique(l$data$sample_size)
        )
    })
  })
  , c = cs
  , theta = map_df(sim_fitted$elite$fitted, function(g) {
    map_df(g, function(l) {
      get_mcmc_estimates(fit.obj = l$fit, params = "theta\\d\\[\\d\\]") %>%
        mutate(
          n_i = unique(l$data$n_i)
          , sample_size = unique(l$data$sample_size)
        )
    })
  })
  , hyperpars = map_df(sim_fitted$elite$fitted, function(g) {
    map_df(g, function(l) {
      map_df("(alpha|beta)[01]", get_mcmc_estimates, fit.obj = l$fit) %>%
        mutate(
          n_i = unique(l$data$n_i)
          , sample_size = unique(l$data$sample_size)
        )
    })
  })
  , postclass_quality = sim_fitted$elite$fitted[[7]][[8]]$data %>%
      select(item, true_class) %>%
      unique() %>%
      mutate(parameter = sprintf("c[%s]", item)) %>%
      right_join(cs) %>%
      group_by(sample_size, n_i, chain, iter) %>%
      summarise(
        n_judgments = n_distinct(item)
        , Accuracy = sum(est == true_class)/n_distinct(item)
        , tp = sum(est == 1 & true_class == 1)
        , fn = sum(est == 0 & true_class == 1)
        , fp = sum(est == 1 & true_class == 0)
        , tn = sum(est == 0 & true_class == 0)
        , TPR = tp / (tp + fn)
        , TNR = tn / (tn + fp)
        , FPR = fp / (tn + fp)
        , FNR = fn / (tp + fn)
        , Precision = tp / (tp + fp)
        , Recall = TPR
        , `F1-score` = 2*((Precision*Recall)/(Precision + Recall))
      ) %>%
      ungroup()
  , created_ts <- Sys.time()
) %>% 
  saveRDS(file.path(file_path, "data", "post_est_antielite_allf8_vary_n_and_n_i.RData"))

