

# set the file path
file_path <- file.path("~", "switchdrive", "Documents", "work", "phd", "methods", "crowd-sourced_annotation")

# load required namespaces

library(dplyr)
library(purrr)
source(file.path(file_path, "code", "R", "get_mcmc_estimates.R"))


test_fitted <- readRDS(file.path(file_path, "fits", "bba_antielite_testf8_vary_n_i.RData"))

# DIC estimates
dics <- map_df(
  test_fitted
  , function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "deviance", use.regex = FALSE) %>%
      mutate(n_i = unique(l$data$n_i)) 
  }
)

# prevalence estimates
pis <- map_df(
  test_fitted
  , function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "pi", use.regex = FALSE) %>%
      mutate(n_i = unique(l$data$n_i)) 
  }
)

# hyperparameter estimates
hyperpars <- map_df(
  test_fitted
  , function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "(alpha|beta)[01]", use.regex = TRUE) %>%
      mutate(n_i = unique(l$data$n_i)) 
  }
)

# instance class estimates
cs <- map_df(
  test_fitted
  , function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "c\\[\\d+\\]") %>%
      mutate(
        item = as.integer(str_extract(parameter, "\\d+")) 
      ) %>% 
      left_join(unique(select(l$data, id_str, item, n_i)))
  }
)

# coder ability parameter estimates
thetas <- map_df(
  test_fitted
  , function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "theta[01]\\[\\d+\\]") %>%
      mutate(
        coder = as.integer(str_replace(parameter, "theta[01]\\[(\\d+)\\]", "\\1"))
        , parameter = str_extract(parameter, "theta[01]")
      ) %>% 
      left_join(unique(select(l$data, f8_worker_id, coder, n_i)))
  }
)


voted <- map_df(
  test_fitted
  , function(l) {
      l$data %>%
        group_by(item, n_i) %>%
        summarise(
          n_pos = sum(judgment)
          , n_neg = sum(!judgment)
          , tie_breaker = rbernoulli(1)
          , voted = case_when(
            n_pos > n_neg ~ 1L
            , n_pos < n_neg ~ 0L
            , TRUE ~ as.integer(tie_breaker)
          )
        )
  }
)




list(
  codings = select(test_fitted[[8]]$data, item, id_str, n_i, coder, f8_worker_id, judgment, elites)
  , dic = select(dics, n_i, chain, iter, parameter, est)
  , pi = select(pis, n_i, chain, iter, parameter, est)
  , c = select(cs, n_i, chain, iter, item, id_str, parameter, est)
  , theta = select(thetas, n_i, chain, iter, coder, f8_worker_id, parameter, est)
  , hyperpars = select(hyperpars, n_i, chain, iter, parameter, est)
  , voted = select(voted, n_i, item, n_pos, n_neg, tie_breaker, voted)
  , created_ts <- Sys.time()
) %>% 
  saveRDS(file.path(file_path, "data", "post_est_antielite_testf8_vary_n_i.RData"))

