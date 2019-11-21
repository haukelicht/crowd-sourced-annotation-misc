# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title:  Estimate intercoder reliability for different values of n[i], the number 
#         of judgements aggregated per item
# Author: Hauke Licht
# File:   est_intercoder_reliability_by_n_i.R
# Date:   2019-03-06
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #

# setup ----

# load required packages 
library(dplyr)
library(icr)
library(irr)
library(purrr)
library(tidyr)

# define file path
file_path <- file.path("~", "switchdrive", "Documents", "work", "phd", "methods", "crowd-sourced_annotation")

# laod simulated codings data
sim_codings <- readRDS(file.path(file_path, "data", "sim_codings_n1000_n_i10.RData")) 
# NOTE: (see <file_path>/jobs/sim_codings_data.R)

# load varying-n[i] annotation model fit data  
s_fitted <- readRDS(file.path(file_path, "fits", "betabinom_by_annotator_vary_n_i.RData"))
# NOTE: n[i] is the number of judgments per item aggregated to obtain posterior classifications
#       eight different beta-binomial by annotator models were fit for values of n[i] = 3,...,10
#       to n[i] randomly sampled judgements per item 
#       (see <file_path>/jobs/fit_beta-binom_classification_by_judgments_per_item.R)

# global parameters
boot_size <- 1000L


# Compute Fleiss's kappa for samples n[i] = 3,...,10 ----

#' Helper function passed to \code[boot]{boot}'s \code{statistic} argument to compute 
#'     Fleiss's kappa
#'     
#' @description Given a Item X Coders judgements/ratings dataframe and row indices,
#'     function returns Fleiss's kappa
#'     
#' @param d a Item X Coders judgements/ratings dataframe
#' 
#' @param idxs row indices used to sample from \code{d} when bootstrapping
#' 
#' @return Fleiss's kappa in the (subsetted) ratings daat as a scalar numeric value 
#' 
#' @importFrom irr kappam.fleiss
compute_fleiss_kappa <- function(d, idxs){
  irr::kappam.fleiss(
    ratings = d[idxs, ]
    , exact = FALSE
    , detail = FALSE
  )$value
}

# compute 
kappas <- map_df(
  .x = s_fitted
  , .id = "s"
  , .f = function(s, boot.size = boot_size) {
    out <- s$samples %>% 
      group_by(item) %>% 
      mutate(row = row_number()) %>% 
      left_join(sim_codings, by = c("item", "coder")) %>% 
      select(-true_class, -coder) %>% 
      # `compute_fleiss_kappa` expects Item X Coder format, hence reshape
      spread(row, judgment, sep = "_") %>% 
      ungroup() %>% 
      select(-item) %>% 
      boot::boot(
        data = .
        , statistic = compute_fleiss_kappa
        , R = boot.size
        , sim = "ordinary"
        , stype = "i"
        , parallel = "multicore"
      )
    
    tibble(kappa = out$t[,1])
  }) %>% 
  mutate(s = as.integer(s) + 2L)


# Compute Krippendorff's alpha for samples n[i] = 3,...,10 ----

kalphas <- map_df(
  .x = s_fitted
  , .id = "s"
  , .f = function(s, boot.size = boot_size) {
    s$samples %>% 
      ungroup() %>% 
      left_join(sim_codings, by = c("item", "coder")) %>% 
      select(-true_class) %>% 
      # icr::krippalpha expects Coder X Item format, hence reshape
      spread(item, judgment, sep = "_") %>% 
      select(-coder) %>% 
      as.matrix() %>% 
      icr::krippalpha(
        data = .
        , metric = "nominal"
        , bootstrap = TRUE
        , nboot = boot.size
        , cores = 2
      ) %>% 
      .$bootstraps %>% 
      tibble(kalpha = .)
  }) %>% 
  mutate(s = as.integer(s) +2L)

# kalphas %>% 
#   group_by(s) %>% 
#   summarise(
#     mean = mean(kalpha)
#     , q5 = quantile(kalpha, .05)
#     , q95 = quantile(kalpha, .95)
#   ) %>% 
#   ggplot(aes(x = s, y = mean)) +
#   geom_violin(
#     data = kalphas
#     , mapping = aes(x = s, y = kalpha, group = s)
#     , color = NA
#     , fill = "grey"
#     , alpha= .75
#   ) +
#   geom_point() +
#   geom_linerange(aes(ymin = q5, ymax = q95)) + 
#   theme_bw() + 
#   labs(
#     title = expression(paste(
#       "Intercoder reliability as function of ", n[i],", the number of judgments aggregated per item."
#     ))
#     , subtitle = "Shaded violins summarize 1000 bootstrapped values, dots indicate mean values, vertical lines 90% intervalls."
#     , x = expression(n[i])
#     , y = expression(paste("Krippendorff's ", alpha))
#   ) + 
#   scale_x_continuous(breaks = 1:8, labels = 3:10)


# Save bootstrapped ICR reliability estimates ----

list(
  codings_data = sim_codings
  , sampled_judgements = map_df(
    .x = s_fitted
    , .id = "s"
    , .f = function(s) s$samples %>% ungroup()
  ) %>% 
    mutate(s = as.integer(s) + 2L)
  , boot_size = boot_size
  , krippendorffs_alphas = kalphas
  , fleiss_kappas = kappas
) %>% 
  saveRDS(., file.path(file_path, "data", "sim_codings_bs_realiability_n_i3-10.RData"))

