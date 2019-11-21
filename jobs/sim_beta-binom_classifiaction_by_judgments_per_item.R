# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title:  Fit Bayesian beta-binomial by annotator models to subsamples of 
#          simulated judgements where subsamples vary in the number of 
#          judgements per item
# Author: Hauke Licht
# File:   sim_beta-binom_classification_by_judgments_per_item.R
# Date:   2019-03-05
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #



# setup ----

# load required namespaces
library(rjags)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)
library(irr)
library(icr)

# set seed
set.seed(1234)

# load DIC modeul
load.module("dic")

# define helper function

#' Obtain parameter estimates from MCMC fit object
#' 
#' @descritpion Given an \code{mcmc.list} object created by fitting a JAGS model,
#'      and a list of parameters of interest, the function returns the 
#'      parameter estimates for each chain and all iterations in a tidy dataframe
#' 
#' @param fit.obj A \code{mcmc.list} containing the parameter estimates for all iterations and all chains obtained by fitting a JAGS model using \code[rjags]{}
#' 
#' @param params A character vector specifiying the parameters of interest.
#'    If \code{use.regex = TRUE} (the default), regular expressions will be expected.
#' 
#' @param use.regex Logical, specifying whether \code{params} is defined using regular expressions
#'    Defaulting to \code{TRUE}
#'
#' @return A tidy dataframe with columns
#'    'parameter' (the paramenter, character),
#'    'chain' (chain counter, integer),
#'    'est' (the estimate), and
#'    'iter' (iteration counter, integer)
#'
#' @note Functionality comparable to the \code[ggmcmc]{ggs} function in the \code[ggmcmc]{ggmcmc} pacakge.
#'
#' @importFrom tibble tibble
#' @importFrom purrr map_df
#' @importFrom dplyr mutate
#' @importFrom dplyr row_number
#' @importFrom dplyr `%>%`
#' 
get_mcmc_estimates <- function(fit.obj, params, use.regex = TRUE){
  
  if (!inherits(fit.obj, "mcmc.list"))
    stop("`fit.obj` must be a 'mcmc.list' object")
  
  # test number of chains (equals No. top-level list elements)
  if ( (n_chains <- length(fit.obj)) == 0)
    stop("`fit.obj` has zero length")
  
  # get parameters contained in fit object
  these_params <- colnames(fit.obj[[1]])
  
  # get index positions of parameters to be obtained
  if (use.regex) {
    idxs <- which(grepl(params, these_params))
  } else {
    if (any((miss_params <- !params %in% these_params)))
      warning(
        "The following `params` are not contained in `fit.obj`: ", 
        paste0("'", params[miss_params], "'", collapse = ", ")
      )
    idxs <- which(these_params %in% params)
  }
  
  # obtain parameter estimates in tidy data frame
  suppressWarnings(
    purrr::map_df(these_params[idxs], function(param, obj = fit.obj){
      purrr::map_df(1:n_chains, function(c, o = obj, p = param) {
        tibble::tibble(
          parameter = param
          , chain = c
          , est = o[[c]][, p]
        ) %>% 
          dplyr::mutate(iter = dplyr::row_number())
      })
    })
  )
}


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
model_file_path <- "./models/beta-binomial_by_annotator.jags"

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


fit_file_path <- file.path(".", "fits", "betabinom_by_annotator_vary_n_i.RData")
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
    
    print(s)
    out <- list()
    
    for (t in 1:n.trials) {
      dat <- if (s < max.s) {
        sim_codings %>% 
          group_by(item) %>% 
          sample_n(s, replace = FALSE)
      } else {
        sim_codings
      }
      
      out$samples <- dat %>% select(item, coder)
      
      out$model <- tryCatch(
        jags.model(
          file = model.file
          , data = get_codings_data(dat)
          , inits = init.vals
          , n.chains = n.chains
        )
        , error = function(err) err)
        
      if (inherits(out$model, "jags"))
          break
      
      if (t == n.trials){
        out$model <- NA
        out$fit <- NA
        warning("Could not initialize model on subsample of ", s, " judgments after ", t, " trials.")
        return(out)
      }
    }
  
    updated <- tryCatch(update(out$model, n.burnin), error = function(err) err)
      
    if (inherits(updated, "error")){
      out$fit <- NA
      warning("Could not update model on subsample of ", s, " judgments.")
      return(out)
    }
    
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
    
    if (inherits(out$fit, "error")){
      out$fit <- NA
      warning("Could not fit model on subsample of ", s, " judgments.")
      return(out)
    }
    
    return(out)
  })
} else {
  s_fitted <- readRDS(fit_file_path)
}

str(s_fitted, 1)

map(s_fitted, str, 1)


map(s_fitted, "fit", .id = "s")


plot(s_fitted[[8]]$fit[, "pi"])
gelman.plot(s_fitted[[8]]$fit[, "pi"])


post_pis <- map_df(.x = s_fitted, .id = "s", .f = function(s){
  get_mcmc_estimates(s$fit, "pi")
})

post_pis %>% 
  group_by(s, parameter) %>% 
  summarise(
    mean = mean(est)
    , q5 = quantile(est, .05)
    , q95 = quantile(est, .95)
  ) %>% 
  ggplot(aes(x = s, y = mean)) +
  geom_violin(
    data = post_pis
    , mapping = aes(x = s, y = est)
    , color = NA, fill = "grey", alpha = .75
  ) +
  geom_point() +
  geom_linerange(aes(ymin = q5, ymax = q95)) +
  scale_x_discrete(labels = 3:10) +
  theme_bw() +
  labs(
    title = expression(paste(
      "Location and dispersion of posterior density ", pi, " as function of ", n[i],", the number of judgments aggregated per item."
    ))
    , subtitle = "Shaded violines summarize posterior values, dots indicate mean, vertical lines 90% credibility intervalls."
    , y = expression(pi)
    , x = expression(n[i])
    , caption = expression(paste(
      "Estimates obtained by fitting beta-binomial by annotator models to randomly sampled subsets of ", 
      n[i], 
      " judgements per item."
    ))
  ) + 
  theme(
    axis.title.y = element_text(angle = 0, vjust = .5)
    , plot.caption = element_text(hjust = 0)
  )


post_cs <- map_df(.x = s_fitted, .id = "s", .f = function(s){
  get_mcmc_estimates(s$fit, "c\\[\\d+\\]")
}) %>% 
  mutate(s = as.integer(s))


post_cs_sds <- post_cs %>% 
  left_join(
    sim_codings %>% 
      select(item, true_class) %>% 
      unique() %>% 
      mutate(parameter = paste0("c[", item, "]"))
  ) %>% 
  group_by(s, item) %>% 
  summarise(sd = sd(est))


$\text{SD}_i = \frac{1}{n_t} \sum_t (c_{it} - \bar{c}_i)^2$, where $t$ indexes the $t^{th}$ iterations

$\text{SD}_i$ is maximized when $\bar{c}_i = .5$, since
- $n_t = C \times T = 3 \times 1000 \forall i \in 1,\ldots, n$
  - $c_{it} \in \{0, 1\}$
  this maximum corresponds to an item that is classified as positive exaclty $(C \times T)/2$ times   

$\text{SD}_i$ is minimized when $\bar{c}_i \in \{0,1\}$, i.e., when i is classified unambigously across iterations and chains.

post_cs_sds %>% 
  group_by(s) %>% 
  summarise(
    mean = mean(sd)
    , q5 = quantile(sd, .05)
    , q95 = quantile(sd, .95)
  ) %>% 
  ggplot(aes(x = s, y = mean)) +
  geom_violin(
    data = post_cs_sds
    , mapping = aes(x = s, y = sd, group = s)
    , color = NA, fill = "grey", alpha = .75
  ) +
  scale_y_continuous(limits = c(0,.6)) +
  scale_x_discrete(labels = 3:10) +
  geom_point() + 
  geom_linerange(aes(ymin = q5, ymax = q95)) +
  theme_bw() +
  labs(
    title = expression(paste(
      "Empirical distribution of item-level posterior classification uncertainty as function of ", n[i],", the number of judgments aggregated per item."
    ))
    , subtitle = "Shaded violins summarize posterior values, dots indicate mean values, vertical lines 90% intervalls."
    , y = expression(paste("item-level classification uncertainty ", SD[i]))
    , x = expression(n[i])
    # , caption = expression(paste(
    #   "classification uncertainty measured as ",
    #   SD[i]==frac(1, D)*sum((c[it]-bar(c)[i])^2, t==1, D),
    #   " where D = 3 x 1000 is the number of estimates aggregated across chains and iterations"
    #   # " where D = 3 x 1000 is the number of estimates aggregated across chains and iterations, and ",
    #   # bar(c)[i]==frac(1, D)*sum(c[it], t==1, D)
    # ))
  ) + 
  theme(plot.caption = element_text(hjust = 0))


# Group by Average Rate of Change in SD on s \in [3,10] 
aroc <- function(data, xmin = 1, xmax = 8){
  fit <- lm(sd ~ poly(s, 2), data = data)
  
  tibble(s = c(xmin, xmax)) %>% 
    cbind(y_hat = predict.lm(fit, newdata = .)) %>% 
    summarise_all(funs(diff)) %>% 
    summarise(aroc = y_hat/s) %>% 
    as.numeric()
}

post_cs_sd_aroc <- post_cs_sds %>% 
  group_by(item) %>% 
  do(aroc = aroc(.)) %>% 
  unnest() %>% 
  mutate(aroc_quantile = findInterval(aroc, quantile(aroc, probs = seq(.1,.9,.06))))


post_cs_sd_adsd <- post_cs_sds %>% 
  group_by(item) %>% 
  # filter(item == 1) %>%
  summarise(avg_delta_sd = mean(sd - lag(sd, order_by = s), na.rm = TRUE)) %>% 
  mutate(adsd_quantile = findInterval(avg_delta_sd, quantile(avg_delta_sd, probs = seq(.1,.9,.06))))


post_cs_sd_dpn <- post_cs_sds %>% 
  group_by(item) %>% 
  # filter(item == 1) %>%
  summarise(avg_delta_is_neg = mean(sd - lag(sd, order_by = s), na.rm = TRUE) < 0)

library(viridis)

post_cs_sds %>% 
  left_join(post_cs_sd_aroc) %>% 
  # filter(s == 1) %>% 
  ggplot(aes(x = jitter(s, amount = .3), y = sd, color = factor(aroc_quantile))) +
  geom_point(size = 1, alpha = .5) +
  scale_x_continuous(limits = c(0,9)) +
  scale_color_viridis_d(direction = -1, option = "D")



post_cs_sds %>% 
  left_join(post_cs_sd_aroc) %>% 
  # filter(s == 1) %>%
  ggplot(aes(x = sd, fill = factor(aroc_quantile))) + 
  geom_density(color = NA, alpha = .5) + 
  scale_fill_viridis_d(direction = -1, option = "D") +
  # scale_fill_brewer(direction = -1, palette = "RdYlGn") +
  coord_flip() +
  facet_grid(cols = vars(s+2), scales = "free_x") + 
  theme_bw()
  

post_cs_sd_adsd %>% 
  filter(adsd_quantile > 11) %>% 
  arrange(desc(avg_delta_sd))

post_cs_sds %>% 
  left_join(post_cs_sd_adsd) %>% 
  filter(adsd_quantile > 11) %>% 
  ggplot(aes(x = sd, fill = factor(adsd_quantile))) + 
  geom_density(color = NA, alpha = 0.5) + 
  # scale_fill_viridis_d(direction = 1, option = "C") +
  scale_fill_brewer(direction = -1, palette = "RdYlGn") +
  coord_flip() +
  # facet_grid(cols = vars(s+2), scales = "free_x") + 
  facet_grid(cols = vars(s+2), scales = "fixed") + 
  theme_bw()

post_cs_sds %>% 
  left_join(post_cs_sd_dpn) %>% 
  ggplot(aes(x = sd, fill = factor(avg_delta_is_neg))) + 
  geom_density(color = NA, alpha = 0.5) + 
  scale_fill_viridis_d(direction = -1, option = "D") +
  # scale_fill_brewer(direction = -1, palette = "RdYlGn") +
  coord_flip() +
  # facet_grid(cols = vars(s+2), scales = "free_x") + 
  facet_grid(cols = vars(s+2), scales = "fixed") + 
  theme_bw()
  

post_cs_sds %>% 
  left_join(post_cs_sd_aroc) %>% 
  ggplot(
    aes(
      x = s
      , y = sd
      , group = interaction(s,aroc_quantile)
      , fill = factor(aroc_quantile)
    )
  ) + 
  geom_violin(color = NA, alpha = .5, position="identity") +
  scale_fill_viridis_d(direction = 1, option = "C") +
  theme_bw()

  
  post_cs_sds %>% 
  left_join(post_cs_sd_aroc) %>% 
  group_by(s, aroc_quantile) %>% 
  summarise(
    mean = mean(sd)
    , q5 = quantile(sd, .05)
    , q95 = quantile(sd, .95)
  ) %>% 
  ggplot(aes(x = s, y = mean)) +
  geom_violin(
    data = left_join(post_cs_sds, post_cs_sd_aroc)
    , mapping = aes(
      x = s
      , y = sd
      , group = interaction(s,aroc_quantile)
      , color = factor(aroc_quantile)
    )
    # , color = NA, fill = "grey", alpha = .75
  ) +
  scale_y_continuous(limits = c(0,.6)) +
  scale_x_discrete(labels = 3:10) +
  geom_point() + 
  geom_linerange(aes(ymin = q5, ymax = q95)) +
  theme_bw() 





idxs <- map_df(.x = s_fitted, .id = "s", .f = function(s) s$samples)

foo2 <- idxs %>% 
  ungroup() %>% 
  left_join(sim_codings) %>% 
  mutate(s = as.integer(s))



lfoo2 <- split(x = foo2, f = factor(s))

d <- foo2[foo2$s == 7, ]

s <- s_fitted[[1]]

kalphas <- map_df(
  .x = s_fitted
  , .id = "s"
  , .f = function(s) {
    s$samples %>% 
      ungroup() %>% 
      left_join(sim_codings, by = c("item", "coder")) %>% 
      select(-true_class) %>% 
      spread(item, judgment, sep = "_") %>% 
      select(-coder) %>% 
      as.matrix() %>% 
      icr::krippalpha(
        data = .
        , metric = "nominal"
        , bootstrap = TRUE
        , nboot = 1000
        , cores = 2
      ) %>% 
      .$bootstraps %>% 
      tibble(kalpha = .)
}) %>% 
  mutate(s = as.integer(s))

kalphas %>% 
  group_by(s) %>% 
  summarise(
    mean = mean(kalpha)
    , q5 = quantile(kalpha, .05)
    , q95 = quantile(kalpha, .95)
  ) %>% 
  ggplot(aes(x = s, y = mean)) +
  geom_violin(
    data = kalphas
    , mapping = aes(x = s, y = kalpha, group = s)
    , color = NA
    , fill = "grey"
    , alpha= .75
  ) +
  geom_point() +
  geom_linerange(aes(ymin = q5, ymax = q95)) + 
  theme_bw() + 
  labs(
    title = expression(paste(
      "Intercoder reliability as function of ", n[i],", the number of judgments aggregated per item."
    ))
    , subtitle = "Shaded violins summarize 1000 bootstrapped values, dots indicate mean values, vertical lines 90% intervalls."
    , x = expression(n[i])
    , y = expression(paste("Krippendorff's ", alpha))
  ) + 
  scale_x_continuous(breaks = 1:8, labels = 3:10)
  
