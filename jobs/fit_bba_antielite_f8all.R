

file_path <- file.path("~", "switchdrive", "Documents", "work", "phd", "methods", "crowd-sourced_annotation")


library(dplyr)
library(purrr)
library(tidyr)
library(data.table)
library(stringr)
library(rjags)

helpers <- c(
  "get_mcmc_estimates.R"
  , "simulate_binary_codings.R"
  , "transform_betabin_fit_posthoc.R"
)

{sapply(file.path(file_path, "code", "R", helpers), source); NULL}

# set seed
this_seed <- 1234
set.seed(this_seed)



d <- data.table::fread(file.path(file_path, "exdata", "f1317730.csv")) %>% 
  tbl_df()

names(d) <- str_replace(names(d), "^_", "f8_")

table(d$elites_gold)



# gold judgments ----
gold_judgments <- d %>% 
  filter(elites_gold != "") %>% 
  mutate(judgment = ifelse(elites == "yes", 1, -1)) %>% 
  group_by(f8_unit_id) %>% 
  summarize(
    n_pos = sum(elites == "yes")
    , n_neg = sum(elites == "no")
    , wvote = mean(judgment*f8_trust)
  ) %>% 
  mutate(
    n_judgments = n_pos + n_neg
    , decision = paste0(n_pos, ":", n_neg)
    , prop_judgments = n_pos/n_judgments
    , label = ifelse(wvote >= 0, "yes", "no")
    , wvote = (wvote + 1)/2
  )

n_gold <- nrow(gold_judgments)

gold_judgments_by_class <- d %>% 
  filter(elites_gold != "") %>% 
  group_by(f8_unit_id, elites_gold, elites) %>% 
  summarise(n_judgments = n_distinct(f8_worker_id)) %>% 
  group_by(f8_unit_id, elites_gold) %>% 
  mutate(prop_judgment = n_judgments/sum(n_judgments))


ggplot(gold_judgments_by_class
  , aes(
    y = n_judgments
    , x = as.factor(f8_unit_id)
    , fill = elites
    , group = elites_gold
  )
) +
  geom_bar(
    position = "stack"
    , stat = "identity"
    , alpha = .66
  ) +
  facet_grid(~elites_gold) + 
  coord_flip() + 
  labs(
    x = NULL
    , y = NULL
    , fill = "judged"
  )
    

# free judgments ----
free_judgments <- d %>% 
  filter(elites_gold == "") %>% 
  mutate(judgment = ifelse(elites == "yes", 1, -1)) %>% 
  group_by(f8_unit_id) %>% 
  summarize(
    n_pos = sum(elites == "yes")
    , n_neg = sum(elites == "no")
    , wvote = mean(judgment*f8_trust)
  ) %>% 
  mutate(
    n_judgments = n_pos + n_neg
    , decision = paste0(n_pos, ":", n_neg)
    , label = ifelse(wvote >= 0, "yes", "no")
    , wvote = (wvote + 1)/2
  )

n_free <- nrow(free_judgments)

sum_free_judgments <- free_judgments %>% 
  group_by(decision, label) %>% 
  summarize(n = n())

sum_free_judgments %>% 
  ggplot(aes(x = decision, y = n, fill = label)) + 
    geom_bar(
      position = "stack"
      , stat = "identity"
      , alpha = .66
    ) + 
    coord_flip() +
    labs(
      x = "Decision (yes:no)"
      , y = "Number of instances"
    )


free_judgments %>% 
  ggplot(aes(x = wvote, fill = label)) +
    geom_histogram(alpha = .75, color = NA) +
    facet_grid(~decision, scales = "free") + 
    labs(x = "Weighted mean of votes (yes = 1, no = -1)", y = NULL)



# fit BBA to test set (incl all judgments of gold instances) ----

samp_insts <- sample(free_judgments$f8_unit_id, size = floor(n_free*.1))

test_codings <- rbind(
  filter(d, elites_gold != "")
  , filter(d, elites_gold == "", f8_unit_id %in% samp_insts)
) %>% 
  select(f8_unit_id, f8_worker_id, elites, f8_golden, f8_trust) %>% 
  # reshuffle
  mutate(
    # judgment index
    index = row_number()
    # item index
    , item = group_indices(., f8_unit_id)
    # coder index
    , coder = group_indices(., f8_worker_id)
    # populism indiactors
    , judgment = as.integer(elites == "yes")
  ) %>% 
  sample_frac(1)

# check
test_codings %>% 
  filter(f8_golden) %>% 
  .$item %>% 
  unique() %>% 
  length() == n_gold




(n_judgments <- nrow(test_codings))
(n_coders <- length(unique(test_codings$coder)))
(n_items <- length(unique(test_codings$item)))

mcmc_data_test <- test_codings %>% 
  select(index, item, coder, judgment) %>% 
  get_codings_data()



model_file_path <- file.path(file_path, "models", "beta-binomial_by_annotator.jags")


n_chains = 3
n_burnin = 500
n_iter = 1000

elites_mcmc_model_test <- jags.model(
  file = model_file_path
  , data = mcmc_data_test
  , n.chains = n_chains
)


update(elites_mcmc_model_test, n_burnin)

load.module("dic")

elites_mcmc_fit_test <- coda.samples(
  elites_mcmc_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    , "alpha0", "beta0"
    , "alpha1", "beta1"
  )
  , n.iter = n_iter*20
  , thin = 20
)


bl_dic <- get_mcmc_estimates(fit.obj = elites_mcmc_fit_test, params = "deviance")

gridExtra::grid.arrange(
  ggplot(bl_dic, aes(x = iter, y = est, color = factor(chain))) +
    geom_line(alpha = .5) +
    theme_bw() +
    labs(
      title = "Deviance information criterion (DIC)"
      , subtitle = "Obtained for baseline model specification"
      , y = "DIC"
      , x = "Iteration"
      , color = "Chains:"
    ) +
    theme(
      legend.position = "bottom"
    )
  , ggplot(bl_dic, aes(x = est)) + 
    geom_density(fill = "grey", alpha = .75, color = NA) + 
    theme_bw() +
    labs(
      title = "Distribution of DIC after burn-in"
      , subtitle = "Obtained from 3 chains for 1000 iterations"
      , x = "DIC"
      , y = "Density"
    )
  , ncol = 2
)


coda::gelman.plot(elites_mcmc_fit_test[, "deviance"])

coda::autocorr.plot(elites_mcmc_fit_test[[1]][, "deviance"], auto.layout = F)
coda::autocorr.plot(elites_mcmc_fit_test[[1]][, "deviance"], auto.layout = F)
coda::autocorr.plot(elites_mcmc_fit_test[[3]][, "deviance"], auto.layout = F)


if (mean(elites_mcmc_fit_test[,"pi"][[1]]) > .5){
  elites_mcmc_fit_test_t <- transform_betabin_fit_posthoc(elites_mcmc_fit_test)
}

# prevalence ----
bl_pi <- get_mcmc_estimates(fit.obj = elites_mcmc_fit_test_t, params = "pi")

gridExtra::grid.arrange(
  bl_pi %>% 
    ggplot(aes(x = iter, y = est, color = factor(chain))) +
    geom_line(alpha = .5) +
    theme_bw() +
    labs(
      title = expression(paste("Posterior estimates of ", pi))
      , subtitle = "3 chains obtained from BBA model fitted to simulated data"
      , y = expression(pi)
      , x = "Iteration"
      , color = "Chains:"
    ) +
    theme(
      legend.position = "bottom"
      , axis.title.y = element_text(angle = 0, vjust = .5)
    )
  , bl_pi %>% 
    ggplot(aes(x = est)) + 
    # geom_density(fill = "grey", alpha = .75, color = NA) + 
    geom_histogram(fill = "grey", alpha = .75, color = NA, bins = 200) + 
    scale_x_continuous(limits = c(0,1)) +
    geom_vline(xintercept = pi, color = "red") +
    theme_bw() +
    labs(
      title = expression(paste("Distribution of posterior estimates of ", pi))
      , subtitle = "Red vertical line marks simulated parameter value"
      , caption = "Obtained from BBA model fitted to simulated data"
      , x = expression(pi)
      , y = "Density"
    ) +
    theme(plot.caption = element_text(hjust = 0))
  , ncol = 2
)

# abilities ----
post_thetas <- get_mcmc_estimates(fit.obj = elites_mcmc_fit_test_t, params = "^theta(0|1)\\[\\d+\\]$")

post_thetas %>% 
  group_by(parameter) %>% 
  summarise(
    mean = mean(est)
    , q05 = quantile(est, .05)
    , q95 = quantile(est, .95)
  ) %>% 
  mutate(
    coder = as.factor(gsub(".+\\[(\\d+)\\]$", "\\1", parameter))
    , parameter = gsub("^theta(\\d).+$", "theta[\\1]", parameter, perl = TRUE)
  ) %>% 
  {right_join(
    filter(., parameter == "theta[0]") %>% 
      rename_at(2:4, paste0, "_theta0") %>% 
      select(-parameter)
    , filter(., parameter == "theta[1]") %>% 
      rename_at(2:4, paste0, "_theta1") %>% 
      select(-parameter)
    , by = "coder"
  )} %>% 
  ggplot(
    aes(
      x = mean_theta0, xmin = q05_theta0, xmax = q95_theta0,
      y = mean_theta1, ymin = q05_theta1, ymax = q95_theta1,
      group = coder, 
      alpha = .5
    )
  ) +
  # geom_errorbarh() +
  # geom_errorbar() +
  geom_point() +
  scale_x_continuous(limits = 0:1) +
  scale_y_continuous(limits = 0:1) +
  guides(alpha = FALSE) +
  labs(
    x = expression(theta[0])
    , y = expression(theta[1])
  ) +
  theme(
    axis.title.y = element_text(angle = 0, vjust = .5)
  )

post_thetas %>%
  group_by(parameter) %>%
  summarize(est = mean(est)) %>%
  mutate(param = if_else(grepl("^theta0", parameter), "bar(theta)[0]", "bar(theta)[1]")) %>%
  ggplot(aes(x = est)) +
  geom_density(color = NA, fill = "grey", alpha = .75) +
  facet_grid(cols = vars(param), labeller = label_parsed) +
  labs(x = NULL, y = NULL)


mean_thetas <- post_thetas %>%
  group_by(parameter) %>% 
  summarize(est = mean(est)) %>% 
  mutate(
    coder = as.integer(gsub(".+\\[(\\d+)\\]$", "\\1", parameter))
    , param = gsub("^theta(\\d).+", "theta[\\1]", parameter)
  ) %>% 
  select(-parameter) 
  

mean_thetas_trust <- test_codings %>% 
  select(coder, f8_worker_id, f8_trust) %>% 
  unique() %>% 
  left_join(mean_thetas)


ggplot(mapping = aes(x = f8_trust, y = est, alpha = .75)) +
  geom_abline(slope = 1, intercept = 0, show.legend = FALSE) +
  geom_point(data = mean_thetas_trust %>% filter(param == "theta[0]")) +
  geom_point(data = mean_thetas_trust %>% filter(param == "theta[1]")) +
  scale_x_continuous(limits = c(.4, 1)) +
  scale_y_continuous(limits = c(.4, 1)) +
  facet_grid(~param, labeller = label_parsed) + 
  guides(alpha = FALSE) +
  labs(
    x = "FigureEight coder trust score"
    , y = "Mean posterior estimate"
  )

Specificties systematically higher than trust scores
Sensitivities systematically lower than trust scores

# ability hyperparameters ----
post_hyperpars <- get_mcmc_estimates(elite_mcmc_fit_t, "^(alpha|beta)(0|1)$")

post_hyperpars %>%
  mutate(
    param = sub(".+(\\d)$", "theta[\\1]", parameter)
    , parameter = sub("\\d$", "", parameter)
  ) %>%
  spread(parameter, est) %>%
  mutate(
    iqscale = 1/(alpha + beta)^2
    , mean = alpha/(alpha + beta)
  ) %>%
  ggplot(aes(x = mean, y = iqscale)) +
  geom_point(alpha = .5, size = .5) +
  facet_grid(~param, labeller = label_parsed) +
  labs(
    x = expression(paste("Mean, ", alpha[.]/(alpha[.]+beta[.])))
    , y = expression(paste("Inverse of squared scale, ", 1/(alpha[.]+beta[.])^2))
  )



# posterior class labels ---- 

post_c <- get_mcmc_estimates(fit.obj = elites_mcmc_fit_test_t, params = "c\\[\\d+\\]")



gold_agreement <- test_codings %>% 
  filter(f8_golden) %>% 
  left_join(select(d, f8_unit_id, f8_worker_id, elites_gold)) %>% 
  select(item, elites_gold) %>% 
  mutate(gold = as.integer(elites_gold == "yes")) %>% 
  unique() %>% 
  left_join(
    post_c %>% 
      mutate(item = as.integer(str_extract(parameter, "\\d+"))) %>% 
      group_by(item) %>% 
      summarise(label = median(est))
  ) %>% 
  mutate(agrees = gold == label)


# accuracy on gold instances    
gold_agreement %>% 
  group_by(agrees) %>% 
  summarise(n = n(), n_gold = n_gold)

overall_agreement <- test_codings %>% 
  left_join(
    rbind(
      gold_judgments %>% select(f8_unit_id, label)
      , free_judgments %>% select(f8_unit_id, label)
    )
  ) %>% 
  select(item, f8label = label) %>% 
  mutate(f8label = as.integer(f8label == "yes")) %>% 
  unique() %>% 
  left_join(
    post_c %>% 
      mutate(item = as.integer(str_extract(parameter, "\\d+"))) %>% 
      group_by(item) %>% 
      summarise(label = median(est))
  ) %>% 
  mutate(agrees = f8label == label)

overall_agreement %>% 
  group_by(agrees) %>% 
  summarise(n = n())


# Fit BBA model to all judgments ----

all_codings <- d %>% 
  select(f8_unit_id, f8_worker_id, elites, f8_golden, f8_trust) %>% 
  # reshuffle
  mutate(
    # judgment index
    index = row_number()
    # item index
    , item = group_indices(., f8_unit_id)
    # coder index
    , coder = group_indices(., f8_worker_id)
    # populism indiactors
    , judgment = as.integer(elites == "yes")
  ) 

(n_judgments <- nrow(test_codings))
(n_coders <- length(unique(test_codings$coder)))
(n_items <- length(unique(test_codings$item)))

mcmc_data <- all_codings %>% 
  select(index, item, coder, judgment) %>% 
  get_codings_data()

elites_mcmc_model <- jags.model(
  file = model_file_path
  , data = mcmc_data
  , n.chains = n_chains
)

update(elites_mcmc_model, n_burnin)


elites_mcmc_fit <- coda.samples(
  elites_mcmc_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    , "alpha0", "beta0"
    , "alpha1", "beta1"
  )
  , n.iter = n_iter*10
  , thin = 10
)


bl_dic <- get_mcmc_estimates(fit.obj = elites_mcmc_fit, params = "deviance")

gridExtra::grid.arrange(
  ggplot(bl_dic, aes(x = iter, y = est, color = factor(chain))) +
    geom_line(alpha = .5) +
    theme_bw() +
    labs(
      title = "Deviance information criterion (DIC)"
      , subtitle = "Obtained for baseline model specification"
      , y = "DIC"
      , x = "Iteration"
      , color = "Chains:"
    ) +
    theme(
      legend.position = "bottom"
    )
  , ggplot(bl_dic, aes(x = est)) + 
    geom_density(fill = "grey", alpha = .75, color = NA) + 
    theme_bw() +
    labs(
      title = "Distribution of DIC after burn-in"
      , subtitle = "Obtained from 3 chains for 1000 iterations"
      , x = "DIC"
      , y = "Density"
    )
  , ncol = 2
)


coda::gelman.plot(elites_mcmc_fit[, "deviance"])

coda::autocorr.plot(elites_mcmc_fit[[1]][, "deviance"], auto.layout = F)
coda::autocorr.plot(elites_mcmc_fit[[1]][, "deviance"], auto.layout = F)
coda::autocorr.plot(elites_mcmc_fit[[3]][, "deviance"], auto.layout = F)


saveRDS(list(elites_mcmc_fit), file.path(file_path, "fits", "betabinom_by_annotator_allf8.RData"))




dt <- data.table::fread(file.path(file_path, "exdata", "f1382750.csv")) %>% 
  tbl_df()


names(dt) <- str_replace(names(dt), "^_", "f8_")

# inspect n of n_i
dt %>% 
  group_by(f8_golden, f8_unit_id) %>% 
  summarise(n_i = n_distinct(f8_worker_id)) %>% 
  group_by(f8_golden, n_i) %>% 
  summarise(n = n_distinct(f8_unit_id))

# gold judgments ----
gold_judgments <- dt %>% 
  filter(elites_gold != "") %>% 
  mutate(judgment = ifelse(elites == "yes", 1, -1)) %>% 
  group_by(f8_unit_id) %>% 
  summarize(
    n_pos = sum(elites == "yes")
    , n_neg = sum(elites == "no")
    , wvote = mean(judgment*f8_trust)
  ) %>% 
  mutate(
    n_judgments = n_pos + n_neg
    , decision = paste0(n_pos, ":", n_neg)
    , prop_judgments = n_pos/n_judgments
    , label = ifelse(wvote >= 0, "yes", "no")
    , wvote = (wvote + 1)/2
  )

n_gold <- nrow(gold_judgments)

gold_judgments_by_class <- dt %>% 
  filter(elites_gold != "") %>% 
  group_by(f8_unit_id, elites_gold, elites) %>% 
  summarise(n_judgments = n_distinct(f8_worker_id)) %>% 
  group_by(f8_unit_id, elites_gold) %>% 
  mutate(prop_judgment = n_judgments/sum(n_judgments))


ggplot(gold_judgments_by_class
       , aes(
         y = n_judgments
         , x = as.factor(f8_unit_id)
         , fill = elites
         , group = elites_gold
       )
) +
  geom_bar(
    position = "stack"
    , stat = "identity"
    , alpha = .66
  ) +
  facet_grid(~elites_gold) + 
  coord_flip() + 
  labs(
    x = NULL
    , y = NULL
    , fill = "judged"
  )


# free judgments ----
free_judgments <- dt %>% 
  filter(elites_gold == "") %>% 
  mutate(judgment = ifelse(elites == "yes", 1, -1)) %>% 
  group_by(f8_unit_id) %>% 
  summarize(
    n_pos = sum(elites == "yes")
    , n_neg = sum(elites == "no")
    , wvote = mean(judgment*f8_trust)
  ) %>% 
  mutate(
    n_judgments = n_pos + n_neg
    , decision = paste0(n_pos, ":", n_neg)
    , label = ifelse(wvote >= 0, "yes", "no")
    , wvote = (wvote + 1)/2
  )

n_free <- nrow(free_judgments)

sum_free_judgments <- free_judgments %>% 
  group_by(decision, label) %>% 
  summarize(n = n())

sum_free_judgments %>% 
  ggplot(aes(x = decision, y = n, fill = label)) + 
  geom_bar(
    position = "stack"
    , stat = "identity"
    , alpha = .66
  ) + 
  coord_flip() +
  labs(
    x = "Decision (yes:no)"
    , y = "Number of instances"
  )


free_judgments %>% 
  ggplot(aes(x = wvote, fill = label)) +
  # geom_density(alpha = .75, color = NA) +
  geom_histogram(alpha = .75, color = NA, bins = 15) +
  facet_grid(~decision, scales = "free") + 
  labs(x = "Weighted mean of votes (yes = 1, no = -1)", y = NULL)


# Fit BBA model to all judgments ----

all_codings <- dt %>% 
  select(f8_unit_id, f8_worker_id, elites, f8_golden, f8_trust) %>% 
  # reshuffle
  mutate(
    # judgment index
    index = row_number()
    # item index
    , item = group_indices(., f8_unit_id)
    # coder index
    , coder = group_indices(., f8_worker_id)
    # populism indiactors
    , judgment = as.integer(elites == "yes")
  ) 

(n_judgments <- nrow(all_codings))
(n_coders <- length(unique(all_codings$coder)))
(n_items <- length(unique(all_codings$item)))

mcmc_data <- all_codings %>% 
  select(index, item, coder, judgment) %>% 
  get_codings_data()

elites_mcmc_model <- jags.model(
  file = model_file_path
  , data = mcmc_data
  , n.chains = n_chains
)

update(elites_mcmc_model, n_burnin)

load.module("dic")

elites_mcmc_fit <- coda.samples(
  elites_mcmc_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    , "alpha0", "beta0"
    , "alpha1", "beta1"
  )
  , n.iter = n_iter*10
  , thin = 10
)


bl_dic <- get_mcmc_estimates(fit.obj = elites_mcmc_fit, params = "deviance")


gridExtra::grid.arrange(
  ggplot(bl_dic, aes(x = iter, y = est, color = factor(chain))) +
    geom_line(alpha = .5) +
    theme_bw() +
    labs(
      title = "Deviance information criterion (DIC)"
      , subtitle = "Obtained for baseline model specification"
      , y = "DIC"
      , x = "Iteration"
      , color = "Chains:"
    ) +
    theme(
      legend.position = "bottom"
    )
  , ggplot(bl_dic, aes(x = est)) + 
    geom_density(fill = "grey", alpha = .75, color = NA) + 
    theme_bw() +
    labs(
      title = "Distribution of DIC after burn-in"
      , subtitle = "Obtained from 3 chains for 1000 iterations"
      , x = "DIC"
      , y = "Density"
    )
  , ncol = 2
)


coda::gelman.plot(elites_mcmc_fit[, "deviance"])

coda::autocorr.plot(elites_mcmc_fit[[1]][, "deviance"], auto.layout = F)
coda::autocorr.plot(elites_mcmc_fit[[1]][, "deviance"], auto.layout = F)
coda::autocorr.plot(elites_mcmc_fit[[3]][, "deviance"], auto.layout = F)





# join codings of test set instances ----

glimpse(dt)
all(dt$f8_unit_id %in% d$f8_unit_id)
any(dt$f8_unit_id %in% d$f8_unit_id)

all(dt$id_str %in% d$id_str)
any(dt$id_str %in% d$id_str)

these_cols <- c(
  "f8_golden"
  , "f8_worker_id"
  , "elites"
  , "elites_gold"
  , "id_str"
)

test_instance_id_strs <- unique(dt$id_str)



d_s <- d %>% 
  select(!!these_cols) %>% 
  # keep only judgments of instances also in test set
  filter(id_str %in% unique(dt$id_str)) %>% 
  mutate(round = 1L)

# get IDs of workers that judged instances items in first run 
these_workers <- d_s %>% 
  filter(!f8_golden) %>% 
  .$f8_worker_id


dta <- rbind(
  # keep only gold judgments of workers that also judged test instances
  filter(d_s, f8_worker_id %in% these_workers)
  , select(dt, !!these_cols) %>% mutate(round = 2L)
) %>% 
  arrange(id_str)


# inspect n of n_i
dta %>% 
  # group_by(f8_golden, id_str) %>% 
  # filter(max(round) == 1)
  group_by(f8_golden, id_str) %>% 
  summarise(n_i = n_distinct(f8_worker_id)) %>% 
  group_by(f8_golden, n_i) %>% 
  summarise(n = n_distinct(id_str))

flag <- dta %>% 
  group_by(f8_golden, id_str) %>% 
  mutate(
    n_i = n_distinct(f8_worker_id)
    # , min_round = min(round) 
    # , max_round = max(round) 
  )

# items that are only coded three times are all from 1st round
flag %>% 
  filter(n_i == 3) %>% 
  .$round %>% 
  unique()

# 2 instances are judged by only three coders
flag %>% 
  filter(n_i == 9) %>% 
  select(id_str) %>% 
  unique()

# this is due to two coders coding in both rounds
dta %>% 
  # group_by(f8_golden, id_str) %>% 
  # filter(max(round) == 1)
  group_by(f8_golden, id_str) %>% 
  summarise(
    n_coders = n_distinct(f8_worker_id)
    , n_i = n()
  ) %>% 
  group_by(f8_golden, n_coders, n_i) %>% 
  summarise(n = n_distinct(id_str)) %>% 
  # filter(n_coders == 9)
  filter(n_coders != n_i) 

# fit models to n_i = 3...10 ----

# clean
dtac <- dta %>% 
  group_by(f8_golden, id_str) %>% 
  mutate(
    n_coders = n_distinct(f8_worker_id)
    , n_i = n()
  ) %>% 
  filter(n_coders > 9) %>% 
  # filter(n_coders != n_i) %>% 
  # arrange(id_str, f8_worker_id) %>% 
  group_by(id_str, f8_worker_id) %>% 
  mutate(n_judged_i = n()) %>% 
  # filter(n_judged_i > 1) %>% 
  mutate(n_dist_judgments = n_distinct(elites)) %>% 
  # filter(n_dist_judgments > 1) %>% 
  filter( # ... keep all judgments ...
    # ... of coders that judged each instance only once ...
    n_judged_i == 1 
    # ... or that multiply coded an instance with the same judgment ... 
    | n_dist_judgments == 1 
    # .... and otherwise keep judgment from 2nd round
    | round == 2
  ) 
  

dtacc <- dtac %>% 
  select(!!these_cols) %>% 
  # removes superflous identical judgment of coders who judged an instance twice 
  unique()

# inspect
dtacc %>% 
  group_by(f8_golden, id_str) %>% 
  summarise(
    n_coders = n_distinct(f8_worker_id)
    , n_i = n()
  ) %>% 
  group_by(f8_golden, n_coders, n_i) %>% 
  summarise(n = n_distinct(id_str)) %>% 
  filter(n_coders != n_i) 
  # was successful

dtacc %>% 
  group_by(f8_golden, id_str) %>% 
  summarise(
    n_coders = n_distinct(f8_worker_id)
    , n_i = n()
  ) %>% 
  group_by(f8_golden, n_coders, n_i) %>% 
  summarise(n = n_distinct(id_str))


valid_workers <- dtacc %>% 
  filter(!f8_golden) %>% 
  .$f8_worker_id %>% 
  unique()
  
  
dtr <- dtacc %>% 
  # get rid of judgments of worker that havn't valid (n_i=10) instances
  filter(f8_worker_id %in% valid_workers)

# resample 
set.seed(1234)

these_instances <- unique(dt[!dt$f8_golden, "id_str"])[[1]]
these_gs <- unique(dt[dt$f8_golden, "id_str"])[[1]]

test_codings <- dtr %>% 
  filter(id_str %in% these_instances) %>% 
  group_by(id_str) %>% 
  # reshuffle
  sample_n(10) %>% 
  # generate judgment id
  mutate(position = row_number()) %>% 
  ungroup() %>% 
  # create
  select(-f8_golden)

# one judgment was gold in 2nd round but not in 1st
test_codings %>% 
  group_by(id_str) %>% 
  summarise(
    n_coders = n_distinct(f8_worker_id)
    , n_i = n()
  ) %>% 
  group_by(n_coders, n_i) %>% 
  summarise(n = n_distinct(id_str))

# create sliced datasets
test_splits <- map(3:10, function(n, d = test_codings) {
  out <- d %>%
    filter(position <= n) %>%
    select(-position) %>%
    mutate() %>% 
    mutate(
      n_i = n
      # # judgment index
      # index = row_number()
      # item index
      , item = group_indices(., id_str)
      # coder index
      , coder = group_indices(., f8_worker_id)
      # populism indiactors
      , judgment = as.integer(elites == "yes")
    )
  
  list(data = out)
})


# fix number of chains
n_chains <- 3L
# read grid with burnin parameters
n_burnin <- read.csv(file.path(file_path, "config", "burnin_grid.csv"), row.names = 1)[[3]]
# read grid with thinning parameters
thin <- read.csv(file.path(file_path, "config", "thinning_grid.csv"), row.names = 1)[[3]]


for (n in 4:8) {
  
  # initialize ...
  model <- jags.model(
    file = model_file_path
    , data = get_codings_data(test_splits[[n]]$data)
    , n.chains = n_chains
  )
  
  # ... update ...
  update(model, n_burnin[n])
  
  # ..- and fit model
  test_splits[[n]]$fit <- coda.samples(
    model
    , variable.names = model_parameters
    , n.iter = 1000L*thin[n]
    , thin = thin[n]
  )
  
}

saveRDS(test_splits, file.path(file_path, "fits", "bba_antielite_testf8_vary_n_i.RData"))










