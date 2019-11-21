# setup ----

library(rjags)
library(purrr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggridges)

set.seed(1234)

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

# create JAGS model inputs ----

# data
codings_data <- list() 

codings_data$N <- length(unique(sim_codings$item))
codings_data$M <- length(unique(sim_codings$coder))
codings_data$J <- length(sim_codings$judgment)
codings_data$y <- sim_codings %>% 
  arrange(item, coder) %>% 
  select(item, coder, judgment)

# Fit baseline model ----

model_file_path <- "./models/beta-binomial_by_annotator.jags"
system(paste("less", model_file_path))

load.module("dic")

stime <- Sys.time()

baseline_model <- jags.model(
  file = model_file_path
  , data = codings_data
  , inits = init_vals
  , n.chains = 3
)

update(baseline_model, 500)


baseline_fit <- coda.samples(
  baseline_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    , "alpha0", "beta0"
    , "alpha1", "beta1"
  )
  , n.iter = 1000
  , thin = 1
)
etime <- Sys.time()

etime-stime

# save to disk
# list(
#   model = baseline_model
#   , fit = baseline_fit
#   , runtime = end-start
# ) %>% 
#   saveRDS("./fits/betabinom_by_annotator_baseline.RData")

# Inspect baseline fit ---- 

# nicely converged already only a couple of iterations
plot(baseline_fit[, "deviance"])
gelman.plot(baseline_fit[, "deviance"])

# very few autocorrelation
autocorr.plot(baseline_fit[, "deviance"], ask = F)

# Inspect baseline parameter estimates ----

# Prevalence

# chains mix nicely
plot(baseline_fit[, "pi"], main = expression(pi))
# nice and quick convergence
gelman.plot(baseline_fit[, "pi"], main = expression(pi))
# but posterior density has mass in region [.78, .81]
ggplot(
  tibble(x = unlist(baseline_fit[, "pi"]))
  , aes(x = x)
) + 
  geom_density(color = "grey", fill = "lightgrey", alpha = .75) + 
  scale_x_continuous(limits = c(0,1)) +
  geom_vline(xintercept = pi, color = "red") +
  theme_bw() +
  labs (
    title = expression(paste("Marginal posterior density of ", pi))
    , subtitle = "Red line indicates 'true' value"
    , y = "Density"
    , x = ""
  )


# Coder abilities

# since we have `r n_coder` posterior density distributions for each $theta_0, theta_1$,
# we'll omit visual inspection of chains and autocorrelation and instead rely on the DIC to judge convergence.

# put we can plot the posterior densities by coder and parameter 

# gives n_coders*n_chains*n_iterations = 60,000 rows
post_thetas <- map_df(1:n_coders, function(j) 
  tibble(
    coder = j
    , theta0 = unlist(baseline_fit[, sprintf("theta0[%s]", j)] )
    , theta1 = unlist(baseline_fit[, sprintf("theta1[%s]", j)] )
  )
)

# post_thetas %>% 
#   gather(parameter, value, -coder) %>% 
#   # filter(parameter == "theta1") %>%
#   ggplot(aes(y = value, x = factor(coder), group = coder)) +
#   geom_violin() +
#   facet_grid(rows = vars(parameter))

post_thetas %>% 
  gather(parameter, value, -coder) %>% 
  ggplot(
    aes(
      y = factor(coder)
      , x = value
      , group = coder
      # , fill = factor(parameter)
    )
  ) +
  geom_density_ridges(
    scale = .5
    , color = NA
  ) +
  geom_point(
    shape = 3
    , color = "red"
    , data = tibble(
      coder = rep(1:n_coders, 2)
      , parameter = rep(paste0("theta", 0:1), each = n_coders)
      , value = c(theta0, theta1)
    )
    , aes(y = factor(coder), x = value, group = coder)
  ) +
  scale_x_continuous(limits = c(0,1))+
  facet_grid(rows = vars(parameter)) + 
  labs(
    title = "Marginal posterior densities of ability parameters by coder"
    , subtitle = "Red crosses mark 'true' values"
    , x = ""
    , y = "Coder"
    , fill = "Parameter"
  ) + 
  coord_flip() +
  theme_bw()
# We can clearly see how ability parameter estimates are consistently flipped relative to the 'true' simulation parameter values
# just as in the case of prevalence

post_thetas %>% 
  gather(param, value, -coder) %>% 
  group_by(coder, param) %>% 
  summarise(
    mean = mean(value)
    , q5 = quantile(value, .05)
    , q95 = quantile(value, .95)
  ) %>% 
  arrange(param, coder) %>% 
  ungroup() %>% 
  mutate(true_val = c(sim_theta0, sim_theta1)) %>% 
  ggplot(aes(x = true_val, y = mean)) +
    # geom_smooth(method='lm', formula=y~x, se = FALSE) +
    geom_point(size = .5) +
    geom_linerange(aes(ymin = q5, ymax = q95), width = .1) +
    scale_x_continuous(limits = 0:1) +
    scale_y_continuous(limits = 0:1) +
    geom_abline(slope = 1, intercept = 0, size = .5, color = "grey") +
    facet_grid(~param) +
    theme_bw() + 
    labs(
      title = "Mean posterior coder abilities vs. 'true' parameter values."
      , subtitle = "Vertical lines overlaying mean values depict 90% credibility intervall, the diagonal indicates perfect correspondence."
      , x = "Simulated parameter value"
      , y = "Posterior value"
    )

# The data pulls ability parameters towards low values, but the informative beta priors pull in the opposite direction
# we can also see how the more informative (more dispersed) prior on $\theta_0$ ($f_{B}(40,8)$ results in higher dispersion 
# of posteriors 


# Item classes 

# gives n_items*n_chains*n_iterations = 1000*3*1000 = 3,000,000 rows
post_c <- map_df(1:n_items, function(i) 
  tibble(item = i, c = unlist(baseline_fit[, sprintf("c[%s]", i)]))
)


post_c %>% 
  group_by(item) %>% 
  summarise(pr_pos = sum(c)/n()) %>% 
  mutate(est_class = as.integer(pr_pos > .5)) %>% 
  left_join(tibble(item = 1:n_items, true_class = sim_item_classes)) %>% 
  mutate(aggreement = true_class == est_class) %>% 
  summarise(accuracy = sum(aggreement)/n())
# Finally, we can see that the model-based classification of items is almost perfectly contradicting true class memberships
# In sum, evaluation of the baseline model suggests that we would want to constraint the model to induce convergence on true values and not mirroring ones

# somehow, all estimates are flipped

# posterior alpha0
post_alpha0 <- unlist(baseline_fit[, "alpha0"])
post_alpha0_dens <- density(post_alpha0)
plot(
  post_alpha0_dens
  , main = expression(paste("Marginal posterior density of",~alpha[0],~"and mean value."))
  , xlab = expression(alpha[0])
  , axes = FALSE
)
axis(1)
axis(2, las = 2)
y1 <- post_alpha0_dens$y[which(abs(post_alpha0_dens$x - mean(post_alpha0)) == min(abs(post_alpha0_dens$x - mean(post_alpha0))))]
segments(
  x0 = mean(post_alpha0)
  , x1 = mean(post_alpha0)
  , y0 = 0
  , y1 = y1
)
text(x = mean(post_alpha0), y = y1, round(mean(post_alpha0),3), pos = 3, offset = .90)


post_beta0 <- unlist(baseline_fit[, "beta0"])
post_beta0_dens <- density(post_beta0)
plot(
  post_beta0_dens
  , main = expression(paste("Marginal posterior density of",~beta[0],~"and mean value."))
  , xlab = expression(beta[0])
  , axes = FALSE
)
axis(1)
axis(2, las = 2)
y1 <- post_beta0_dens$y[which(abs(post_beta0_dens$x - mean(post_beta0)) == min(abs(post_beta0_dens$x - mean(post_beta0))))]
segments(
  x0 = mean(post_beta0)
  , x1 = mean(post_beta0)
  , y0 = 0
  , y1 = y1
)
text(x = mean(post_beta0), y = y1, round(mean(post_beta0),3), pos = 3, offset = .90)


post_p <- rbeta(300, post_alpha0, post_beta0)

plot(density(post_p))
hist(post_p)
rb
# When we take the sampled alpha0 and beta0 parameters and use them pairwise to randomly generate probabilites from a Beta distribution,
# we get the following distribution of probabilities




post_alpha1 <- unlist(baseline_fit[, "alpha1"])
post_alpha1_dens <- density(post_alpha1)
plot(
  post_alpha1_dens
  , main = expression(paste("Posterior density of",~alpha[1],~"and mean value."))
  , xlab = expression(alpha[1])
  , axes = FALSE
)
axis(1)
axis(2, las = 2)
y1 <- post_alpha1_dens$y[which(abs(post_alpha1_dens$x - mean(post_alpha1)) == min(abs(post_alpha1_dens$x - mean(post_alpha1))))]
segments(
  x0 = mean(post_alpha1)
  , x1 = mean(post_alpha1)
  , y0 = 0
  , y1 = y1
)
text(x = mean(post_alpha1), y = y1, round(mean(post_alpha1),3), pos = 3, offset = .10)
# less well in recovering the simulation parameter alpha0, which was set to 


post_beta1 <- unlist(baseline_fit[, "beta1"])
post_beta1_dens <- density(post_beta1)
plot(
  post_beta1_dens
  , main = expression(paste("Posterior density of",~beta[1],~"and mean value."))
  , xlab = expression(beta[1])
  , axes = FALSE
)
axis(1)
axis(2, las = 2)
y1 <- post_beta1_dens$y[which(abs(post_beta1_dens$x - mean(post_beta1)) == min(abs(post_beta1_dens$x - mean(post_beta1))))]
segments(
  x0 = mean(post_beta1)
  , x1 = mean(post_beta1)
  , y0 = 0
  , y1 = y1
)
text(x = mean(post_beta1), y = y1, paste(round(mean(post_beta1),3), sprintf("(%s)", beta1)), pos = 3, offset = .20)


post_p <- rbeta(300, post_alpha0, post_beta0)

plot(density(post_p))
hist(post_p)
rb



plot(baseline_fit[, grepl"theta0[1]"])
gelman.plot(baseline_fit[, "theta0[1]"])

plot(baseline_fit[, "c[10]"])
gelman.plot(baseline_fit[, "c[10]"])


curve(dbeta(x, 2, 8), 0, 1, main = expression(x%~%f[B](2,8)))

curve(dbeta(x*2, 2, 2), 0, 1, ylab = main = expression(x%~%f[B](2,2)/2))
curve(dbeta(x*2, 1.5, 1.5), 0, 1, main = expression(x%~%f[B](1.5,1.5)/2))
curve(dbeta(x*2, 1.1, 1.1), 0, 1, main = expression(x%~%f[B](1.1,1.1)/2))

codings$true_class[1]

sum_baseline_fit <- summary(baseline_fit)





# Fit baseline model with starting values ----


# initial values

# create intial values for theta0 and theta1 parameters

init_vals <- lapply(1:3, function(chain) {
  
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

model_file_path <- "./models/beta-binomial_by_annotator.jags"
system(paste("less", model_file_path))

load.module("dic")

baseline_w_svals_model <- jags.model(
  file = model_file_path
  , data = codings_data
  , inits = init_vals
  , n.chains = 3
)

update(baseline_w_svals_model, 500)

baseline_w_svals_fit <- coda.samples(
  baseline_w_svals_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    , "alpha0", "beta0"
    , "alpha1", "beta1"
  )
  , n.iter = 1000
  , thin = 1
)

# save to disk
list(
  model = baseline_w_svals_model
  , fit = baseline_w_svals_fit
  , runtime = end-start
) %>%
  saveRDS("./fits/betabinom_by_annotator_w_start_vals.RData")

plot(baseline_w_svals_fit[, "deviance"])
gelman.plot(baseline_w_svals_fit[, "deviance"])

plot(baseline_w_svals_fit[, "pi"])
gelman.plot(baseline_w_svals_fit[, "pi"])





# Fit adapted model ----- 

# Fit simplified model (multilevel structure of pi omitted) ----


model_file_path <- "./models/beta-binomial_by_annotator_simple.jags"
system(paste("less", model_file_path))

stime <- Sys.time()

simplified_model <- jags.model(
  file = model_file_path
  , data = codings_data
  , inits = init_vals
  , n.chains = 3
)

update(simplified_model, 500)


simplified_fit <- coda.samples(
  simplified_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    # , "alpha0", "beta0"
    # , "alpha1", "beta1"
  )
  , n.iter = 1000
  , thin = 1
)
etime <- Sys.time()

etime-stime

# save to disk
list(
  model = simplified_model
  , fit = simplified_fit
  , runtime = etime-stime
) %>%
  saveRDS("./fits/betabinom_by_annotator_simplified.RData")

# Inspect fit 

# nicely converged already only a couple of iterations
plot(simplified_fit[, "deviance"])
gelman.plot(simplified_fit[, "deviance"])

# very few autocorrelation
autocorr.plot(simplified_fit[, "deviance"], ask = F)


# nicely converged already only a couple of iterations
plot(simplified_fit[, "pi"])
gelman.plot(simplified_fit[, "pi"])



# 


# Fit baseline model with starting values ----


# initial values

# create intial values for theta0 and theta1 parameters

init_vals <- lapply(1:3, function(chain) {
  
  out <- list()
  out[["theta0"]] <- rbeta(n_coders, 8, 2)
  out[["theta1"]] <- rbeta(n_coders, 8, 2)
  out[["pi"]] <- rbeta(1, 1, 8)
  out[[".RNG.name"]] <- "base::Wichmann-Hill"
  out[[".RNG.seed"]] <- 1234
  
  return(out)
})

# inspect
map(init_vals, "theta0")
map(init_vals, "theta1")
map(init_vals, "pi")

model_file_path <- "./models/beta-binomial_by_annotator.jags"
system(paste("less", model_file_path))

load.module("dic")

baseline_w_svals_model <- jags.model(
  file = model_file_path
  , data = codings_data
  , inits = init_vals
  , n.chains = 3
)

update(baseline_w_svals_model, 500)

baseline_w_svals_fit <- coda.samples(
  baseline_w_svals_model
  , variable.names = c(
    "deviance"
    , "pi"
    , "c"
    , "theta0", "theta1"
    , "alpha0", "beta0"
    , "alpha1", "beta1"
  )
  , n.iter = 1000
  , thin = 1
)

# save to disk
list(
  model = baseline_w_svals_model
  , fit = baseline_w_svals_fit
  , runtime = end-start
) %>%
  saveRDS("./fits/betabinom_by_annotator_w_start_vals.RData")

plot(baseline_w_svals_fit[, "deviance"])
gelman.plot(baseline_w_svals_fit[, "deviance"])

plot(baseline_w_svals_fit[, "pi"])
gelman.plot(baseline_w_svals_fit[, "pi"])


# Inspect baseline fit ---- 

# nicely converged already only a couple of iterations
plot(baseline_fit[, "deviance"])
gelman.plot(baseline_fit[, "deviance"])

# very few autocorrelation
par(mfrow = c(3,1))
autocorr.plot(baseline_fit[, "deviance"], ask = F)

# Inspect baseline parameter estimates ----

# Prevalence

# chains mix nicely
plot(baseline_fit[, "pi"], main = expression(pi))
# nice and quick convergence
gelman.plot(baseline_fit[, "pi"], main = expression(pi))
# but posterior density has mass in region [.78, .81]
ggplot(
  tibble(x = unlist(baseline_fit[, "pi"]))
  , aes(x = x)
) + 
  geom_density(color = "grey", fill = "lightgrey", alpha = .75) + 
  scale_x_continuous(limits = c(0,1)) +
  geom_vline(xintercept = pi, color = "red") +
  theme_bw() +
  labs (
    title = expression(paste("Marginal posterior density of ", pi))
    , subtitle = "Red line indicates 'true' value"
    , y = "Density"
    , x = ""
  )


# Coder abilities

# since we have `r n_coder` posterior density distributions for each $theta_0, theta_1$,
# we'll omit visual inspection of chains and autocorrelation and instead rely on the DIC to judge convergence.

# put we can plot the posterior densities by coder and parameter 

# gives n_coders*n_chains*n_iterations = 60,000 rows
post_thetas <- map_df(1:n_coders, function(j) 
  tibble(
    coder = j
    , theta0 = unlist(baseline_fit[, sprintf("theta0[%s]", j)] )
    , theta1 = unlist(baseline_fit[, sprintf("theta1[%s]", j)] )
  )
)

# post_thetas %>% 
#   gather(parameter, value, -coder) %>% 
#   # filter(parameter == "theta1") %>%
#   ggplot(aes(y = value, x = factor(coder), group = coder)) +
#   geom_violin() +
#   facet_grid(rows = vars(parameter))

post_thetas %>% 
  gather(parameter, value, -coder) %>% 
  ggplot(
    aes(
      y = factor(coder)
      , x = value
      , group = coder
      # , fill = factor(parameter)
    )
  ) +
  geom_density_ridges(
    scale = .5
    , color = NA
  ) +
  geom_point(
    shape = 3
    , color = "red"
    , data = tibble(
      coder = rep(1:n_coders, 2)
      , parameter = rep(paste0("theta", 0:1), each = n_coders)
      , value = c(theta0, theta1)
    )
    , aes(y = factor(coder), x = value, group = coder)
  ) +
  scale_x_continuous(limits = c(0,1))+
  facet_grid(rows = vars(parameter)) + 
  labs(
    title = "Marginal posterior densities of ability parameters by coder"
    , subtitle = "Red crosses mark 'true' values"
    , x = ""
    , y = "Coder"
    , fill = "Parameter"
  ) + 
  coord_flip() +
  theme_bw()
# We can clearly see how ability parameter estimates are consistently flipped relative to the 'true' simulation parameter values
# just as in the case of prevalence

post_thetas %>% 
  gather(param, value, -coder) %>% 
  group_by(coder, param) %>% 
  summarise(
    mean = mean(value)
    , q5 = quantile(value, .05)
    , q95 = quantile(value, .95)
  ) %>% 
  arrange(param, coder) %>% 
  ungroup() %>% 
  mutate(true_val = c(sim_theta0, sim_theta1)) %>% 
  ggplot(aes(x = true_val, y = mean)) +
  # geom_smooth(method='lm', formula=y~x, se = FALSE) +
  geom_point(size = .5) +
  geom_linerange(aes(ymin = q5, ymax = q95), width = .1) +
  scale_x_continuous(limits = 0:1) +
  scale_y_continuous(limits = 0:1) +
  geom_abline(slope = 1, intercept = 0, size = .5, color = "grey") +
  facet_grid(~param) +
  theme_bw() + 
  labs(
    title = "Mean posterior coder abilities vs. 'true' parameter values."
    , subtitle = "Vertical lines overlaying mean values depict 90% credibility intervall, the diagonal indicates perfect correspondence."
    , x = "Simulated parameter value"
    , y = "Posterior value"
  )

# The data pulls ability parameters towards low values, but the informative beta priors pull in the opposite direction
# we can also see how the more informative (more dispersed) prior on $\theta_0$ ($f_{B}(40,8)$ results in higher dispersion 
# of posteriors 


# Item classes 

# gives n_items*n_chains*n_iterations = 1000*3*1000 = 3,000,000 rows
post_c <- map_df(1:n_items, function(i) 
  tibble(item = i, c = unlist(baseline_fit[, sprintf("c[%s]", i)]))
)


post_c %>% 
  group_by(item) %>% 
  summarise(pr_pos = sum(c)/n()) %>% 
  mutate(est_class = as.integer(pr_pos > .5)) %>% 
  left_join(tibble(item = 1:n_items, true_class = sim_item_classes)) %>% 
  mutate(aggreement = true_class == est_class) %>% 
  summarise(accuracy = sum(aggreement)/n())
# Finally, we can see that the model-based classification of items is almost perfectly contradicting true class memberships
# In sum, evaluation of the baseline model suggests that we would want to constraint the model to induce convergence on true values and not mirroring ones
