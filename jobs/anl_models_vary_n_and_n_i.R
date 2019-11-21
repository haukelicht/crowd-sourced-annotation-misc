
people_deviances <- map_df(people_dat, function(g) {
  map_df(g, function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "deviance") %>% 
      mutate(
        n_i = unique(l$data$n_i)
        , sample_size = unique(l$data$sample_size)
      )
  })
})

# all chains mix nicely
ggplot(people_deviances, aes(x = iter, y = est, color = factor(chain))) +
  geom_line(alpha = .5) +
  facet_grid(rows = vars(sample_size), cols = vars(n_i), switch = "y", scales = "free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 180))

for(g in seq_along(people_dat)) {
  for (n in 1:8) {
    
    nm <- sprintf("shrinkage_w_%s_jugments_for_each_of_%s_items", n+2, seq(200, 500, 50)[g])
    
    png(file.path(file_path, "plots", paste0(nm, ".png")))
    gelman.plot(people_dat[[g]][[n]]$fit[,"deviance"], main = gsub("_", " ", nm))
    dev.off()
  } 
}

for(g in seq_along(people_dat)) {
  for (n in 1:8) {
    
    nm <- sprintf("autocorrelation_w_%s_jugments_for_each_of_%s_items", n+2, seq(200, 500, 50)[g])
    
    png(file.path(file_path, "plots", paste0(nm, ".png")))
    autocorr.plot(people_dat[[g]][[n]]$fit[,"deviance"], ask = FALSE, main = gsub("_", " ", nm), lag.max = 50)
    dev.off()
  } 
}


# inspect Pis

people_pis <- map_df(people_dat, function(g) {
  map_df(g, function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "pi") %>% 
      mutate(
        n_i = unique(l$data$n_i)
        , sample_size = unique(l$data$sample_size)
      )
  })
})

ggplot(people_pis, aes(x = iter, y = est, color = factor(chain))) +
  geom_line(alpha = .5) +
  facet_grid(rows = vars(sample_size), cols = vars(n_i), switch = "y", scales = "free") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 180))

ggplot(people_pis, aes(x = est)) +
  geom_density(color = NA, fill = "grey", alpha = .75) +
  facet_grid(cols = vars(sample_size), rows = vars(n_i), switch = "y") + 
  theme_bw() + 
  theme(strip.text.y = element_text(angle = 180))

# inspect c[.]s
people_cs <- map_df(people_dat, function(g) {
  map_df(g, function(l) {
    get_mcmc_estimates(fit.obj = l$fit, params = "c\\[\\d+\\]") %>% 
      mutate(
        n_i = unique(l$data$n_i)
        , sample_size = unique(l$data$sample_size)
      )
  })
})

people_cs_sum <- people_cs %>% 
  group_by(sample_size, n_i, parameter) %>% 
  summarise(pcu = sd(est)) %>% 
  group_by(sample_size, n_i) %>% 
  summarise(
    mean_pcu = mean(pcu)
    , sd_pcu = sd(pcu)
    , lwr_pcu = quantile(pcu, .05)
    , upr_pcu = quantile(pcu, .95)
  )

people_cs_sum %>% 
  ggplot() +
  geom_violin(
    data = people_cs %>% 
      group_by(sample_size, n_i, parameter) %>% 
      summarise(pcu = sd(est))
    , mapping = aes(x = n_i, y = pcu, group = n_i)
    , color = NA, fill = "grey", alpha = .76 
  ) +
  geom_point(aes(x = n_i, y = mean_pcu)) +
  # geom_line(aes(x = n_i, y = mean_pcu)) +
  geom_linerange(aes(x = n_i, ymin = lwr_pcu, ymax = upr_pcu)) +
  facet_grid(rows = vars(sample_size)) + 
  theme_bw()

people_cs_sum %>% 
  ggplot(
    aes(
      x = n_i
      , y = mean_pcu
      , ymin = lwr_pcu
      , ymax = upr_pcu
      , color = factor(sample_size)
      , group = sample_size
      , alpha = .75
    )
  ) +
  geom_line() +
  geom_point(size = 3.5, color = "white", alpha = 1) +
  geom_point(size = 3, shape = 1) +
  geom_point(size = 1) +
  scale_alpha(guide = FALSE) +
  theme_bw() +
  labs(
    title = "Change in average posterior classification uncertainty"
    # , subtitle = bquote("People-centrism classification simulated with"~
    #     pi==.(round(people_param_means$pi, 3))~","~
    #     theta[j0]~"~"~f[Beta](.(round(people_param_means$alpha0, 3)),.(round(people_param_means$beta0, 3)))~
    #     "and"~
    #     theta[j1]~"~"~f[Beta](.(round(people_param_means$alpha1, 3)),.(round(people_param_means$beta1, 3)))
    # )
    , x = expression(paste("No. judgement aggregated per item, ", n[i]))
    , y = expression(paste("Average posterior classification uncertainty, ", bar(SD)(c[i])))
    , color = "No. items"
  )

# inspect accuracy and bias
people_classification_quality <- sim_people_codings$codings %>% 
  select(item, true_class) %>% 
  unique() %>% 
  mutate(parameter = sprintf("c[%s]", item)) %>% 
  right_join(people_cs) %>% 
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
    , F1score = 2*((Precision*Recall)/(Precision + Recall))
  ) %>% 
  ungroup()

get_q05 <- function(x, na.rm) quantile(x, .05, na.rm = na.rm)
get_q95 <- function(x, na.rm) quantile(x, .95, na.rm = na.rm)

people_classification_quality_sum <- people_classification_quality %>% 
  group_by(sample_size, n_i) %>% 
  summarise_at(
    vars(Accuracy, ends_with("R", ignore.case = FALSE), Precision, Recall, F1score)
    , funs(mean, sd, q05 = get_q05, q95 = get_q95, .args = list(na.rm = TRUE))
  )  %>% 
  gather(metric, value, -sample_size, -n_i) %>% 
  separate(metric, c("metric", "statistic"), sep = "_") %>% 
  spread(statistic, value)

these_metrics <- c(
  "Accuracy"
  , "TNR"
  # , "TPR"
  , "Precision"
  , "Recall"
  , "F1score"
)

(p_people_sum_performance_distribution <- people_classification_quality_sum %>% 
  filter(metric %in% these_metrics) %>% 
  ggplot() +
  facet_grid(
    cols = vars(sample_size)
    , rows = vars(metric)
    , scales = "free_y"
    , switch = "y"
  ) +
  geom_point(aes(x = n_i, y = mean)) +
  geom_linerange(aes(x = n_i, ymin = q05, ymax = q95), size = .2) +
  scale_x_continuous(limits = c(2.7, 10.3)) +
  theme_bw() + 
  labs(
    title = "Distribution of performance metrics in simulated people-centrism classification"
    , y = NULL
    , x = expression(n[i])
    # , caption = "Panels separated by the number of items coded from left to right, and performance metrics from top to bottom"
  ) #+ 
  # theme(strip.text.y = element_text(angle = 90))
)


(p_people_sum_performance_change <- people_classification_quality_sum %>% 
  filter(metric %in% these_metrics) %>% 
  ggplot(
    aes(
      x = n_i
      , y = mean
      , color = factor(sample_size)
      , group = sample_size
      , alpha = .75
    )
  ) +
  facet_grid(cols = vars(metric), scales = "free_y") + 
  geom_line() +
  geom_point(size = 3.5, color = "white", alpha = 1) +
  geom_point(size = 3, shape = 1) +
  geom_point(size = 1) +
  scale_alpha(guide = FALSE) +
  theme_bw() +
  labs(
    title = "Change in performance of posterior classifications"
    , subtitle = bquote("People-centrism classification simulated with"~
                          pi==.(round(people_param_means$pi, 3))~","~
                          theta[j0]~"~"~f[Beta](.(round(people_param_means$alpha0, 3)),.(round(people_param_means$beta0, 3)))~
                          "and"~
                          theta[j1]~"~"~f[Beta](.(round(people_param_means$alpha1, 3)),.(round(people_param_means$beta1, 3)))
    )
    , x = expression(paste("No. judgement aggregated per item, ", n[i]))
    , y = NULL
    , color = "No. items"
  )
)
