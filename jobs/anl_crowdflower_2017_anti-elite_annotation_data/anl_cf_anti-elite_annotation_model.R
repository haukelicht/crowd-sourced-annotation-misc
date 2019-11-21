# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title: Analyse crowd-sourced dataset on anti-elite messages in political
#   communication on Twitter and Facebook
# Author: Hauke Licht
# Date: 2018-07-12
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #

library(dplyr)
library(tidyr)
library(purrr)
library(ggplot2)
library(viridisLite)
theme_set(theme_minimal())

set.seed(12345)

plot_dir <- "./jobs/anl_crowdflower_2017_anti-elite_annotation_data/plots"
save_plot <- partial(ggsave, path = plot_dir, device = "png", width = 7, height = 7)



# load and prer data ----
foo <- data.table::fread("~/switchdrive/Documents/work/phd/projects/populism/crowd_sourcing_populism/codings/crowdflower-data-test1.csv")
dat <- data.table::fread("./exdata/CF-task1-codings-clean.csv", stringsAsFactors = FALSE)

str(foo)

str(dat)

annotations <- dat %>% 
  filter(filter == "ok") %>% 
  mutate_at(vars(1,3,5), Vectorize(function(x) if (x == "") NA_character_ else x)) %>% 
  mutate(
    index = row_number(),
    item = group_indices(., id) ,
    annotator = group_indices(., X_worker_id),
    elites = elites == "yes",
    exclusionary = exclusionary == "yes",
    people = people == "yes",
    populist = people & elites,
    right_populist = people & elites & exclusionary
  ) %>% 
  tbl_df()

str(annotations)
# source Bob Carpenter's implementation of the David-Skene EM-alogorithm for annotation modeling
source("./code/R/compute_maxexpec.R")

# Anti-elitism ----

# fit annotation model to anti-elite annotations 
elite_annotations <- annotations %>%
  rename(label = elites) %>% 
  filter(!is.na(label)) %>% 
  select(index, item, annotator, label) 

elites_fit <- compute_maxexpec(elite_annotations, deflt.accur = .6,verbose = F)

# majority-vote based estimate of a prevalence (18.92%) closely aligns with
# the parametric prevalence estimate of anti-elite messages (18.39%) obtained
# from the annotation model

elites_fit$est_class_prevl

elites_fit$est_class_probs %>% 
  select(-majority_vote) %>% 
  spread(category, est_prob, sep = "_") %>% 
  mutate(is_antielite = category_TRUE >= category_FALSE) %>% 
  group_by() %>% 
  summarize(prevalence = sum(is_antielite)/n())
# assigning messages to classes based on estimated probabilities of class membership,
# however, yields a somewhat higher 'observed' prevalence: 21.7% messages are 
# classified as anti-elite based on the model, in contrast to the 18.92% of messages
# a majority of annotators classified as anti-elite


# sample from items (i.e., social media posts)

# NOTE: Because the estimated parameters of the annotation model strongly depend on 
#   the data (only minimal smoohting is applied to account for missing annotations),
#   drawing multiple samples of messages and averaging estimates over these 
#   'bootstrapped' samples allows to examine how estimated prevalence, model-based
#   class assignments, and estimates of annotator qualities (i.e., per-annotator 
#   sensitivity and specificity) vary

# create empty list
boot_samples <- vector(mode = "list", length = 500L)
item_idxs <- unique(elite_annotations$item)
ssize <- ceiling(length(item_idxs)*.9)

# obtain EM estimates on 500 bootstrap samples
for (i in seq_along(boot_samples)) {
  if (i == 1 || i %% 50 == 0)
    cat(i, "... ")

  # sample fraction of message IDs (w/o replacement) 
  bs_idx <- sample(x = item_idxs, size = ssize, replace = FALSE)
  # cosntruct bootstrap sample: all annotations of sampled message 

  
  x <- elite_annotations[(elite_annotations$item %in% bs_idx),]
  # NOTE: The rationale to sample all annotations of a random subset
  #   of messages is that parameter estimates will not change due to 
  #   changing annotations, but due to estimation on a different 'corpus'
  
  # fit model
  est <- compute_maxexpec(x, verbose = FALSE)
  
  est[["bootstrapped_datasets"]] <- bs_idx
  
  boot_samples[[i]] <- est
  
  if (i == max(length(boot_samples)))
    cat("finished!")
}

# examine majority-vote vs. annotation-model based prevalence estimates
ae_sampl_prevalences <- map_df(boot_samples, "est_class_prevl", .id = "boot_sample") %>% 
  gather("meth", "prob", est_prob:prop_maj_winner, -category) %>% 
  mutate(meth = case_when(
    grepl("^est_", meth) ~ "model-based classification (MC)", 
    grepl("winner$", meth) ~ "majority voting (MV)",
    TRUE ~ "proportion of messages")
  ) %>% 
  filter(category == 1, !grepl("^prop", meth)) %>% 
  select(boot_sample, category, meth, prob) %>% 
  arrange(boot_sample, category, meth) %>% 
  tbl_df()

t_diff <- t.test(data = ae_sampl_prevalences, prob ~ meth, alternative= "less")

pp_voted_v_est_prevalence <- ggplot(
  data = ae_sampl_prevalences,
  aes(x = prob, fill = meth)
) + 
  geom_density() +
  scale_fill_viridis_d(alpha = .5) + 
  labs(
    title = "Estimated prevalence of anti-elite statements",
    subtitle = sprintf("Based on 500 bootstrap samples of %s statements each", ssize),
    x = "Est. prevalence",
    y = "",
    fill = "Decision rule"
  ) +
  geom_vline(xintercept = t_diff$estimate) + 
  annotate(geom = "text"
    , x = t_diff$estimate[1] - .001 , y = 2
    , label = paste("hat(italic(mu))[MV] ==", round(t_diff$estimate[1], 3))
    , parse = TRUE
    , hjust = 1
  ) +
  annotate(geom = "text"
    , x = t_diff$estimate[2] + .001 , y = 2
    , label = paste("hat(italic(mu))[MC] ==", round(t_diff$estimate[2], 3))
    , parse = TRUE
    , hjust = 0
  ) +
  labs(
    caption=bquote(
      "One-tailed"~
        italic(t)*"-test"~of~H[0]*": "~
        hat(italic(mu))[MV] >= hat(italic(mu))[MC]~
        "gives"~italic(t) %~~% .(round(t_diff$statistic, 3))*
        ", "~
        italic(p) %~~% .(round(t_diff$p.value, 5))
    )
  ) 
  

# distributions of bootstrapped point estimates tightly overlap
pp_voted_v_est_prevalence
save_plot(filename = "voted_v_est_prevalence", plot = pp_voted_v_est_prevalence) 

# identical mean, only one percentage point difference in medians
ae_sampl_prevalences %>% 
  group_by(meth) %>%
  summarize(
    mean = mean(prob),
    median = quantile(prob, .5),
    p025 = quantile(prob, .025),
    p975 = quantile(prob, .975)
  )

all.equal(
  boot_samples[[1]]$majority_votes
  , boot_samples[[2]]$majority_votes
)
# note: difference in majority-based classification induced by tie-breaking in 
# even-numbered sets of codings 


# examine uncertainty in model-based categorization of statements
ae_classified <- map_df(boot_samples, "est_class_probs", .id = "boot_sample") %>% 
  filter(category) %>% 
  group_by(item) %>% 
  mutate(
    est_prob_sd = sd(est_prob)
    , prop_mv = sum(majority_vote)/n()
  ) %>% 
  arrange(item)

# ambiguou instances uncovered by bootstrapping majority-winners 
hist(ae_classified$prop_mv)

ae_classified %>% 
  group_by(prop_mv) %>% 
  summarise(n = n_distinct(item))




pp_estimated_class_membership_vs_majority_vote <- ggplot(
  annotations %>% 
    group_by(item)  %>% 
    summarize(
      n = n()
      , n_elites = sum(elites)
    ) %>% 
    left_join(ae_classified) %>% 
    filter(est_prob_sd > .01) %>%
    mutate(
      voted = n_elites/n
      , voted = case_when(
          voted == 1 ~ "pure positive",
          voted > .5 ~ "impure positive",
          voted == .5 ~ "tie",
          voted == 0 ~ "pure negative",
          voted < .5 ~ "impure negative",
          TRUE ~ NA_character_
      )
    )
  ,
  aes(x = reorder(factor(item), est_prob), y = est_prob, 
      fill = voted
      # fill = majority_vote
      )
) + 
  geom_violin(
    scale = "width"
    # , trim = FALSE
    , adjust = .75
    , colour = NA#"grey"
    , width = .75
  ) + 
  # geom_boxplot(
  #   notch = FALSE,
  #   # notch = TRUE,
  #   # notchwidth = 1,
  #   outlier.shape = NA,
  #   show.legend = FALSE
  #   # show.legend = TRUE
  #   , fill = NA
  # ) +
  stat_summary(
    aes(group = factor(item))
    , fun.y = mean
    , geom = "point"
    , fill = "black"
    , shape = 1, size = 1
    , position = position_dodge(width = .7)
  ) + 
  stat_summary(
    fun.data = function(x, n = 1) {
      tibble(
        y = mean(x),
        sd = sd(x),
        ymin = y - n*sd,
        ymax = y + n*sd)
    }
    , geom = "errorbar"
    , colour = "black"
    , width = .25
    , size = .25
    , position = position_dodge(width = .7)
  ) + 
  scale_fill_viridis_d(alpha = .75) +
  coord_flip() +
  labs(
    title = "Anti-elite statements: Uncertainty in model- vs. voting based classification",
    # subtitle = sprintf("Uncertainty in model-based classification computed based on 500 bootstrap samples of %s statements each", ssize),
    y = "Est. probability of being an anti-elite statement (model-based)",
    x = "", # "Message ID",
    fill = "voting-based classification",
    caption = paste(
      ""
      , "Model-based classifications obtained by applying Dawid-Skene expectation-maximation algorithm to each of 500 resampled judgement sets."
      , "Uncertainty in model-based classification displayed as the variability in estimated class-membership probabities across reasampled jugement sets"
      , "Majority-vote classifications are generated by aggregating judgements at the level of text."
      , "Statments with a standard deviation in est. class probabilities > .01 not displayed."
      , sep = "\n"
    )
  ) + 
  theme(
    legend.position = "top"
    , plot.caption = element_text(hjust = 0)
  )

pp_estimated_class_membership_vs_majority_vote
save_plot(
  filename = "estimated_class_membership_vs_majority_vote.png", 
  plot = pp_estimated_class_membership_vs_majority_vote) 

# examine annotator performances ----
anno_perf <- map_df(boot_samples, "est_annotator_params", .id = "boot_sample") %>% 
  mutate(
    # boot_sample = rep(1:(nrow(.)/12), each = 12),
    metric = case_when(
      category & labeled ~ "true-positive rate (sensitivity)",
      !category & labeled ~ "false-positive rate (1 - sensitivity)",
      category & !labeled ~ "false-negative rate (1 - specificity)",
      !category & !labeled ~ "true-negative rate (specificity)"
    )
  ) %>% 
  select(boot_sample, annotator, metric, est_prob)

pp_annotator_biases <- ggplot(
  data = anno_perf %>% 
    filter(grepl("^true", metric)) %>% 
    arrange(desc(metric))
  , aes(x = factor(annotator), 
        y = est_prob,
        color = factor(annotator), 
        fill = factor(annotator))
) +
  # geom_boxplot(
  #   notch = TRUE,
  #   notchwidth = 1,
  #   outlier.shape = NA,
  #   show.legend = FALSE
  # ) +
  geom_violin(scale = "width") +
  scale_fill_viridis_d(alpha = .5) + 
  scale_color_viridis_d(alpha = .7) + 
  facet_grid(cols = vars(metric)) + 
  labs(
    title = "Estimated true-positive and true-negative rates of annotators",
    subtitle = sprintf("in 500 bootstrap samples of %s messages each", ssize),
    y = "Est. probability of agreement with model-based classification",
    x = "Annotator ID")

pp_annotator_biases
save_plot(filename = "annotator_biases", plot = pp_annotator_biases)

# We can observe different biases: 
#   - annotator 1 is strongly biased against detecting anti-elite messages 
#     (virtually always classifies negative instances correctly, but classifies
#      only about 4 out of 10 positive instances correctly)
#   - annotator 2, in contrast, is biased (though less severly) towards detecting 
#     anti-elite messages (detects about 8 out of 10 negative instances correctly,
#      but virtually always classifies positive instances correctly)
#   - annotator 3, in turn, is relatively unbiased, and correctly classifies 
#     about 9 out of 10 positive as well as negative instances
# Note: Positive/negative instances refers to model-based classifications of messages
#   as anti-elite and not anti-elite, respectively.