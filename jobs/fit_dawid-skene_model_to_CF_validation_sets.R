# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title: Fit Dawid-Skene annotation model to crowd-sourced validation dataset on 
#         populist instances among political statements on Twitter and Facebook.
# Author: Hauke Licht
# Date: 2019-02-19
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #

# setup ----
library(dplyr)
library(tidyr)
library(purrr)

source("./code/R/compute_maxexpec.R")

set.seed(12345)


# load and prep data ----
dat <- data.table::fread("./exdata/CF-task1-codings-clean.csv", stringsAsFactors = FALSE)
table(dat$filter)

dat %>% 
  # filter(filter == "language") %>% 
  # filter(filter == "notpoli") %>% 
  # filter(filter == "spam") %>% 
  # filter(filter == "language\nok") %>% 
  filter(filter == "spam\nok") %>% 
  head(1)


annotations <- dat %>% 
  # keep only judgements that are 'ok' (i.e., not spam, non-political, or non-comprehended)
  filter(grepl("ok$", filter)) %>% 
  # replace "" with NA in string vectors
  mutate_if(is.character, Vectorize(function(x) if (x == "") NA_character_ else x)) %>%
  mutate(
    # judgement index
    index = row_number(),
    # item index
    item = group_indices(., id) ,
    # annotator index
    annotator = group_indices(., X_worker_id),
    # populism indiactors
    elites = elites == "yes",
    exclusionary = exclusionary == "yes",
    people = people == "yes",
    populist = people & elites,
    right_populist = people & elites & exclusionary
  ) %>% 
  tbl_df()

unique(annotations$X_worker_id)
length(unique(annotations$id))



# NOTE: no translations available (in orginal test, coders were asked to copy-paste 
# statements they could not understand into Google Translate)
str(annotations)

# fit models ----

# Anti-elitism 
elite_annotations <- annotations %>%
  rename(label = elites) %>% 
  filter(!is.na(label)) %>% 
  select(index, item, annotator, label) 

# fit annotation model to anti-elite annotations 
elites_fit <- compute_maxexpec(elite_annotations, deflt.accur = .6,verbose = F)


# People-centrism
people_annotations <- annotations %>%
  rename(label = people) %>% 
  filter(!is.na(label)) %>% 
  select(index, item, annotator, label) 

# fit annotation model to anti-elite annotations 
people_fit <- compute_maxexpec(people_annotations, deflt.accur = .6,verbose = F)


# Exclusionism
exclu_annotations <- annotations %>%
  rename(label = exclusionary) %>% 
  filter(!is.na(label)) %>% 
  select(index, item, annotator, label) 

# fit annotation model to anti-elite annotations 
exclusionary_fit <- compute_maxexpec(exclu_annotations, deflt.accur = .6,verbose = F)


# item-level bootstrap  ----

# NOTE: Because the estimated parameters of the annotation model strongly depend on 
#   the data (only minimal smoohting is applied to account for missing annotations),
#   drawing multiple samples of messages and averaging estimates over these 
#   'bootstrapped' samples allows to examine how estimated prevalence, model-based
#   class assignments, and estimates of annotator qualities (i.e., per-annotator 
#   sensitivity and specificity) vary

all.equal(
  unique(elite_annotations$item)
  , unique(people_annotations$item)
)

all.equal(
  unique(elite_annotations$item)
  , unique(exclu_annotations$item)
)

all.equal(
  unique(people_annotations$item)
  , unique(exclu_annotations$item)
)

all.equal(
  unique(annotations$item)
  , unique(exclu_annotations$item)
)

# get item indices
item_idxs <- unique(annotations$item)
# set number of items per bootsrap sample
ssize <- ceiling(length(item_idxs)*.9)
# create empty list
boot_samples <- vector(mode = "list", length = 500L)

# obtain EM estimates on 500 bootstrap samples
for (i in seq_along(ae_boot_samples)) {
  if (i == 1 || i %% 50 == 0)
    cat(i, "... ")
  
  # sample fraction of item IDs (w/o replacement) 
  bs_idx <- sample(x = item_idxs, size = ssize, replace = FALSE)
  
  # NOTE: The rationale to sample all annotations of a random subset
  #   of messages is that parameter estimates will not change due to 
  #   changing annotations, but due to estimation on a different 'corpus'
  
  # fit models
  est <- list(
    "sampled_items" = bs_idx
    , "anti_elitism" = elite_annotations %>% 
      filter(item %in% bs_idx) %>% 
      compute_maxexpec(verbose = FALSE)
    , "people_centrism" = people_annotations %>% 
      filter(item %in% bs_idx) %>% 
      compute_maxexpec(verbose = FALSE)
    , "exclusionary" = exclu_annotations %>% 
      filter(item %in% bs_idx) %>% 
      compute_maxexpec(verbose = FALSE)
  )
  
  boot_samples[[i]] <- est
  
  if (i == max(length(boot_samples)))
    cat("finished!")
}


# create output object -----


out <- list(
  data = annotations
  , anti_elitism = list(
    fit = elites_fit
    , bootstraps = map(boot_samples, "anti_elitism")
  )
  , people_centrism = list(
    fit = people_fit
    , bootstraps = map(boot_samples, "people_centrism")
  )
  , exclusionary = list(
    fit = exclusionary_fit
    , bootstraps = map(boot_samples, "exclusionary")
  )
)

saveRDS(out, "./exdata/em_fits_CF_validation_set.RData")
