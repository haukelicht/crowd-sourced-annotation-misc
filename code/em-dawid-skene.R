# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title: 
# Author: Hauke Licht
# Date: 2018-07-11
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #


# 
library(dplyr)
library(tidyr)


dat <- read.csv("./exdata/CF-task1-codings-clean.csv", stringsAsFactors = FALSE)


str(dat)

head (dat,20)
dat %>% 
  group_by(elites) %>% 
  summarize()

table(dat$elites == "")

as.integer(dat$elites)

annotations <- dat %>% 
  filter(elites != "") %>%
  mutate(
    index = row_number(),
    item = group_indices(., id) ,
    annotator = group_indices(., X_worker_id),
    label = ifelse(elites == "yes", 1, 2)
  )  %>% 
  select(index, item, annotator, label)

table(annotations$item)

ii <- annotations$item
jj <- annotations$annotator
y <- annotations$label

N <- length(ii)
K <- length(unique(y))
J <- length(unique(jj))
I <- max(ii)


# ii <- data[[1]]
# jj <- data[[2]]
# y <- data[[3]]
# 
# N <- length(ii)
# K <- max(y)
# J <- max(jj)
# I <- max(ii)

# print data sizes
print(paste("K=",K, " N=",N, " J=",J, "I=",I))

##### EM ALGORITHM #####
# parameters
min.relative.diff = 1E-8
max.iters = 100

### INITIALIZATION

# create array with J K-by-K matrices defining annotators' accuracies
theta_hat <- array(NA, c(J, K, K))
# NOTE: theta_hat[j,k,k'] gives the probability that annotator j assigns the 
#       label k' to an item whose true category is k

# set annotator accuracies for first iteration step
for (j in 1:J) {
  # assuming that each annotator has .7 precision on each category, ...
  # ... set off-diagonal (error-rates)
  theta_hat[j,,] <- 0.3/(K-1)
  # ... and precision
  diag(theta_hat[j,,]) <- 0.7
}

theta_hat[1,,]

# set initial class prevalences/probabilites 
pi_hat <- array(1/K, K)

### EM ITERATIONS
iter_step <- 1
last_log_posterior = - Inf

# create items-by-categories array to store estimated class probabilities of items 
E_z <- array(1/K, c(I, K))
# NOTE: E_z[i,k] gives the probability that item i's true category is k

for (iter_step in 1:max.iters) {
  
# --- Expectation Step --- #
  # set estimated class probabilities to current estimates of category prevalences
  for (i in 1:I)
    E_z[i,] <- pi_hat
  
  # for each annotation ...
  for (n in 1:N) {
    # update estimated category probabilities by taking into 
    # account annotators' accuracies/error-rates
    E_z[ii[n], ] <- E_z[ii[n], ] * theta_hat[jj[n], , y[n]]
  }
  
  # rescale so that class probabilities sum to one
  E_z <- E_z / rowSums(E_z)
    
# --- Maximization Step --- #
  beta <- 0.01 
  # add beta smoothing on pi_hat
  pi_hat <- rep(beta, K)
  for (i in 1:I)
    pi_hat <- pi_hat + E_z[i,]
  
  # ensure that probabilities sum to one
  pi_hat <- pi_hat / sum(pi_hat)

  # add alpha smoothing for theta_hat
  alpha <- 0.01
  count <- array(alpha, c(J, K, K))
  # for each annotation ...
  for (n in 1:N) {
    # ... add alpha prior to estimated class probabilities
    count[jj[n], , y[n]] <- count[jj[n], , y[n]] + E_z[ii[n], ]
  }
  
  # for each annotator ...
  for (j in 1:J) {
    # ... and for each category ...
    for (k in 1:K)
      # ensure that individual error-rates sum to one 
      theta_hat[j,k,] <- count[j,k,] / sum(count[j,k,])
  }

  # create array 
  p <- array(0, c(I, K))
  
  # ... and assign smoothed estimated class prevalences
  for (i in 1:I)
    p[i,] <- pi_hat
  
  # for each annotation ...
  for (n in 1:N) {
    # ... and for each category ...
    for (k in 1:K) {
      # ... update 
      p[ii[n], k] <- p[ii[n], k] * theta_hat[jj[n], k, y[n]]
    }
  }

  # compute log of posterior probability of the data given current parameter estimates  
  log_posterior <- sum(log(rowSums(p)))

  # progress or terminate
  if (iter_step == 1) {
    print(
      paste("iteration nr. ", iter_step,
            ": log posterior = ", log_posterior, 
            sep = "")
    )
  } else {
    diff <- log_posterior - last_log_posterior
    relative_diff <- abs(diff / last_log_posterior)
    print(
      paste("iteration nr. ", iter_step,
            ": log posterior = ", log_posterior, 
            ", relative difference = ", relative_diff,
            sep = "")
    )
    if (relative_diff < min.relative.diff) {
      print("Converged!")
      break
    }
  }
  last_log_posterior <- log_posterior
}


# VOTED PREVALENCE AS A SANITY CHECK compare to estimates of pi ----
voted_prevalence <- rep(0,K)
for (k in 1:K)
  voted_prevalence[k] <- sum(y == k)
(voted_prevalence <- voted_prevalence / sum(voted_prevalence))

pi_out <- array(0, c(K, 2))
dimnames(pi_out) <- list(NULL, c("category", "prob"))
pos <- 1
for (k in 1:K) {
  pi_out[pos,] <- c(k,pi_hat[k])
  pos <- pos + 1
}
pi_out

theta_out <- array(0, c(J*K*K, 4))
dimnames(theta_out) <- list(NULL,c("annotator", "reference", "response", "prob"))

pos <- 1
for (j in 1:J) {
  for (ref in 1:K) {
    for (resp in 1:K) {
      theta_out[pos,] <- c(j, ref, resp, theta_hat[j, ref , resp])
      pos <- pos + 1
    }
  }
}

foo <- theta_out %>% 
  as.data.frame() %>% 
  mutate(prob = round(prob,5)) %>% 
  split(., .$annotator, drop = T)

df <- foo[[1]]



get_conf_matrix <- function(df) {
  out <- matrix(df$prob, length(unique(df$reference)), byrow = TRUE)
  rownames(out) <- c("category(true)", "category(false)")
  colnames(out) <- c("labeled(true)", "labeled(false)")
  out
}

lapply(foo, get_conf_matrix)
# hli: Fitting the annotation model reveals that the three annotators vary in 
#  quality: 
#   - Annotator 1 has the highest specificity (correctly classified 99.99% 
#      of non-anti-elitist instances), but suffers from low precision/sensitivity
#      (classified only 39.07% of anti-elitist instances correctly); this is 
#      evidence of a bias inducing her to miss actual anti-elitist messages.
#   - Annotator 2 has both high sensitivity (99.92%) and specificity (91.1%), and 
#      is thus only mildely biased towards detecting anti-elitism.
#   - Annotator 3 is, in contrast to Annotator 1, biased towards detecting 
#      anti-elitism in actualy non-anti-elitist messages (specificity = 77.8%, 
#      sensitivity = 91.54%)
# 




z_out <- array(0, c(I*K, 3))
dimnames(z_out) <- list(NULL, c("item", "category", "prob"))


library(tidyr)

z_out <- E_z %>% 
  as.data.frame() %>% 
  mutate(item = row_number()) %>% 
  gather(., category, prob, V1:V2) %>% 
  arrange(item, category) %>% 
  mutate(
    category = as.integer(sub("V", "", category)),
    prob = round(prob, 5)
  )


foo <- z_out %>% 
  group_by(item) %>% 
  filter(max(prob) < 0.90) %>% 
  ungroup()


library(ggplot2)

ggplot() + 
  aes(x = reorder(item, -prob), y = prob) +
  geom_bar(data = foo %>% filter(category == 1), width = .5, stat = "identity") + 
  geom_bar(data = foo %>% 
             filter(category == 2) %>% 
             mutate(prob = -prob), 
           width = .5, 
           stat = "identity") + 
  geom_hline(yintercept=0) +
  scale_y_continuous(name = "Probability of 'anti-elite' message content",
                     labels = paste0(seq(-100, 100, 20), "%"),
                     limits = c(-1, 1),
                     breaks = seq(-1, 1, .2)) +
  scale_x_discrete(name = "Message ID") +
  coord_flip()




