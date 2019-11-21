#' @title
#'
#' @description
#'
#' @note Code adapted from Bob Carpenter's implementation of 
#'    of Dawid and Skene (1979) "Maximum Likelihood Estimation
#     of Observer Error-Rates Using the EM Algorithm"
#'    (see \url{https://github.com/bob-carpenter/anno/blob/master/R/em-dawid-skene.R})
#'
#' @param x a data frame object with columns 
#'     \itemize{
#'       \item 'index' (integer indexing judgement of item $i$ by annotator $j$),
#'       \item 'item' (integer indexing items $i = 1, ...., n$), 
#'       \item 'annotator' (integer indexing annotators $j = 1, ...., m$), and
#'       \item 'label' (value $y_{ij}$ assigned to item $i$ by annotator $j$, where $y \in \mathcal{K}$, the set of outcome categories)
#'     }
#'
#'
#'
#'
#'
#'
#' @return
compute_maxexpec <- function(
  x,
  item.col,
  annotator.col,
  annotation.col,
  deflt.accur = 0.7,
  beta.prior = 0.01,
  alpha.prior = 0.01,
  min.relative.diff = 1e-8,
  max.iters = 100,
  verbose = TRUE
){
  
  # Define internal helpers -----
  if (verbose)
    echo <- match.fun("cat")
  else
    echo <- function(...) {}
  
  # Test input ----
  stopifnot(is.data.frame(x))

  req_colns <- c("item", "annotator", "label")
  
  if (!all(req_colns %in% colnames(x))) {
    stop("input to argument `x' must be a data.frame object with columns:", paste("'", req_colns, "'", sep = ""))
  }
  
  # Setup data structure ----
  
  `%>%` <- magrittr::`%>%`
  
  x <- x %>% 
    dplyr::mutate(
      `_item` = dplyr::group_indices(., item) ,
      `_annotator` = dplyr::group_indices(., annotator),
      `_label` = dplyr::group_indices(., label)
    ) 

  # item vector and number of items 
  ii <- x$`_item`
  I <- length(unique(ii))
  
  # annotator vector and number of annotators
  jj <- x$`_annotator`
  J <- length(unique(jj))
  
  # label vector and number of categories (unique labels)
  y <- x$`_label`
  K <- length(unique(y))
  
  # number of annotations
  N <- length(ii)
  
  # test if complete panel design (i.e., all annotators labeled all items)
  complete_panel_design <- TRUE
  for (i in 1:I) {
    if (nrow(x[x$`_item` == i,]) != J)
      complete_panel_design <- FALSE
  }
  
  # report data sizes (if verbose)
  echo("Data input: \n")
  echo(" ", I, "items\n")
  echo(" ", J, "annotators\n")
  echo(" ", K, "categories\n")
  echo(sprintf("  annotations made in %scomplete panel design", ifelse(complete_panel_design, "a ", "an in")), "\n")
  echo("\n")
  
  ##### EM ALGORITHM #####
  
  # initialize parameter estimates
  
  # create array with J K-by-K matrices defining annotators' confusion matrices
  theta_hat <- array(NA, c(K, K, J))
  
  # set annotator accuracies for first iteration step
  for (j in 1:J) {
    theta_hat[,,j] <- (1-deflt.accur)/(K-1) # set off-diagonal (error-rates)
    diag(theta_hat[,,j]) <- deflt.accur # and diagonal
  }
  # NOTE: theta_hat[k,k',j] gives the probability that annotator j assigns the 
  #       label k' to an item whose true category is k
  
  # set initial class prevalences/probabilites 
  pi_hat <- array(1/K, K)

  # create items-by-categories array to store estimated class probabilities of items 
  E_z <- array(1/K, c(I, K))
  # NOTE: E_z[i,k] gives the probability that item i's true category is k
  
  # EM ITERATIONS ----
  
  # setup for first iteration step
  last_log_posterior = - Inf
  
  for (iter_step in 1:max.iters) {
    
    # --- Expectation Step --- #
    # set estimated class probabilities to current estimates of category prevalences
    for (i in 1:I)
      E_z[i,] <- pi_hat
    
    # for each annotation ...
    for (n in 1:N) {
      # update estimated category probabilities by taking into account annotator qualities
      E_z[ii[n], ] <- E_z[ii[n], ] * theta_hat[, y[n], jj[n]]
    }
    
    # rescale so that class probabilities sum to one
    E_z <- E_z / rowSums(E_z)
    
    # --- Maximization Step --- #
    # add beta smoothing on pi_hat
    pi_hat <- rep(beta.prior, K)
    for (i in 1:I)
      pi_hat <- pi_hat + E_z[i,]
    
    # ensure that probabilities sum to one
    pi_hat <- pi_hat / sum(pi_hat)
    
    # add alpha smoothing for theta_hat
    count <- array(alpha.prior, c(K, K, J))
    # for each annotation ...
    for (n in 1:N) {
      # ... add alpha prior to estimated class probabilities
      count[, y[n], jj[n]] <- count[, y[n], jj[n]] + E_z[ii[n], ]
    }
    
    # for each annotator ...
    for (j in 1:J) {
      # ... and for each category ...
      for (k in 1:K)
        # ensure that individual error-rates sum to one 
        theta_hat[k,,j] <- count[k,,j] / sum(count[k,,j])
    }
    
    # create array for individual class posterior probabilities 
    p <- matrix(rep(pi_hat, I), ncol = 2, byrow = TRUE)
    
    # for each annotation ...
    for (n in 1:N) {
      # ... and for each category ...
      for (k in 1:K) {
        # ... update 
        p[ii[n], k] <- p[ii[n], k] * theta_hat[k, y[n], jj[n]]
      }
    }
    
    # compute log of posterior probability of the data given current parameter estimates  
    log_posterior <- sum(log(rowSums(p)))
    
    # decide wether to progress
    if (iter_step == 1) {
      echo("Estimation results:\n")
      echo(
        paste("  iteration nr. ", iter_step,
              ": log posterior = ", round(log_posterior, 7), 
              "\n", sep = "")
      )
    } else {
      diff <- log_posterior - last_log_posterior
      relative_diff <- abs(diff / last_log_posterior)
      echo(
        paste("  iteration nr. ", iter_step,
              ": log posterior = ", round(log_posterior, 7), 
              ", relative difference = ", relative_diff,
              "\n", sep = "")
      )
      if (relative_diff < min.relative.diff) {
        echo("Converged!\n")
        break
      }
    }
    last_log_posterior <- log_posterior
  }
  
  # create data frames for function output ----
  
  # compute proportion of labels 
  prop_labels <- rep(0,K)
  for (k in 1:K)
    prop_labels[k] <- sum(y == k) / N
  
  # compute majority votes
  majority_votes <- x %>% 
    dplyr::group_by(`_item`, item, label) %>% 
    dplyr::summarize(votes = n_distinct(item)) %>% 
    dplyr::group_by(`_item`, item, label, votes) %>% 
    dplyr::mutate(
      tie_breaker = runif(1),
      majority_vote = label  
    ) %>% 
    dplyr::group_by(item) %>% 
    dplyr::filter(votes == max(votes)) %>% 
    dplyr::slice(which.max(tie_breaker)) %>% 
    dplyr::select(item, majority_vote) %>% 
    ungroup()

  z_out <- tibble::tibble(
    this_item = rep(1:I, each=K),
    cat_label = rep(1:K, times=I),
    est_prob = as.vector(t(E_z))
  ) %>% 
    dplyr::left_join(., unique(x[,c("item", "_item")]), by = c("this_item" = "_item")) %>%
    dplyr::left_join(., unique(x[,c("label", "_label")]), by = c("cat_label" = "_label")) %>% 
    dplyr::rename(category = label) %>% 
    dplyr::left_join(., majority_votes, by = "item") %>%
    dplyr::mutate(majority_vote = majority_vote == category) %>% 
    dplyr::select(item, category, est_prob, majority_vote)
  
  prop_maj_winner <- z_out %>% 
    dplyr::group_by(category) %>% 
    dplyr::summarize(prop_maj_winner = sum(majority_vote)/I)
  
  pi_out <- tibble::tibble(
    cat_label = sort(unique(y)), 
    est_prob = pi_hat, 
    prop_labels = prop_labels
  ) %>% 
    dplyr::left_join(., unique(x[,c("label", "_label")]), by = c("cat_label" = "_label")) %>% 
    dplyr::rename(category = label) %>% 
    dplyr::select(category, est_prob, prop_labels) %>% 
    dplyr::left_join(., prop_maj_winner, by = "category")
  
  theta_out <- tibble::tibble(
    this_annotator = rep(1:J, each=K*K),
    cat_label = rep(rep(1:K, times=K), times=J),
    an_labeled = rep(rep(1:K, each=K), times=J),
    est_prob = as.vector(theta_hat)
  ) %>% 
    dplyr::left_join(., unique(x[,c("annotator", "_annotator")]), by = c("this_annotator" = "_annotator")) %>%
    dplyr::left_join(., unique(x[,c("label", "_label")]), by = c("an_labeled" = "_label")) %>% 
    dplyr::rename(labeled = label) %>% 
    dplyr::left_join(., unique(x[,c("label", "_label")]), by = c("cat_label" = "_label")) %>% 
    dplyr::rename(category = label) %>% 
    dplyr::select(annotator, category, labeled, est_prob)
  
  # return output as list
  list(
    # estimated class probabilites
    est_class_probs = z_out,
    # majority votes
    majority_votes = majority_votes,
    # estimated class prevalence
    est_class_prevl = pi_out,
    # estimated annotator parameters (accuracies and erro rates)
    est_annotator_params = theta_out
  )
}

