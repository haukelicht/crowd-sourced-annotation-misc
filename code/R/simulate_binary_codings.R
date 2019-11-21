#' Simulate binary codings data
#' 
#' @description Function simulates codings for a binary classification task
#' 
#' @param n.items scalar integer value, specifying the number of items to be coded 
#'    Defaults to 1000.
#'    
#' @param n.coders scalar integer value, specifying the number of coders sampled
#'    Defaults to 20.
#'    
#' @param missing.rate scalar numeric value \math{\in (0, 1)}, determining the number of judgments
#'     sampled per item. If \code{n.coders} is 10, and \code{missing.rate} is set to .5, for example,
#'     then \code{missing.rate} \math{\times} \code{n.coders} = 5 judgments from 5 different 
#'     coders per item will be sampled at random.
#'     Note that the \code{missing.rate} \math{\times} \code{n.coders} is rounded to the next bigger 
#'     integer value. 
#' 
#' @param pi scalar numeric value \math{\in (0, 1)}, specifying the true prevalence of positive items
#'    Defaults to 0.2.
#'    
#' @param alpha0 scalar numeric value, specifiying the first shape parameter governing the 
#'    hyperdistribution of coders' specificities (i.e., true-negative detection abilities)
#'    Defaults to 40.
#'    
#' @param beta0 scalar numeric value, specifiying the second shape parameter governing the 
#'    hyperdistribution of coders' specificities (i.e., true-negative detection abilities)
#'    Defaults to 8.
#'    
#' @param alpha1 scalar numeric value, specifiying the first shape parameter governing the 
#'    hyperdistribution of coders' sensitivities (i.e., true-positive detection abilities)
#'    Defaults to 20.
#'    
#' @param beta1 scalar numeric value, specifiying the second shape parameter governing the 
#'     hyperdistribution of coders' sensitivities (i.e., true-positive detection abilities)
#'     Defaults to 8.
#'    
#' @note Coders abilities are drawn from two Beta-distributions, governed by 
#'      \code{alpha0} and \code{beta0}, and \code{alpha1} and \code{beta1}, respectively
#'
#' @return A list with two elements: 
#'     (1) a \code[tibble]{tibble} object with columns 
#'        'item' (integer item indices), 
#'        'coder' (integer coder indices), 
#'        'judgment' (simulated judgments \math{\in \{0, 1\}}), and
#'        'true_class' (items' true classes \math{\in \{0, 1\}}); and
#'     (2) a list with two elements:
#'        (i) 'theta0' (a numeric vector of length \code{n.coders} with coders' simulated specificities), and
#'        (ii) 'theta1' (a numeric vector of length \code{n.coders} with coders' simulated sensitivities).
#'        
#' @importFrom purrr rbernoulli
#' @importFrom purrr map_df
#' @importFrom purrr map_int
#' @importFrom tibble tibble
#' 
#' @export
#' 
#' @example 
#' \donotrun{
#'   sim <- simulate_binary_codings(n.items = 10, n.coders = 3)
#' }
simulate_binary_codings <- function(
  n.items = 1000
  , n.coders = 20
  , missing.rate = .5
  , pi = .2
  , alpha0 = 40
  , beta0 = 8
  , alpha1 = 20
  , beta1 = 8
){
  
  # sample items' classes
  sim_item_classes <- as.integer(purrr::rbernoulli(n = n.items, p = pi))
  
  # sample 20 coders
  sim_theta0 <- rbeta(n = n.coders, shape1 = alpha0, shape2 = beta0)
  sim_theta1 <- rbeta(n = n.coders, shape1 = alpha1, shape2 = beta1)
  
  # for each item i
  sim_codings <- purrr::map_df(1:n.items, function(i){
    
    # sample coders who judge the item's class
    idxs <- sample(1:n.coders, ceiling(n.coders*(1-missing.rate)))
    
    # for each coder in the sample, 
    y_i <- purrr::map_int(idxs, function(idx) {
      # generate coder j's judgment given item i' s 'true' class:
      #  - if the 'true' class is positive, j classifies it as positive with probability theta_{j1}
      #  - if the 'true' class is negative, j classifies it as negative with probability theta_{j0}
      purrr::rbernoulli(1, p = sim_item_classes[i]*sim_theta1[idx] + (1 - sim_item_classes[i])*(1 - sim_theta0[idx]) )
    })
    
    # gather data and return 
    tibble::tibble(
      # the item index (ID)
      item = i
      # coders indices (IDs)
      , coder = idxs
      # simulated judgments 
      , judgment = y_i
      # item i's 'true' class
      , true_class = sim_item_classes[i] 
    )
  })
  
  list(
    codings = sim_codings
    , abilities = list(theta0 = sim_theta0, theta1 = sim_theta1)
  )
}
