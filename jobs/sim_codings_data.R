# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #  
#
# Title:  Sumulate codings data
# Author: Hauke Licht
# File:   sim_codings_data.R
# Date:   2019-03-03
#
# +~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~+~ #

req_packages <- c("purrr", "tibble")

idx <- !c("purrr", "tibble") %in% loadedNamespaces()

for (pkg in req_packages[idx])
  require(as.character(pkg))


# sample items' classes
sim_item_classes <- as.integer(rbernoulli(n = n_items, p = pi))

# sample 20 coders
sim_theta0 <- rbeta(n = n_coders, shape1 = alpha0, shape2 = beta0)
sim_theta1 <- rbeta(n = n_coders, shape1 = alpha1, shape2 = beta0)

# for each item i
sim_codings <- purrr::map_df(1:n_items, function(i){
  
  # sample ~10 coders who judge the item's class
  idxs <- sample(1:n_coders, ceiling(n_coders*(1-missingness_rate)))
  
  # for each coder in the sample, 
  y_i <- map_int(idxs, function(idx) {
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
