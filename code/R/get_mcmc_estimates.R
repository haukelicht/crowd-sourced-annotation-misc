
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
