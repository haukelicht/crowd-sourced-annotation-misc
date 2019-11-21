#' Function creating JAGS model data input from codings data
#' 
#' @description Given a codings dataframe with columngs 'item', 'coder', and 'judgement', 
#'     function returns a list object that can be passed to \code{data} argument 
#'     of \[rjags]{jags.model} function fitting a beta-binomial by annotator model
#'     
#' @param Given codings dataframe with columngs 'item', 'coder', and 'judgement'
#' 
#' @return A list object with elements 
#'     'N' (the number of items, scalar integer),
#'     'M' (the number of coders, scalar integer),
#'     'J' (the number of judgements, scalar integer), and
#'     'y' (dataframe with columns 'item', 'coder' and 'judgment').
get_codings_data <- function(codings.df) {
  
  out <- list() 
  
  out$N <- length(unique(codings.df$item))
  out$M <- length(unique(codings.df$coder))
  out$J <- length(codings.df$judgment)
  out$y <- codings.df %>% 
    arrange(item, coder) %>% 
    select(item, coder, judgment)
  
  return(out)
}
