#' Uses a named list (or vector of values) to calculate the quantiles of
#' specified probabilities
#'
#' If `value_list` is a vector, the function creates a data frame (tibble) with
#' two columns ("probability" and "quantile_value") and one row. If `value_list`
#' is a named list, the function creates a data frame with three columns
#' ("group_var", "probability", and "quantile_value"), and one row per
#' list-element
#'

#' @param value_list numeric vector or named list to calculate quantiles.
#' @param probs As specified in `quantile`
#'
#' @export
#'
#'
list_to_quantiles <- function(value_list, probs = c(0.01, 0.98)){
  require(tidyverse)
  if(is.vector(value_list) & !is.list(value_list)) {
    q <- quantile(value_list, probs = probs) %>%
      enframe(name = "probability", value = "quantile_value") %>%
      mutate(probability = as.numeric(str_remove(probability, "\\%")),
             probability = probability/100)
  } else if(is.vector(value_list) & is.list(value_list)) {
    q <- lapply(value_list, function(x) quantile(x, probs = probs)) %>%
      unlist %>%
      tibble(var_name = names(.), quantile_value = .) %>%
      separate(col = var_name, into = c("group_var", "probability"), sep = "\\.",
               extra = "merge") %>%
      mutate(probability = as.numeric(str_remove(probability, "\\%")),
             probability = probability/100)
  } else{
    stop("value_list not a list a list of vector")
  }
  return(q)
}
