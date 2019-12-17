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

#' Calculates correlation between two (symmetrical) matrices using `stats::cor`
#'
#' @param m1,m2 The two matrices. Function only keeps row and colnames that exist in m1 and m2.
#' @param zeros_remove Should all zeros be removed? (set to NA before running `cor`).
#'
#' @inheritParams stats::cor
#'

calc_matrix_correlation <- function(m1, m2, zeros_remove = FALSE,
                                    method = c("pearson", "kendall", "spearman")){
  method <- match.arg(method)
  all_names <- list(sort(union(dimnames(m1)[[1]], dimnames(m2)[[1]])))
  # Create empty matrices than includes all Individuals from m1 and m2
  m1_new <- m2_new <- matrix(nrow = length(all_names[[1]]),
                             ncol = length(all_names[[1]]),
                             dimnames = c(all_names, all_names))

  # Fill these matrcies with values from m1 and m2
  m1_new[rownames(m1_new) %in% rownames(m1),
         colnames(m1_new) %in% colnames(m1)] <- m1[rownames(m1) %in% rownames(m1_new),
                                                   colnames(m1) %in% colnames(m1_new)]

  m2_new[rownames(m2_new) %in% rownames(m2),
         colnames(m2_new) %in% colnames(m2)] <- m2[rownames(m2) %in% rownames(m2_new),
                                                   colnames(m2) %in% colnames(m2_new)]
  # Remove all zeros if argument is set to TRUE
  if(zeros_remove){
    m1_new[m1_new == 0] <- NA
    m2_new[m2_new == 0] <- NA
  }
  # Calculate mantel coefficient (as done in vegan::mantel)
  statistic <- cor(as.dist(m1_new), as.dist(m2_new),
                   use = "pairwise.complete.obs",
                   method = method)
  return(list(statistic = statistic,
              method = method))
}
