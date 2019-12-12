#' Function to calculate the simple ratio association index for several time period.
#'
#' So far, only works with Red Colobus Nearest Neighbour data
#'
#'
#' @param dyadic_obs_list List created with `get_dyadic_observations`
#' @param i_j_observed_sum_min Minimum number of times individual i and j were observed while both adult and present
#' @param new_dim_names A vector with two elements for column and row names (default is Name1, Name2)
#'
#' @export
#'
#' @examples
#'
#'

get_simple_ratios <- function(dyadic_obs_list, i_j_observed_sum_min, new_dim_names = c("Name1", "Name2")){
  # Create empty list for simple ratio index and copy attributes from source list
  simple_ratio_list <- vector("list", length = length(dyadic_obs_list))
  attributes(simple_ratio_list) <- attributes(dyadic_obs_list)
  attr(simple_ratio_list, "i_j_observed_sum_min") <- i_j_observed_sum_min

  # Calculate simple ratio for observed and all permuted data
  for(i in 1:length(dyadic_obs_list)){
    cat(names(dyadic_obs_list)[[i]], "\n")

    ### 1. i_j_observed_sum (required to calculate the ratios below)
    # Get the number of times an individual was observed while the other
    # individual was in the group, Then make this matrix symmetrical: this is the
    # sum of the number of time both individuals were recorded regardless of
    # whether they were together (nearest neighbours). But this is the max. number
    # of times they could have been observed together.
    # With regard to association indices, this is d, which is a function of
    # x (i and j associated), yij (i and j observed but not associated) and
    # yi/yj (number of times only individual i/j observed). See, e.g. Whitehead 2008
    i_j_observed_sum <- dyadic_obs_list[[i]]$n_Ind_while_NN_present
    i_j_observed_sum <- i_j_observed_sum + t(i_j_observed_sum)

    # Set all values below minimum number of observations to NA
    i_j_observed_sum[i_j_observed_sum < i_j_observed_sum_min] <- NA

    ### 2. simple_ratio_observed
    i_j_together_observed <- dyadic_obs_list[[i]]$n_Ind_NN_together_observed
    # Replace NAs by 0, and summarize both sides of matrixs
    i_j_together_observed[is.na(i_j_together_observed)] <- 0
    i_j_together_observed <- i_j_together_observed + t(i_j_together_observed)
    # Calculate simple ratio index
    simple_ratio_list[[i]]$simple_ratio_observed <- i_j_together_observed/i_j_observed_sum

    ### 3. simple_ratio_permuted
    i_j_together_permuted <- dyadic_obs_list[[i]]$n_Ind_NN_together_permuted
    # Replace NAs by 0, and summarize both sides of matrixs
    i_j_together_permuted[is.na(i_j_together_permuted)] <- 0
    for(j in 1:dim(i_j_together_permuted)[[3]]){
      i_j_together_permuted[,,j] <- i_j_together_permuted[,,j] + t(i_j_together_permuted[,,j])
    }
    # Create empty arry for the simple ratios
    simple_ratio_list[[i]]$simple_ratio_permuted <- i_j_together_permuted
    simple_ratio_list[[i]]$simple_ratio_permuted[,,] <- NA

    #  Calculate simple ratio index for all permutations and fill the array
    for(j in 1:dim(i_j_together_permuted)[[3]]){
      simple_ratio_list[[i]]$simple_ratio_permuted[,,j] <- i_j_together_permuted[,,j]/i_j_observed_sum
    }
    # Rename the Individual dimensions
    names(dimnames(simple_ratio_list[[i]]$simple_ratio_observed)) <- new_dim_names
    names(dimnames(simple_ratio_list[[i]]$simple_ratio_permuted))[1:2] <- new_dim_names
  }
  return(simple_ratio_list)
}
