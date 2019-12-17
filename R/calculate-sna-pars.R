#' Calculate various individual social network parameters based on DSI values
#'
#' Determines Strength, Betweenness, Eigenvector centrality, Clustering
#' coefficient (cc), Reach. Works with weighted, but undirected edges
#'
#' @param edge_table Table with two individual and one weight column
#' @param indA_col,indB_col,weight_col Specify names of columns for the two
#'   individuals and the weight to be used in the SNA
#'
#' @export
#'
#' @examples
#'

calc_sna_pars_dyad <- function(edge_table, indA_col = "IndA", indB_col = "IndB", weight_col = "DSI") {
  df <- edge_table[,c(indA_col, indB_col)]
  df$weight <- edge_table[,weight_col][[1]]

  # Check if none of the dyads is listes more than once
  if(any(
    (df %>%
     group_by(!!sym(indA_col), !!(sym(indB_col))) %>%
     summarize(n_dyad = n()))$n_dyad > 1)) {
    warning("at least one dyads has more than one value in the edge_table")
  }

  # igraph is masking too many other functions from the tidyverse, thus don't attach the package
  # create igraph object
  net <- igraph::graph_from_data_frame(d = df, vertices = NULL, directed = F)

  # remove edges with no weight (only keep positive edges) because these edges cause problems with some igraph-function (e.g. betweeness)
  net.pos <- net - igraph::E(net)[igraph::E(net)$weight == 0]

  net.weights <- igraph::E(net.pos)$weight
  net.vertices <- igraph::V(net.pos)

  vertex_strength <- igraph::strength(net.pos, weights = net.weights)
  # Use inverse of weights because in igraph, weights are considered as costs
  vertex_betweenness <- igraph::betweenness(net.pos, directed = FALSE, weights = 1 / net.weights)
  vertex_centrality <- igraph::eigen_centrality(net.pos, directed = FALSE, scale = FALSE, weights = net.weights)
  # Use "weigted" for transitivity because "local" assumes weight of 1 for all edges
  vertex_cc <- igraph::transitivity(net.pos, type = "weighted", vids = net.vertices, weights = net.weights)
  # For reach, use function from CePa (see below)
  # and use inverse because low reach means max(shortest.path) is small --> high reach
  vertex_reach <- 1/reach(graph = net.pos, weights = (1/net.weights))

  # Create dataframes for each parameter to make joining command easier to read
  vertex_strength_df <- data_frame(Ind = names(vertex_strength), strength = vertex_strength)
  vertex_betweenness_df <- data_frame(Ind = names(vertex_betweenness), betweenness = vertex_betweenness)
  vertex_centrality_df <- data_frame(Ind = names(vertex_centrality$vector), centrality = vertex_centrality$vector)
  vertex_cc_df <- data_frame(Ind = igraph::as_ids(net.vertices), cc = vertex_cc)
  vertex_reach_df <- data_frame(Ind = names(vertex_reach), reach = vertex_reach)

  # Then join the dataframe
  network_df <- vertex_strength_df %>%
    left_join(., vertex_betweenness_df, by = "Ind") %>%
    left_join(., vertex_centrality_df, by = "Ind") %>%
    left_join(., vertex_cc_df, by = "Ind") %>%
    left_join(., vertex_reach_df, by = "Ind")

  return(network_df)
}

#' Calculate various individual social network parameters based on DSI values
#' for multiple years and groups
#'
#' Determines Strength, Betweenness, Eigenvector centrality, Clustering
#' coefficient (cc), Reach. Works with weighted, but undirected edges. Works
#' with several groups and years
#'
#' @param edge_table Table with two individual and one weight column. Also
#'   requires columns "GroupCode" and "YearOf" or will fail.
#' @inheritParams calc_sna_pars_dyad
#'
#' @export
#'
#' @examples
#'

calc_sna_pars_dyad_mymg <- function(edge_table, indA_col = "IndA", indB_col = "IndB", weight_col = "DSI") {

  groups <- unique(edge_table$GroupCode)
  years <- unique(edge_table$YearOf)

  # Make sure that sna_par_df doesn't exist in the local environment (can probably be removed)
  if(exists("sna_par_df")) rm(sna_par_df)

  # Then, calculate sna parameters separately for each year and group
  for(i in groups){
    for(j in years){
      temp_df <- edge_table[edge_table$GroupCode == i & edge_table$YearOf == j, ]
      if(nrow(temp_df) == 0) next
      sna_par_df_temp <- calc_sna_pars_dyad(temp_df)
      sna_par_df_temp$GroupCode <- i
      sna_par_df_temp$YearOf <- j
      if(exists("sna_par_df")){
        sna_par_df <- bind_rows(sna_par_df, sna_par_df_temp)
        message(paste(i, j, sep = " - ", " added to df"))
      } else {
        sna_par_df <- sna_par_df_temp
        message(paste(i, j, sep = " - ", " used to create df"))
      }}}
return(sna_par_df)
}

# `reach` function from the CePa package:
# https://github.com/jokergoo/CePa
#
# The largest reach centrality is calculated as ``max(d(w, v))`` where ``d(w, v)`` is the
# length of the shortest path from node ``w`` to node ``v``.
#
# == author
# Zuguang Gu <z.gu@dkfz.de>
#
# == example
# require(igraph)
# pathway = barabasi.game(200)
# reach(pathway)
reach <-  function(graph, weights = igraph::E(graph)$weight, mode=c("all", "in", "out")) {
  mode = match.arg(mode)[1]

  sp = igraph::shortest.paths(graph, weights = weights, mode = mode)
  s = apply(sp, 1, function(x) {
    if(all(x == Inf)) {
      return(0)
    }
    else {
      return(max(x[x != Inf]))
    }
  })
  return(s)
}



#' Calculate individual strength
#'
#' Determines Strength
#'
#' @param rate_table Table with two individual ("IndA", "IndB") and one "Rate" column. Each dyad
#'   should be listed twice (once in each direction)
#'
#' @export
#'
#' @examples
#'

get_strength_from_rates <- function(rate_table){
  df <- select(rate_table, from = IndA, to = IndB, weight = Rate)
  inds <- unique(c(df$from, df$to))

  net <- igraph::graph_from_data_frame(d = df, vertices = inds, directed = T)

  # remove edges with no weight
  net.pos <- net - igraph::E(net)[igraph::E(net)$weight == 0]
  net.weights <- igraph::E(net.pos)$weight
  net.vertices <- igraph::V(net.pos)

  # To calculate strength, use "out". The 'rate' indicates how much each individual received/gave.
  vertex_strength <- igraph::strength(net.pos,
                                      weights = net.weights,
                                      mode = "out")

  return(tibble(Ind = names(vertex_strength), strength = vertex_strength))
}


#' Calculate individual strength for multiple years and groups
#'
#' Determines Strength
#' Works with several groups and years
#'
#' @param rate_table Table with two individual ("IndA", "IndB") and one "Rate" column. Each dyad
#'   should be listed twice (once in each direction) per year. Also requires columns "GroupCode" and "YearOf" or will fail.
#'
#' @export
#'
#' @examples
#'
get_strength_from_rates_mymg <- function(rate_table){
  groups <- unique(rate_table$GroupCode)
  years <- unique(rate_table$YearOf)

  # Make sure that strength_df doesn't exist in the local environment (can probably be removed)
  if(exists("strength_par_df")) rm(strength_par_df)

  # Then, calculate sna parameters separately for each year and group
  for(i in groups){
    for(j in years){
      temp_df <- rate_table[rate_table$GroupCode == i & rate_table$YearOf == j, ]
      if(nrow(temp_df) == 0) next
      strength_df_temp <- get_strength_from_rates(temp_df)
      strength_df_temp$GroupCode <- i
      strength_df_temp$YearOf <- j
      if(exists("strength_df")){
        strength_df <- bind_rows(strength_df, strength_df_temp)
        message(paste(i, j, sep = " - ", " added to df"))
      } else {
        strength_df <- strength_df_temp
        message(paste(i, j, sep = " - ", " used to create df"))
      }}}
  return(strength_df)
}




#' Calculate multiple estimates of Social Differentiation according to Whitehead
#' 2008 from a list of time periods, each with observed and permuted data.
#'
#' For details, see `get_S_from_matrix()` and `get_S_from_array()`. Difference
#' is that this functions takes a list of all time periods, each with observed
#' (in a matrix) and permuted data (in a 3d-array).
#'
#'
#'
#' @param dyadic_obs_list List created by `get_dyadic_observations()`
#' @param perm.start.layer Integer specifying the first permutation to be
#'   included
#' @param perm.by.layer Integer specifying the the permutations to be inlcluded
#'   (`by` argument in the `seq()` function)
#'
#' @inheritParams get_S_from_matrix
#' @inheritParams get_S_from_array
#'
#' @export
#'
#' @examples
#'

get_S_from_timeperiod_list <- function(dyadic_obs_list,
                                       perm.start.layer = 1, perm.by.layer = NULL,
                                       initial.values, lower, upper,
                                       use.parallel = FALSE, cores = 1) {

  ### Create empty list for all time periods, copy attributes from original
  ### list,  add additional attributes.
  S_list <- vector("list", length = length(dyadic_obs_list))
  attributes(S_list) <- attributes(dyadic_obs_list)
  attr(S_list, "perm.start.layer") <- perm.start.layer
  attr(S_list, "perm.by.layer") <- perm.by.layer

  ### Calculate S (and other parameters) for all time periods, including
  ### observed and permuted data
  for(i in 1:length(names(dyadic_obs_list))){
    cat("\n", names(dyadic_obs_list)[[i]])

    ### Prepare matrixes/arrays:

    # d (i_j_sum_matrix), x observed (i_j_together_observed_matrix), and x
    # permuted (i_j_together_permuted_array)

    # d
    i_j_sum_matrix <- dyadic_obs_list[[i]]$n_Ind_while_NN_present
    # Replace NAs by 0 and summarize both sides of matrix to make it symmetrical
    i_j_sum_matrix[is.na(i_j_sum_matrix)] <- 0
    i_j_sum_matrix <- i_j_sum_matrix + t(i_j_sum_matrix)

    # x observed
    i_j_together_observed_matrix <- dyadic_obs_list[[i]]$n_Ind_NN_together_observed
    # Replace NAs by 0 and summarize both sides of matrix to make it symmatrical
    i_j_together_observed_matrix[is.na(i_j_together_observed_matrix)] <- 0
    i_j_together_observed_matrix <- i_j_together_observed_matrix + t(i_j_together_observed_matrix)

    # x permuted
    # 1. Limit to a smaller number of permutations to have a reasonable
    # computation time: Only start with matrix perm.start.layer, and then only
    # keep every matrix perm.by.layer
    i_j_together_permuted_array <- dyadic_obs_list[[i]]$n_Ind_NN_together_permuted
    perms_n <- dim(i_j_together_permuted_array)[[3]]
    if(is.null(perm.by.layer)){
      perms_included <- perm.start.layer:perms_n
    } else {
      perms_included <- seq(from = perm.start.layer,
                            to = perms_n,
                            by = perm.by.layer)
    }
    i_j_together_permuted_array <- i_j_together_permuted_array[,,perms_included]
    # 2. Replace NAs by 0 and summarize both sides of matrix to make it symmetrical
    i_j_together_permuted_array[is.na(i_j_together_permuted_array)] <- 0
    # 3. Summarize both sides of matrix to make it symmetrical
    for(j in 1:dim(i_j_together_permuted_array)[[3]]){
      i_j_together_permuted_array[,,j] <- i_j_together_permuted_array[,,j] + t(i_j_together_permuted_array[,,j])
    }

    ### Get S, mu, logLik for observed data

    # This function takes symmetrical matrices and only uses one half of this
    # matrix to estimate social differentiation
    cat("\ncalculate S for observed relationships")
    S_list[[i]]$observed <- get_S_from_matrix(i_j_sum = i_j_sum_matrix,
                                              i_j_together = i_j_together_observed_matrix,
                                              initial.values = initial.values,
                                              lower = lower, upper = upper)

    ### Get S, mu, logLik for permuted data

    # This function takes an array of symmetrical matrices and does the calculation on each matrix
    cat("\ncalculate S for permuted relationships")
    S_list[[i]]$permuted <- get_S_from_array(i_j_sum = i_j_sum_matrix,
                                             i_j_together_array = i_j_together_permuted_array,
                                             initial.values = initial.values,
                                             lower = lower, upper = upper,
                                             use.parallel = use.parallel,
                                             cores = cores)

    # Add the names of used permutationss
    S_list[[i]]$permuted$Permutation_i <- dimnames(i_j_together_permuted_array)$Permutation_i
  }
  return(S_list)
}


#' Calculate Estimate of Social Different according to Whitehead 2008
#'
#' This function uses a max. likelhood to separate the variance of true association indices and the sampling variance.
#' It's described in "Precision and power in the analysis of social structure using associations" (Whitehead, 2008, Animal Behaviour) and
#' "Whitehead, H. 2008. Analyzing Animal Societies: Quantitative Methods for Vertebrate Social Analysis. Chicago: University of Chicago Press."
#'
#' @param i_j_sum Matrix with values i_j_sum[ij] used as denominator of the estimated association index (d in the Whitehead paper)
#' @param i_j_together Matrix with values i_j_together[ij] number of observations of individuals i and j together (x in the Whitehead paper)
#' @param initial.values Initial values for the parameters to be optimized over (`par` argument in `optim``)
#' @param lower lower bounds for parameters
#' @param upper upper bounds for parameters
#'
#' @export
#'
#' @examples
#'
get_S_from_matrix <- function(i_j_sum, i_j_together, initial.values = c(0.5, 0.5),
                              lower = c(0.01, 0.01), upper = c(1, 10)){
  # Check if d and x are both matrices
  if(!is.matrix(i_j_sum)) stop("i_j_sum is not a matrix")
  if(!is.matrix(i_j_together)) stop("i_j_together is not a matrix")

  # Change i_j_sum (d) and i_j_together (x) into vectors of 'distances'
  dat <- data.frame(dvec = as.vector(as.dist(i_j_sum)),
                    xvec = as.vector(as.dist(i_j_together)))

  res <- try(optim(par = initial.values,
                   calc_LL_S,
                   data = dat,
                   method = "L-BFGS-B",
                   lower = lower, upper = upper),
             silent = TRUE)

  if(class(res) == "try-error"){
    S_mu_logLik <- list(mu = NA,
                        S = NA,
                        logLik = NA,
                        error_message = attr(res, "condition")$message)
  } else {
    S_mu_logLik <- list(mu = res$par[1],
                        S = res$par[2],
                        logLik = res$value,
                        error_message = NA)
  }

  return(S_mu_logLik)
}


#' Calculate multiple estimates of Social Differentiation according to Whitehead
#' 2008 from an array.
#'
#' For details, see `get_S_from_matrix()`. Difference is that this functions
#' takes an 3d-array, where each slice is an matrix of observed association,
#' e.g. for permuted networks. Can use parallization, which is highly
#' recommended.
#'
#' @param i_j_together_array Array with matrices with values i_j_together[ij] number of
#'   observations of individuals i and j together (x in the Whitehead paper)
#' @param use.parallel TRUE/FALSE
#' @param cores Number of cores to be used for parallization
#'
#' @inheritParams get_S_from_matrix
#'
#' @export
#'
#' @examples
#'

get_S_from_array <- function(i_j_sum, i_j_together_array,
                             initial.values = c(0.5, 0.5),
                             lower = c(0.01, 0.01), upper = c(1, 10),
                             use.parallel = FALSE, cores = 1){
  # Check if i_j_sum is a matrix and i_j_together_array a 3d array
  if(!is.matrix(i_j_sum)) {
    stop("i_j_sum is not a matrix")}
  if(!is.array(i_j_together_array) | length(dim(i_j_together_array)) != 3) {
    stop("i_j_together_array is not a 3d array")}

  # Prepare (empty) list
  S_list <- list(mu = NA, S = NA, logLik = NA, error_message = NA)

  ### Calculate mu, S, and logLik for all layers of the array
  # Without parallelization
  if(!use.parallel){
    cat("\nWithout parallelization")
    for(i in 1:dim(i_j_together_array)[[3]]){
      cat("\n", i, " out of ", dim(i_j_together_array)[[3]])
      i_j_together_matrix <- i_j_together_array[,,i]
      temp_S <- get_S_from_matrix(i_j_sum = i_j_sum,
                                  i_j_together = i_j_together_matrix,
                                  initial.values = initial.values,
                                  lower = lower,
                                  upper = upper)

      S_list$mu[[i]] <- temp_S$mu
      S_list$S[[i]] <- temp_S$S
      S_list$logLik[[i]] <- temp_S$logLik
      S_list$error_message[[i]] <- temp_S$error_message

    }
    # With parallization
  } else {
    cat("\nWith parallelization - progress can be looked up in log.txt\n")
    require(doParallel)
    require(foreach)
    require(doSNOW)
    cl <- makeCluster(cores)
    registerDoSNOW(cl)

    # To get updates while running the loop, create a log-files
    writeLines(c(""), "log.txt")
    # Then run the loop
    temp_S <- foreach(i = 1:dim(i_j_together_array)[[3]],
                      .combine = "rbind",
                      .packages = c("daginR")) %dopar%  {
                        # Write progress in log-file (which can be opened during loop)
                        cat(paste0("\n", i, " out of ", dim(i_j_together_array)[[3]]),
                            file = "log.txt", append = TRUE)

                        # Extract matrix from array how often individuals were together
                        i_j_together_matrix <- i_j_together_array[,,i]

                        # Estimate mu and S
                        temp_S <- get_S_from_matrix(i_j_sum = i_j_sum,
                                                    i_j_together = i_j_together_matrix,
                                                    initial.values = initial.values,
                                                    lower = lower, upper = upper)
                        return(temp_S)
                      }

    stopCluster(cl)
    # Add data from parallel processing to list
    S_list$mu <- unlist(temp_S[,"mu"], use.names = FALSE)
    S_list$S <- unlist(temp_S[,"S"], use.names = FALSE)
    S_list$logLik <- unlist(temp_S[,"logLik"], use.names = FALSE)
    S_list$error_message <- unlist(temp_S[,"error_message"], use.names = FALSE)
  }
  return(S_list)
}

#' Function to calculate likelihood of mu and S (for `get_S``)
#'
#' For details, see `get_S`
#'
#' @param par vector with the two parameters mu and S
#' @param data Data frame with columns dvec and xvec
#'
#' @export
#'

calc_LL_S <- function(par, data){
  mu <- par[1]
  S <- par[2]

  # Define the integrand function
  integrand <- function(alpha){
    alpha^x * (1 - alpha)^(d - x) * dbeta(alpha, shape1 = beta1, shape2 = beta2)
  }

  # Define beta1 and beta2
  beta1 <- mu *     ((1-mu) / (mu * S^2) - 1)
  beta2 <- (1-mu) * ((1-mu) / (mu * S^2) - 1)

  # If one of the parameters is 0 or negative, the beta distribution is not defined.
  # In that case, set the logLik to a very large number
  if(beta1 <= 0 | beta2 <= 0){
    logLik <- 1e10
  } else {
    # Else, calculate the likelihood with the beta parameters
    # Remove all lines with d = 0 (i.e. individual never recorded while both were present)
    data <- data[data$dvec != 0,]
    ## Then solve integral numerically from 0 to 1 for all dyads (d and x)
    n <- length(data$dvec)
    integrated <- rep(NA, n)

    for(i in 1:n){
      d <- data$dvec[i]
      x <- data$xvec[i]
      integrated[i] <- (integrate(integrand, lower = 0, upper = 1, rel.tol = 1e-8, abs.tol = 0))$value
    }

    # If any result here is 0, the product of all integrated values would be 0.
    # Then, the log of the product would approach -Inf
    # Therefore, set the logLik to very large number.
    if(any(integrated == 0)){
      logLik <- 1e10
      # If not, calculate the negative sum of log(integrated) (which is the same as log(prod(integrated)))
    } else {
      logLik <- -sum(log(integrated))
    }
  }
  return(logLik)
}

#' Function to calculate Coefficients of variation (CV) of association indices
#' Whitehead (2008) suggests to use his estimate of S, instead, but this
#' approach is also commonly used.
#'
#' This function uses a list of the Simple Ratio Association Index (or any other
#' index) for all included time periods. It then calculate the CV for the
#' observed data and for matrices of permuted social networks
#'
#' @param simple_ratio_list List created with `get_simple_ratios()``
#' @param perm.by.layer In case not all layers (matrices) of permuted data
#'   should be included. This is the `by` argument of the `seq()` function.
#' @param zeros.rm Should all zeros be removed (i.e. dyads that were never
#'   observed together). Could be problematic with permuted networks.
#' @param na.rm Should NAs be removed. This argument is probably redudant as all
#'   NAs are set to 0 (which is correct in this case). But if set to FALSE, it
#'   controls if all matrices were well constructed.
#'
#' @export
#'

get_CVs_from_timeperiod_list <- function(simple_ratio_list, perm.by.layer = NULL,
                                         zeros.rm = FALSE, na.rm = FALSE) {
  # Create empty list
  cv_list <- vector("list", length = length(simple_ratio_list))
  attributes(cv_list) <- attributes(simple_ratio_list)
  attr(cv_list, "perm.by.layer") <- perm.by.layer
  attr(cv_list, "na_removed") <- na.rm
  attr(cv_list, "zeros_removed") <- zeros.rm

  # Calculate CV for all time periods and observed and permuted data
  for(i in 1:length(names(simple_ratio_list))){
    cat("\n", names(simple_ratio_list)[[i]])

    # Define function to calculate CV from vector of values
    matrix_to_cv <- function(symmetrical_matrix, zeros.rm = FALSE, na.rm = FALSE){
      vec <- as.vector(as.dist(symmetrical_matrix))
      if(zeros.rm) vec <- vec[vec != 0]
      if(na.rm) vec <- vec[!is.na(vec)]
      return(sd(vec)/mean(vec))
    }

    ### 1. Calculate CV for observed data
    matrix_observed <- simple_ratio_list[[i]]$simple_ratio_observed
    cv_list[[i]]$observed <- list(cv = matrix_to_cv(matrix_observed, zeros.rm, na.rm))

    ### 2. Calculate CV for permuted data
    # Create vector of permutations to be included
    perms_n <- dim(simple_ratio_list[[i]]$simple_ratio_permuted)[[3]]
    if(is.null(perm.by.layer)){
      perms_included <- 1:perms_n
    } else {
      perms_included <- seq(from = 1, to = perms_n, by = perm.by.layer)
    }
    array_permuted <- simple_ratio_list[[i]]$simple_ratio_permuted[,,perms_included]
    cv_list[[i]]$permuted <- list(cv = apply(array_permuted, MARGIN = 3, FUN = matrix_to_cv,
                                             zeros.rm, na.rm))
  }
  return(cv_list)
}


#' Get matrix correlations from timeperiod list
#'
#' Function to calculate correlation coefficients of all two consecutive
#' networks within a list of timeperiods, each with observed and permuted
#' association indices.
#'
#'
#' @param simple_ratio_list List created with `get_simple_ratios()`
#' @inheritParams calc_matrix_correlation
#' @inheritParams stats::cor
#'
#' @export
#'

get_matrix_cors_from_timeperiod_list <- function(simple_ratio_list, zeros_remove = FALSE,
                                                 method = c("pearson", "kendall", "spearman")){
  method <- match.arg(method)
  timeperiods <- length(simple_ratio_list)
  # Create empty list for correlation values
  corr_list = vector("list", length = (timeperiods - 1))
  # Copy attributes from simple_ratio list (excluding timeperiod information)
  attributes_to_copy <- names(attributes(simple_ratio_list))
  attributes_to_copy <- attributes_to_copy[!attributes_to_copy %in% c("names", "timeperiod_start", "timeperiod_end")]
  attributes(corr_list) <- attributes(simple_ratio_list)[attributes_to_copy]
  # Set additional attributes
  corr_list <- structure(corr_list,
                         corr_method = method,
                         names = paste0("t", 1:(timeperiods - 1), "_t", 2:timeperiods),
                         t1_start = attr(simple_ratio_list, "timeperiod_start")[1:(timeperiods - 1)],
                         t1_end = attr(simple_ratio_list, "timeperiod_end")[1:(timeperiods - 1)],
                         t2_start = attr(simple_ratio_list, "timeperiod_start")[2:timeperiods],
                         t2_end = attr(simple_ratio_list, "timeperiod_end")[2:timeperiods])


  # For all consecutive periods of time periods, calculate the correlation between networks
  for(i in 1:(timeperiods - 1)){
    cat("\nProcess", names(corr_list)[[i]])
    # Get the pair of timeperiods i and i+1
    t1 <- simple_ratio_list[[i]]
    t2 <- simple_ratio_list[[i + 1]]
    # Calculate correlation of observed networks
    corr_list[[i]]$observed <- as.numeric(calc_matrix_correlation(t1$simple_ratio_observed,
                                                                  t2$simple_ratio_observed,
                                                                  zeros_remove = zeros_remove,
                                                                  method = method)$statistic)
    # Calculate correlations of permuted networks
    corr_list[[i]]$permuted <- sapply(1:dim(t1$simple_ratio_permuted)[[3]],
                                      FUN = function(x)
                                        calc_matrix_correlation(t1$simple_ratio_permuted[,,x],
                                                                t2$simple_ratio_permuted[,,x],
                                                                zeros_remove = zeros_remove,
                                                                method = method)$statistic)
  }
  return(corr_list)
}

#' Calculates the Partner Preference Index (PPI) based on two network matrices.
#'
#' This function calculate the Partner Preference Index (PPI) as described by
#' Silk, Cheney and Seyfarth (2013, Evolutionary Anthropology):
#'
#' PPI = (2S -  U)/(2S - S - X) with S = number of partner rank slots being
#' evaluated (here `n_Slots`), U = Number of different top partners that the
#' individual had in those years, and X = Number of top partners present in t1,
#' but not t2:
#'
#' @param m1,m2 The two matrices.
#' @param n_Slots number of top partners evaluated (S in Silk et al. 2013)
#' @param if FALSE, returns a single value. If TRUE, a dataframe with details
#'
#' @export
#'
calc_ppi <- function(m1, m2, n_Slots = 3, details = FALSE){
  # Determine distinct individuals in t1 and t2
  # Only evaluate partner preference for individuals who are present in t1 and t2
  # (Not possible to evaluate preference for other individuals)
  inds_t1 <- rownames(m1)
  inds_t2 <- rownames(m2)
  inds_to_evaluate <- intersect(inds_t1, inds_t2)

  # For each individual who is present in t1 and t2 determine ppi"
  # Create empty dataframe
  cols_to_create <- c("ind",
                      paste0("top_t1_", 1:n_Slots),
                      paste0("top_t2_", 1:n_Slots),
                      "n_top_t1t2",
                      "n_top_t1_not_present_t2")
  df <- data.frame(matrix(ncol = length(cols_to_create), nrow = length(inds_to_evaluate)))
  colnames(df) <- cols_to_create

  # For the calculation of PPI = = (2S -  U)/(2S - S - X)
  # S = n_Slots
  # U = union(top_t1, top_t2) with top_t1, top_t2 = top_n_partners in t1/t2
  # X = n_top_t1_not_present_in_t2
  for(i in seq_along(inds_to_evaluate)){
    ind = inds_to_evaluate[i]
    top_t1 <- sort(m1[ind,], decreasing = TRUE)[1:n_Slots]
    top_t2 <- sort(m2[ind,], decreasing = TRUE)[1:n_Slots]
    # Check if n_Slots top partners were available
    if(length(top_t1) != n_Slots | length(top_t2) != n_Slots) stop("Length of top partners not equal to S")

    top_t1t2 <- union(names(top_t1), names(top_t2))
    top_t1_not_present_t2 <- setdiff(names(top_t1), inds_t2)

    vals <- c(ind, names(top_t1), names(top_t2), length(top_t1t2),
              length(top_t1_not_present_t2))
    if(length(names(df)) != length(vals)) stop("(internal): colnames differ between template and new df")
    for(j in 1:length(vals)){
      df[i,j] <- vals[j]
    }
  }
  # Now, calculate PPI
  df$n_top_t1t2 <- as.numeric(df$n_top_t1t2)
  df$n_top_t1_not_present_t2 <- as.numeric(df$n_top_t1_not_present_t2)
  df$ppi <- (2 * n_Slots - df$n_top_t1t2) /
    (2 * n_Slots - n_Slots - df$n_top_t1_not_present_t2)

  # Return entire dataframe, or only the group average
  if(details) {
    return(df)
  } else {
    mean(df$ppi, na.rm = TRUE)
  }
}

#' Get the Partner Preference Index (PPI) from timeperiod list with association indices.
#'
#' Function to calculate the Partner Preference Index (PPI) of all two consecutive
#' networks within a list of timeperiods, each with observed and permuted
#' association indices.
#'
#'
#' @param simple_ratio_list List created with `get_simple_ratios()`
#'
#' @inheritParams calc_ppi
#'
#'
#' @export
#'
get_ppi_from_timeperiod_list <- function(simple_ratio_list, n_slots = 3){
  timeperiods <- length(simple_ratio_list)
  # Create empty list for correlation values
  ppi_list = vector("list", length = (timeperiods - 1))
  # Copy attributes from simple_ratio list (excluding timeperiod information)
  attributes_to_copy <- names(attributes(simple_ratio_list))
  attributes_to_copy <- attributes_to_copy[!attributes_to_copy %in% c("names", "timeperiod_start", "timeperiod_end")]
  attributes(ppi_list) <- attributes(simple_ratio_list)[attributes_to_copy]
  # Set additional attributes
  ppi_list <- structure(ppi_list,
                        n_slots = n_slots,
                        names = paste0("t", 1:(timeperiods - 1), "_t", 2:timeperiods),
                        t1_start = attr(simple_ratio_list, "timeperiod_start")[1:(timeperiods - 1)],
                        t1_end = attr(simple_ratio_list, "timeperiod_end")[1:(timeperiods - 1)],
                        t2_start = attr(simple_ratio_list, "timeperiod_start")[2:timeperiods],
                        t2_end = attr(simple_ratio_list, "timeperiod_end")[2:timeperiods])


  # For all consecutive periods of time periods, calculate the correlation between networks
  for(i in 1:(timeperiods - 1)){
    cat("\nProcess", names(ppi_list)[[i]])
    # Get the pair of timeperiods i and i+1
    t1 <- simple_ratio_list[[i]]
    t2 <- simple_ratio_list[[i + 1]]
    # Calculate correlation of observed networks
    ppi_list[[i]]$observed <- calc_ppi(t1$simple_ratio_observed, t2$simple_ratio_observed)
    # Calculate correlations of permuted networks
    ppi_list[[i]]$permuted <- sapply(1:dim(t1$simple_ratio_permuted)[[3]],
                                     FUN = function(x)
                                       calc_ppi(t1$simple_ratio_permuted[,,x],
                                                t2$simple_ratio_permuted[,,x]))
  }
  return(ppi_list)
}
