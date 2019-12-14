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

  # Change d and x into vectors of 'distances'
  dat <- data.frame(dvec = as.vector(as.dist(i_j_sum)),
                    xvec = as.vector(as.dist(i_j_together)))

  res <- optim(par = initial.values,
               calc_LL_S,
               data = dat,
               method = "L-BFGS-B",
               lower = lower, upper = upper)

  S_mu_logLik <- list(mu = res$par[1],
                      S = res$par[2],
                      logLik = res$value)

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
  S_list <- list(mu = NA, S = NA, logLik = NA)

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
    }
    # With parallization
  } else {
    cat("\nWith parallelization - progress can be looked up in log.txt")
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
    S_list$mu <- as.numeric(temp_S[,"mu"])
    S_list$S <- as.numeric(temp_S[,"S"])
    S_list$logLik <- as.numeric(temp_S[,"logLik"])
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
