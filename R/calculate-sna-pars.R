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
#' @inheritParams get_strength_from_rates
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

