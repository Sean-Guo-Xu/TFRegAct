# Graph validation and path-confidence summaries used by exact adjustment search.

tf_adjustment_build_checker <- function(dag_edges, exposure, outcome) {
  vertices <- sort(unique(c(dag_edges$tf, dag_edges$target, exposure, outcome)))
  original_graph <- igraph::graph_from_data_frame(
    dag_edges[, c("tf", "target"), drop = FALSE], directed = TRUE,
    vertices = data.frame(name = vertices)
  )
  forbidden <- setdiff(igraph::V(original_graph)$name[
    igraph::subcomponent(original_graph, exposure, mode = "out")
  ], exposure)
  backdoor_edges <- dag_edges[
    dag_edges$tf != exposure, c("tf", "target"), drop = FALSE
  ]
  backdoor_graph <- igraph::graph_from_data_frame(
    backdoor_edges,
    directed = TRUE,
    vertices = data.frame(name = vertices, stringsAsFactors = FALSE)
  )
  if (!igraph::is_dag(backdoor_graph)) {
    stop("The adjustment validator requires a directed acyclic graph.", call. = FALSE)
  }
  graph_nodes <- igraph::V(backdoor_graph)$name
  node_index <- stats::setNames(seq_along(graph_nodes), graph_nodes)
  parents <- lapply(igraph::as_adj_list(backdoor_graph, mode = "in"), as.integer)
  children <- lapply(igraph::as_adj_list(backdoor_graph, mode = "out"), as.integer)
  exposure_index <- unname(node_index[[exposure]])
  outcome_index <- unname(node_index[[outcome]])

  is_valid <- function(adjustment_set) {
    adjustment_set <- unique(as.character(adjustment_set))
    if (exposure %in% adjustment_set || outcome %in% adjustment_set) return(FALSE)
    if (any(adjustment_set %in% forbidden)) return(FALSE)
    if (length(setdiff(adjustment_set, graph_nodes))) return(FALSE)
    recursive_bayes_ball_cpp(
      parents,
      children,
      exposure_index,
      outcome_index,
      as.integer(unname(node_index[adjustment_set]))
    )
  }

  list(
    is_valid = is_valid,
    graph_nodes = graph_nodes,
    backdoor_edge_count = igraph::ecount(backdoor_graph)
  )
}

tf_adjustment_score_paths <- function(
  dag_edges,
  adjustment_set,
  outcome,
  exposure = NULL
) {
  edge_confidence <- suppressWarnings(as.numeric(dag_edges$confidence_score))
  edge_confidence[!is.finite(edge_confidence)] <- 0
  edge_direction <- find_adjustment_dagitty_effect_to_direction(dag_edges$effect)
  scored_edges <- data.frame(
    tf = as.character(dag_edges$tf),
    target = as.character(dag_edges$target),
    confidence = edge_confidence,
    direction = edge_direction,
    stringsAsFactors = FALSE
  )
  scored_edges <- scored_edges[
    order(scored_edges$tf, scored_edges$target, -scored_edges$confidence),
    , drop = FALSE
  ]
  scored_edges <- scored_edges[
    !duplicated(scored_edges[, c("tf", "target")]), , drop = FALSE
  ]
  # Adjustment relevance must be scored on the back-door graph. Otherwise a
  # common cause could receive credit through common_cause -> exposure ->
  # outcome, which is the causal route whose effect is being estimated.
  if (!is.null(exposure)) {
    scored_edges <- scored_edges[scored_edges$tf != exposure, , drop = FALSE]
  }
  graph_nodes_input <- sort(unique(c(
    dag_edges$tf, dag_edges$target, adjustment_set, outcome
  )))
  graph <- igraph::graph_from_data_frame(
    scored_edges[, c("tf", "target"), drop = FALSE],
    directed = TRUE,
    vertices = data.frame(name = graph_nodes_input, stringsAsFactors = FALSE)
  )
  graph_nodes <- igraph::V(graph)$name
  distance_to_outcome <- as.numeric(igraph::distances(
    graph, v = graph_nodes, to = outcome, mode = "out"
  )[, 1L])
  names(distance_to_outcome) <- graph_nodes
  best_sum <- stats::setNames(rep(-Inf, length(graph_nodes)), graph_nodes)
  next_node <- stats::setNames(rep(NA_character_, length(graph_nodes)), graph_nodes)
  next_direction <- stats::setNames(integer(length(graph_nodes)), graph_nodes)
  best_sum[[outcome]] <- 0
  finite_distance <- distance_to_outcome[is.finite(distance_to_outcome)]
  max_distance <- if (length(finite_distance)) max(finite_distance) else 0
  if (max_distance >= 1L) {
    for (step in seq_len(as.integer(max_distance))) {
      step_nodes <- names(distance_to_outcome)[distance_to_outcome == step]
      for (node in step_nodes) {
        candidates <- scored_edges[
          scored_edges$tf == node &
            distance_to_outcome[scored_edges$target] == step - 1L,
          , drop = FALSE
        ]
        if (!nrow(candidates)) next
        candidate_score <- candidates$confidence + best_sum[candidates$target]
        best <- order(-candidate_score, candidates$target)[[1L]]
        best_sum[[node]] <- candidate_score[[best]]
        next_node[[node]] <- candidates$target[[best]]
        next_direction[[node]] <- candidates$direction[[best]]
      }
    }
  }

  trace_path <- function(start) {
    node <- start
    path <- node
    directions <- integer(0)
    while (!identical(node, outcome) && !is.na(next_node[[node]])) {
      directions <- c(directions, next_direction[[node]])
      node <- next_node[[node]]
      path <- c(path, node)
    }
    direction <- if (!length(directions) || any(directions == 0L)) {
      0L
    } else {
      as.integer(prod(directions))
    }
    list(path = paste(path, collapse = " -> "), direction = direction)
  }
  traced <- lapply(adjustment_set, trace_path)
  path_length <- distance_to_outcome[adjustment_set]
  path_sum <- best_sum[adjustment_set]
  path_mean <- path_sum / path_length
  details <- data.frame(
    variable = adjustment_set,
    path_length = as.numeric(path_length),
    path_confidence_sum = as.numeric(path_sum),
    path_mean_confidence = as.numeric(path_mean),
    path_direction = vapply(traced, `[[`, integer(1), "direction"),
    path_effect = vapply(traced, function(value) {
      if (value$direction == 1L) "activation" else if (value$direction == -1L) {
        "repression"
      } else {
        "unknown"
      }
    }, character(1)),
    selected_path = vapply(traced, `[[`, character(1), "path"),
    stringsAsFactors = FALSE
  )
  usable <- is.finite(details$path_length) & details$path_length > 0 &
    is.finite(details$path_confidence_sum)
  list(
    details = details,
    path_edge_mean_confidence = if (!length(adjustment_set)) {
      0
    } else if (all(usable) && any(usable)) {
      sum(details$path_confidence_sum) / sum(details$path_length)
    } else {
      -Inf
    },
    total_path_confidence = if (!length(adjustment_set)) {
      0
    } else if (all(usable) && any(usable)) {
      sum(details$path_confidence_sum)
    } else {
      -Inf
    },
    total_path_edges = if (!length(adjustment_set)) {
      0
    } else if (all(usable) && any(usable)) {
      sum(details$path_length)
    } else {
      Inf
    }
  )
}
