# Randomized recursive adjustment-set search. The local DAG is built by the
# existing network layer; this file only searches and validates adjustment sets.

tf_recursive_adjustment_build_checker <- function(dag_edges, exposure, outcome) {
  vertices <- sort(unique(c(dag_edges$tf, dag_edges$target, exposure, outcome)))
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

tf_recursive_adjustment_single_start <- function(
  start_id,
  start_seed,
  initial_set,
  fallback_sets,
  checker,
  distance_matrix,
  nodes,
  allowed_mask,
  max_iterations,
  max_pair_checks,
  max_ancestor_candidates,
  max_pairs_per_ancestor
) {
  set.seed(as.integer(start_seed))
  started <- Sys.time()
  validation_checks <- 1L
  if (!checker$is_valid(initial_set)) {
    valid_fallback <- NULL
    for (candidate in fallback_sets) {
      validation_checks <- validation_checks + 1L
      if (checker$is_valid(candidate)) {
        valid_fallback <- candidate
        break
      }
    }
    if (is.null(valid_fallback)) {
      stop(
        "Neither the outcome-parent nor exposure-parent fallback is a valid adjustment set in the local DAG.",
        call. = FALSE
      )
    }
    initial_set <- valid_fallback
  }

  current <- sort(unique(initial_set))
  lineage <- stats::setNames(lapply(current, function(node) node), current)
  operations <- list()
  operation_index <- 0L
  pair_checks <- 0L
  record_operation <- function(iteration, action, removed, added, before, after) {
    operation_index <<- operation_index + 1L
    operations[[operation_index]] <<- data.frame(
      iteration = iteration,
      action = action,
      removed = paste(removed, collapse = ";"),
      added = if (length(added)) paste(added, collapse = ";") else "",
      size_before = length(before),
      size_after = length(after),
      stringsAsFactors = FALSE
    )
  }

  iteration <- 0L
  repeat {
    iteration <- iteration + 1L
    if (iteration > max_iterations || pair_checks >= max_pair_checks) break
    changed <- FALSE

    for (node in sample(current, length(current))) {
      proposed <- setdiff(current, node)
      validation_checks <- validation_checks + 1L
      if (checker$is_valid(proposed)) {
        before <- current
        current <- proposed
        lineage[[node]] <- NULL
        record_operation(
          iteration, "drop_redundant", node, character(0), before, current
        )
        changed <- TRUE
        break
      }
    }
    if (changed) next
    if (length(current) < 2L) break

    ancestor_nodes <- nodes[allowed_mask]
    coverage <- rowSums(is.finite(
      distance_matrix[ancestor_nodes, current, drop = FALSE]
    ))
    keep <- coverage >= 2L
    ancestor_nodes <- ancestor_nodes[keep]
    coverage <- coverage[keep]
    if (!length(ancestor_nodes)) break
    mean_distance <- vapply(ancestor_nodes, function(ancestor) {
      values <- distance_matrix[ancestor, current]
      mean(values[is.finite(values)])
    }, numeric(1))
    ancestor_nodes <- ancestor_nodes[
      order(-coverage, mean_distance, stats::runif(length(ancestor_nodes)), ancestor_nodes)
    ]
    ancestor_nodes <- utils::head(ancestor_nodes, max_ancestor_candidates)

    for (ancestor in ancestor_nodes) {
      if (pair_checks >= max_pair_checks) break
      descendants <- current[is.finite(distance_matrix[ancestor, current])]
      removable <- setdiff(descendants, ancestor)
      if (length(descendants) < 2L || !length(removable)) next

      proposed <- sort(unique(c(setdiff(current, descendants), ancestor)))
      pair_checks <- pair_checks + 1L
      validation_checks <- validation_checks + 1L
      if (length(proposed) < length(current) && checker$is_valid(proposed)) {
        before <- current
        absorbed <- unique(unlist(
          lineage[unique(c(descendants, ancestor))], use.names = FALSE
        ))
        for (node in descendants) lineage[[node]] <- NULL
        lineage[[ancestor]] <- absorbed
        current <- proposed
        record_operation(
          iteration, "replace_descendants_with_ancestor", removable,
          ancestor, before, current
        )
        changed <- TRUE
        break
      }

      pair_matrix <- utils::combn(descendants, 2L)
      pair_distance <- distance_matrix[ancestor, pair_matrix[1L, ]] +
        distance_matrix[ancestor, pair_matrix[2L, ]]
      pair_order <- order(pair_distance, stats::runif(length(pair_distance)))
      pair_order <- utils::head(pair_order, max_pairs_per_ancestor)
      for (column in pair_order) {
        if (pair_checks >= max_pair_checks) break
        left <- pair_matrix[1L, column]
        right <- pair_matrix[2L, column]
        if (ancestor %in% c(left, right)) next
        proposed <- sort(unique(c(setdiff(current, c(left, right)), ancestor)))
        if (length(proposed) >= length(current)) next
        pair_checks <- pair_checks + 1L
        validation_checks <- validation_checks + 1L
        if (!checker$is_valid(proposed)) next

        before <- current
        absorbed <- unique(c(lineage[[left]], lineage[[right]], lineage[[ancestor]]))
        lineage[[left]] <- NULL
        lineage[[right]] <- NULL
        lineage[[ancestor]] <- absorbed
        current <- proposed
        record_operation(
          iteration, "replace_pair_with_common_ancestor", c(left, right),
          ancestor, before, current
        )
        changed <- TRUE
        break
      }
      if (changed) break
    }
    if (!changed) break
  }

  validation_checks <- validation_checks + 1L
  final_valid <- checker$is_valid(current)
  operations_table <- if (length(operations)) {
    do.call(rbind, operations)
  } else {
    data.frame(
      iteration = integer(0), action = character(0), removed = character(0),
      added = character(0), size_before = integer(0), size_after = integer(0),
      stringsAsFactors = FALSE
    )
  }
  list(
    start_id = as.integer(start_id),
    seed = as.integer(start_seed),
    final_set = sort(current),
    lineage = lineage,
    operations = operations_table,
    pair_checks = pair_checks,
    validation_checks = validation_checks,
    iterations = iteration,
    terminated_by_limit = iteration > max_iterations || pair_checks >= max_pair_checks,
    final_valid = final_valid,
    elapsed_seconds = as.numeric(difftime(Sys.time(), started, units = "secs")),
    backdoor_edges = checker$backdoor_edge_count
  )
}

tf_recursive_adjustment_run_batch <- function(start_ids, context) {
  checker <- tf_recursive_adjustment_build_checker(
    context$dag_edges, context$tf, context$gene
  )
  lapply(start_ids, function(start_id) {
    tf_recursive_adjustment_single_start(
      start_id = start_id,
      start_seed = context$seeds[[start_id]],
      initial_set = context$initial_set,
      fallback_sets = context$fallback_sets,
      checker = checker,
      distance_matrix = context$distance_matrix,
      nodes = context$nodes,
      allowed_mask = context$allowed_mask,
      max_iterations = context$max_iterations,
      max_pair_checks = context$max_pair_checks,
      max_ancestor_candidates = context$max_ancestor_candidates,
      max_pairs_per_ancestor = context$max_pairs_per_ancestor
    )
  })
}

tf_recursive_adjustment_score_paths <- function(
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
    path_edge_mean_confidence = if (all(usable) && any(usable)) {
      sum(details$path_confidence_sum) / sum(details$path_length)
    } else {
      -Inf
    },
    total_path_confidence = if (all(usable) && any(usable)) {
      sum(details$path_confidence_sum)
    } else {
      -Inf
    },
    total_path_edges = if (all(usable) && any(usable)) {
      sum(details$path_length)
    } else {
      Inf
    }
  )
}

tf_recursive_adjustment_search <- function(
  dag_edges,
  exposure,
  outcome,
  n_starts = 8L,
  cores = 4L,
  seed = 123L,
  eligible_nodes = NULL,
  max_iterations = 1000L,
  max_pair_checks = 200000L,
  max_ancestor_candidates = 500L,
  max_pairs_per_ancestor = 20L
) {
  positive_integer <- function(value, name) {
    if (is.null(value) || !length(value)) {
      stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
    }
    value <- suppressWarnings(as.integer(value[[1L]]))
    if (length(value) != 1L || is.na(value) || value < 1L) {
      stop(sprintf("`%s` must be a positive integer.", name), call. = FALSE)
    }
    value
  }
  n_starts <- positive_integer(n_starts, "n_starts")
  cores <- min(positive_integer(cores, "cores"), n_starts)
  if (is.null(seed) || !length(seed)) {
    stop("`seed` must be an integer.", call. = FALSE)
  }
  seed <- suppressWarnings(as.integer(seed[[1L]]))
  if (length(seed) != 1L || is.na(seed)) {
    stop("`seed` must be an integer.", call. = FALSE)
  }
  dag_edges <- dag_edges[
    !duplicated(dag_edges[, c("tf", "target")]), , drop = FALSE
  ]
  graph <- igraph::graph_from_data_frame(
    dag_edges[, c("tf", "target"), drop = FALSE], directed = TRUE
  )
  if (!igraph::is_dag(graph)) {
    stop("`dag_edges` must define a directed acyclic graph.", call. = FALSE)
  }
  nodes <- igraph::V(graph)$name
  if (!(exposure %in% nodes) || !(outcome %in% nodes)) {
    stop("Exposure and outcome must both occur in `dag_edges`.", call. = FALSE)
  }
  distance_matrix <- igraph::distances(graph, mode = "out")
  descendants <- nodes[
    is.finite(distance_matrix[exposure, ]) & distance_matrix[exposure, ] > 0
  ]
  direct_parents <- sort(unique(dag_edges$tf[dag_edges$target == outcome]))
  direct_parents <- setdiff(direct_parents, c(exposure, outcome, descendants))
  exposure_parents <- sort(unique(dag_edges$tf[dag_edges$target == exposure]))
  exposure_parents <- setdiff(
    exposure_parents, c(exposure, outcome, descendants)
  )
  exposure_ancestors <- nodes[is.finite(distance_matrix[, exposure])]
  on_backdoor_trek <- vapply(direct_parents, function(node) {
    any(is.finite(distance_matrix[exposure_ancestors, node]))
  }, logical(1))
  initial_set <- direct_parents[on_backdoor_trek]
  allowed_ancestors <- setdiff(
    unique(as.character(dag_edges$tf)), c(exposure, outcome, descendants)
  )
  if (!is.null(eligible_nodes)) {
    allowed_ancestors <- intersect(allowed_ancestors, as.character(eligible_nodes))
  }

  set.seed(seed)
  seeds <- sample.int(.Machine$integer.max, n_starts)
  context <- list(
    dag_edges = dag_edges,
    tf = exposure,
    gene = outcome,
    seeds = seeds,
    initial_set = initial_set,
    fallback_sets = list(direct_parents, exposure_parents),
    distance_matrix = distance_matrix,
    nodes = nodes,
    allowed_mask = nodes %in% allowed_ancestors,
    max_iterations = positive_integer(max_iterations, "max_iterations"),
    max_pair_checks = positive_integer(max_pair_checks, "max_pair_checks"),
    max_ancestor_candidates = positive_integer(
      max_ancestor_candidates, "max_ancestor_candidates"
    ),
    max_pairs_per_ancestor = positive_integer(
      max_pairs_per_ancestor, "max_pairs_per_ancestor"
    )
  )
  start_ids <- seq_len(n_starts)
  started <- Sys.time()
  if (cores > 1L && n_starts > 1L) {
    cluster <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cluster), add = TRUE)
    tfregact_library <- dirname(find.package("TFRegAct"))
    parallel::clusterExport(
      cluster, "tfregact_library", envir = environment()
    )
    parallel::clusterEvalQ(cluster, {
      .libPaths(c(tfregact_library, .libPaths()))
      library(TFRegAct)
      NULL
    })
    worker_context <- context
    parallel::clusterExport(cluster, "worker_context", envir = environment())
    batches <- split(start_ids, rep(seq_len(cores), length.out = n_starts))
    batch_results <- parallel::parLapply(cluster, batches, function(ids) {
      worker <- utils::getFromNamespace(
        "tf_recursive_adjustment_run_batch", "TFRegAct"
      )
      worker(ids, worker_context)
    })
    all_starts <- unlist(batch_results, recursive = FALSE)
    all_starts <- all_starts[
      order(vapply(all_starts, `[[`, integer(1), "start_id"))
    ]
  } else {
    all_starts <- tf_recursive_adjustment_run_batch(start_ids, context)
  }
  search_seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))

  score_cache <- new.env(parent = emptyenv())
  start_scores <- lapply(all_starts, function(item) {
    key <- paste(item$final_set, collapse = "\r")
    if (!exists(key, envir = score_cache, inherits = FALSE)) {
      assign(
        key,
        tf_recursive_adjustment_score_paths(
          dag_edges, item$final_set, outcome, exposure = exposure
        ),
        envir = score_cache
      )
    }
    get(key, envir = score_cache, inherits = FALSE)
  })
  multistart_summary <- do.call(rbind, lapply(seq_along(all_starts), function(i) {
    item <- all_starts[[i]]
    score <- start_scores[[i]]
    data.frame(
      start_id = item$start_id,
      seed = item$seed,
      final_adjustment_size = length(item$final_set),
      path_edge_mean_confidence = score$path_edge_mean_confidence,
      total_path_confidence = score$total_path_confidence,
      total_path_edges = score$total_path_edges,
      final_valid = item$final_valid,
      operations = nrow(item$operations),
      pair_checks = item$pair_checks,
      validation_checks = item$validation_checks,
      elapsed_seconds = item$elapsed_seconds,
      terminated_by_limit = item$terminated_by_limit,
      stringsAsFactors = FALSE
    )
  }))
  valid <- which(multistart_summary$final_valid)
  if (!length(valid)) {
    stop("None of the randomized starts produced a valid adjustment set.", call. = FALSE)
  }
  rank_order <- order(
    multistart_summary$final_adjustment_size[valid],
    -multistart_summary$path_edge_mean_confidence[valid],
    multistart_summary$start_id[valid]
  )
  winner_index <- valid[rank_order[[1L]]]
  list(
    winner = all_starts[[winner_index]],
    winner_score = start_scores[[winner_index]],
    all_starts = all_starts,
    multistart_summary = multistart_summary,
    initial_set = initial_set,
    direct_parents = direct_parents,
    exposure_parents = exposure_parents,
    n_starts = n_starts,
    cores = cores,
    seed = seed,
    search_seconds = search_seconds,
    path_rule = paste(
      "shortest directed path to outcome; among equally short paths,",
      "maximum cumulative edge confidence"
    )
  )
}
