# Exact minimum-cardinality adjustment under the package's Pearl back-door
# criterion. See van der Zander, Liskiewicz & Textor (2019), Artificial
# Intelligence 270:1-40, doi:10.1016/j.artint.2018.12.006.
# Forbidden/unmeasured nodes stay in the graph: they cannot be cut, not erased.
tf_minimum_adjustment_search <- function(dag_edges, exposure, outcome,
                                         eligible_nodes = NULL) {
  started <- proc.time()[["elapsed"]]
  if (!is.data.frame(dag_edges) || !all(c("tf", "target") %in% names(dag_edges))) {
    stop("`dag_edges` must contain tf and target columns.", call. = FALSE)
  }
  for (value in list(exposure, outcome)) {
    if (!is.character(value) || length(value) != 1L || is.na(value) || !nzchar(value)) {
      stop("Exposure and outcome must be nonempty single strings.", call. = FALSE)
    }
  }
  if (exposure == outcome) stop("Exposure and outcome must differ.", call. = FALSE)
  dag_edges$tf <- as.character(dag_edges$tf)
  dag_edges$target <- as.character(dag_edges$target)
  if (anyNA(dag_edges[, c("tf", "target")]) ||
      any(!nzchar(dag_edges$tf)) || any(!nzchar(dag_edges$target))) {
    stop("Graph endpoints must be nonempty and nonmissing.", call. = FALSE)
  }
  if (!"confidence_score" %in% names(dag_edges)) dag_edges$confidence_score <- rep(0, nrow(dag_edges))
  if (!"effect" %in% names(dag_edges)) dag_edges$effect <- rep("unknown", nrow(dag_edges))
  # Canonical ordering makes the selected representative reproducible on ties.
  dag_edges <- dag_edges[order(dag_edges$tf, dag_edges$target,
                               -dag_edges$confidence_score), , drop = FALSE]
  dag_edges <- dag_edges[!duplicated(dag_edges[, c("tf", "target")]), , drop = FALSE]
  nodes <- sort(unique(c(dag_edges$tf, dag_edges$target, exposure, outcome)))
  graph <- igraph::graph_from_data_frame(dag_edges[, c("tf", "target")],
    directed = TRUE, vertices = data.frame(name = nodes))
  if (!igraph::is_dag(graph)) stop("`dag_edges` must define a DAG.", call. = FALSE)
  forbidden <- setdiff(igraph::as_ids(igraph::subcomponent(graph, exposure, mode = "out")), exposure)
  backdoor <- igraph::delete_edges(graph, igraph::E(graph)[
    igraph::ends(graph, igraph::E(graph), names = TRUE)[, 1L] == exposure])
  ancestors <- sort(unique(c(
    igraph::as_ids(igraph::subcomponent(backdoor, exposure, mode = "in")),
    igraph::as_ids(igraph::subcomponent(backdoor, outcome, mode = "in")))))
  h <- igraph::induced_subgraph(backdoor, ancestors)
  h_nodes <- igraph::V(h)$name
  allowed <- setdiff(h_nodes, c(exposure, outcome, forbidden))
  if (!is.null(eligible_nodes)) {
    eligible_nodes <- as.character(eligible_nodes)
    if (anyNA(eligible_nodes) || any(!eligible_nodes %in% nodes)) {
      stop("`eligible_nodes` must contain existing, nonmissing graph nodes.", call. = FALSE)
    }
    allowed <- intersect(allowed, eligible_nodes)
  }
  # A[u,v] means u -> v. A A' connects co-parents of every retained child.
  adjacency <- igraph::as_adjacency_matrix(h, sparse = TRUE)
  moral <- adjacency + Matrix::t(adjacency) + Matrix::tcrossprod(adjacency)
  entries <- Matrix::summary(moral)
  undirected <- entries[entries$i < entries$j & entries$x > 0, c("i", "j"), drop = FALSE]
  undirected <- undirected[order(undirected$i, undirected$j), , drop = FALSE]
  n <- length(h_nodes)
  inside <- seq_len(n)
  outside <- n + inside
  big <- length(allowed) + 1
  internal_capacity <- ifelse(h_nodes %in% allowed, 1, big)
  split_edges <- rbind(cbind(inside, outside),
    cbind(n + undirected$i, undirected$j),
    cbind(n + undirected$j, undirected$i))
  flow_graph <- igraph::make_empty_graph(2L * n, directed = TRUE)
  flow_graph <- igraph::add_edges(flow_graph, as.vector(t(split_edges)))
  capacity <- c(internal_capacity, rep(big, 2L * nrow(undirected)))
  cut <- igraph::min_cut(flow_graph,
    source = n + match(exposure, h_nodes), target = match(outcome, h_nodes),
    capacity = capacity, value.only = FALSE)
  if (cut$value >= big) {
    stop("No valid back-door adjustment set exists among the allowed nodes.", call. = FALSE)
  }
  cut_ids <- as.integer(cut$cut)
  if (any(cut_ids > n) || any(internal_capacity[cut_ids] != 1)) {
    stop("Internal error: minimum cut contains a protected edge.", call. = FALSE)
  }
  final_set <- sort(h_nodes[cut_ids])
  checker <- tf_adjustment_build_checker(dag_edges, exposure, outcome)
  if (abs(cut$value - length(final_set)) > 1e-8 || !checker$is_valid(final_set)) {
    stop("Internal error: minimum-cut result failed back-door validation.", call. = FALSE)
  }
  # C^0 is descriptive only; it does not restrict the exact optimizer.
  direct_parents <- sort(setdiff(dag_edges$tf[dag_edges$target == outcome],
                                 c(exposure, outcome, forbidden)))
  exposure_parents <- sort(setdiff(dag_edges$tf[dag_edges$target == exposure],
                                   c(exposure, outcome, forbidden)))
  # Multi-source traversal from ancestors(T), avoiding an all-pairs matrix.
  seen <- rep(FALSE, length(nodes))
  queue <- match(igraph::as_ids(igraph::subcomponent(graph, exposure, mode = "in")), nodes)
  seen[queue] <- TRUE
  children <- igraph::as_adj_list(graph, mode = "out")
  cursor <- 1L
  while (cursor <= length(queue)) {
    next_nodes <- as.integer(children[[queue[[cursor]]]])
    next_nodes <- next_nodes[!seen[next_nodes]]
    seen[next_nodes] <- TRUE
    queue <- c(queue, next_nodes)
    cursor <- cursor + 1L
  }
  initial_set <- intersect(direct_parents, nodes[seen])
  score <- tf_adjustment_score_paths(dag_edges, final_set, outcome, exposure)
  elapsed <- proc.time()[["elapsed"]] - started
  operations <- data.frame(iteration = integer(), action = character(),
    removed = character(), added = character(), size_before = integer(), size_after = integer())
  winner <- list(start_id = 1L, seed = NA_integer_, final_set = final_set,
    lineage = stats::setNames(lapply(final_set, identity), final_set),
    operations = operations, pair_checks = 0L, ancestor_checks = 0L,
    validation_checks = 1L, iterations = 1L, terminated_by_limit = FALSE,
    final_valid = TRUE, elapsed_seconds = elapsed,
    backdoor_edges = checker$backdoor_edge_count)
  # One row retained under legacy names for downstream/API compatibility.
  summary <- data.frame(start_id = 1L, seed = NA_integer_,
    final_adjustment_size = length(final_set),
    path_edge_mean_confidence = score$path_edge_mean_confidence,
    total_path_confidence = score$total_path_confidence,
    total_path_edges = score$total_path_edges, final_valid = TRUE, operations = 0L,
    pair_checks = 0L, ancestor_checks = 0L, validation_checks = 1L,
    elapsed_seconds = elapsed, terminated_by_limit = FALSE)
  list(winner = winner, winner_score = score, all_starts = list(winner),
    multistart_summary = summary, initial_set = initial_set,
    direct_parents = direct_parents, exposure_parents = exposure_parents,
    n_starts = 1L, cores = 1L, seed = NA_integer_, search_seconds = elapsed,
    algorithm = "minimum_vertex_cut_v1", minimum_cardinality = length(final_set),
    optimality_certified = TRUE, cut_value = cut$value,
    allowed_nodes = sort(allowed), forbidden_nodes = sort(forbidden),
    ancestor_nodes = n, moral_edges = nrow(undirected),
    confidence_optimized = FALSE,
    path_rule = "shortest directed back-door-graph path; highest confidence among equally short paths",
    optimization_scope = "minimum cardinality under Pearl back-door criterion and eligible_nodes restriction")
}

# Compatibility entry point: restart/seed/core/iteration options no longer
# affect this deterministic single-query optimization. Outer gene parallelism
# is handled by query_TF_target_adjustment_sets().
tf_adjustment_search <- function(
  dag_edges, exposure, outcome, n_starts = 1L, cores = 1L, seed = 123L,
  eligible_nodes = NULL, max_iterations = 1000L, max_pair_checks = 200000L,
  max_ancestor_candidates = 500L, max_pairs_per_ancestor = 20L
) {
  result <- tf_minimum_adjustment_search(dag_edges, exposure, outcome, eligible_nodes)
  result$legacy_options <- list(n_starts = n_starts, cores = cores, seed = seed,
    max_iterations = max_iterations, max_pair_checks = max_pair_checks,
    max_ancestor_candidates = max_ancestor_candidates,
    max_pairs_per_ancestor = max_pairs_per_ancestor)
  result
}

# Legacy internal entry point; delegates to the exact solver, not a heuristic.
tf_recursive_adjustment_search <- tf_adjustment_search
