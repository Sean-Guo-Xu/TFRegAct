#!/usr/bin/env Rscript

find_adjustment_dagitty_get_arg <- function(args, flag, default = NULL) {
  hit <- grep(paste0("^", flag, "="), args, value = TRUE)
  if (length(hit) == 0) {
    return(default)
  }
  sub(paste0("^", flag, "="), "", hit[[1]])
}

find_adjustment_dagitty_install_if_missing <- function(pkg) {
  .tfregact_require(pkg)
}

find_adjustment_dagitty_load_bundle_file <- function(path) {
  env <- new.env(parent = emptyenv())
  load(path, envir = env)
  if (!exists("dagitty_network_bundle", envir = env, inherits = FALSE)) {
    stop(sprintf("Bundle file does not contain `dagitty_network_bundle`: %s", path))
  }
  get("dagitty_network_bundle", envir = env, inherits = FALSE)
}

find_adjustment_dagitty_extract_bundle <- function(network = NULL, edge_file = NULL) {
  if (is.null(edge_file) || !length(edge_file) || !nzchar(as.character(edge_file[[1]]))) {
    edge_file <- getOption("TFRegAct.network_path", Sys.getenv("TFREGACT_NETWORK_PATH", unset = ""))
  }
  if (!is.null(network)) {
    if (is.list(network) && !is.null(network$dagitty_bundle)) {
      return(network$dagitty_bundle)
    }

    if (is.list(network) && !is.null(network$files) && !is.null(network$files$dagitty_rdata)) {
      bundle_path <- as.character(network$files$dagitty_rdata[[1]])
      if (nzchar(bundle_path) && file.exists(bundle_path)) {
        return(find_adjustment_dagitty_load_bundle_file(bundle_path))
      }
    }

    if (is.character(network) && length(network) == 1 && grepl("\\.RData$", network, ignore.case = TRUE) && file.exists(network)) {
      return(find_adjustment_dagitty_load_bundle_file(network))
    }
  }

  if (!is.null(edge_file) && nzchar(edge_file) && grepl("\\.RData$", edge_file, ignore.case = TRUE) && file.exists(edge_file)) {
    return(find_adjustment_dagitty_load_bundle_file(edge_file))
  }

  default_bundle <- "tf_union_output/TF_Full_Map.RData"
  if (file.exists(default_bundle)) {
    return(find_adjustment_dagitty_load_bundle_file(default_bundle))
  }

  NULL
}

find_adjustment_dagitty_extract_edge_file <- function(network = NULL, edge_file = NULL) {
  if (is.null(edge_file) || !length(edge_file) || !nzchar(as.character(edge_file[[1]]))) {
    edge_file <- getOption("TFRegAct.network_path", Sys.getenv("TFREGACT_NETWORK_PATH", unset = ""))
  }
  if (!is.null(network)) {
    if (is.character(network) && length(network) == 1 && nzchar(network) && file.exists(network)) {
      if (grepl("\\.RData$", network, ignore.case = TRUE)) {
        return(NULL)
      }
      return(normalizePath(network, winslash = "/", mustWork = TRUE))
    }

    if (is.list(network) && !is.null(network$files) && !is.null(network$files$csv)) {
      csv_path <- as.character(network$files$csv[[1]])
      if (nzchar(csv_path) && file.exists(csv_path)) {
        return(normalizePath(csv_path, winslash = "/", mustWork = TRUE))
      }
    }

    stop("`network` must be a valid edge-file path or a build_tf_network() result with `files$csv`.")
  }

  if (!is.null(edge_file) && nzchar(edge_file) && file.exists(edge_file)) {
    if (grepl("\\.RData$", edge_file, ignore.case = TRUE)) {
      return(NULL)
    }
    return(normalizePath(edge_file, winslash = "/", mustWork = TRUE))
  }

  default_file <- "tf_union_output/tf_gene_merged_human_weighted_clean.csv"
  if (file.exists(default_file)) {
    return(normalizePath(default_file, winslash = "/", mustWork = TRUE))
  }

  if (exists("build_tf_network", mode = "function")) {
    build_res <- build_tf_network()
    if (is.list(build_res) && !is.null(build_res$files) && !is.null(build_res$files$csv)) {
      csv_path <- as.character(build_res$files$csv[[1]])
      if (nzchar(csv_path) && file.exists(csv_path)) {
        return(normalizePath(csv_path, winslash = "/", mustWork = TRUE))
      }
    }
  }

  stop("No usable TF network file was found. Supply `edge_file` explicitly or set options(TFRegAct.network_path = '<path>').")
}

find_adjustment_dagitty_read_edges <- function(path) {
  ext <- tolower(tools::file_ext(path))
  if (ext == "csv") {
    df <- utils::read.csv(path, stringsAsFactors = FALSE, check.names = FALSE)
  } else {
    df <- utils::read.delim(path, stringsAsFactors = FALSE, check.names = FALSE)
  }

  need_cols <- c("tf", "target")
  if (!all(need_cols %in% colnames(df))) {
    stop("Edge file must contain at least columns: tf, target")
  }

  keep_cols <- intersect(c("tf", "target", "effect", "confidence_score", "supporting_databases"), colnames(df))
  df <- df[, keep_cols, drop = FALSE]
  if (!("effect" %in% colnames(df))) {
    df$effect <- "unknown"
  }
  if (!("confidence_score" %in% colnames(df))) {
    df$confidence_score <- NA_integer_
  }
  if (!("supporting_databases" %in% colnames(df))) {
    df$supporting_databases <- ""
  }

  df$tf <- trimws(as.character(df$tf))
  df$target <- trimws(as.character(df$target))
  df$effect <- trimws(as.character(df$effect))
  df$confidence_score <- suppressWarnings(as.integer(df$confidence_score))
  df$supporting_databases <- trimws(as.character(df$supporting_databases))

  if ("supporting_databases" %in% colnames(df)) {
    db_tokens <- strsplit(df$supporting_databases, "\\s*;\\s*")
    grouped_score <- vapply(db_tokens, function(x) {
      vals <- trimws(x)
      vals <- vals[vals != ""]
      has_trrust <- any(vals == "TRRUST")
      has_regnetwork <- any(vals == "RegNetwork")
      has_chip <- any(vals %in% c("ChEA", "ENCODE"))
      has_low <- any(vals %in% c("JASPAR"))
      4L * as.integer(has_trrust) +
        3L * as.integer(has_regnetwork) +
        2L * as.integer(has_chip) +
        1L * as.integer(has_low)
    }, integer(1))
    df$confidence_score <- as.integer(grouped_score)
  }

  df <- df[!is.na(df$tf) & !is.na(df$target) & df$tf != "" & df$target != "", , drop = FALSE]
  df <- df[df$tf != df$target, , drop = FALSE]
  unique(df)
}

find_adjustment_dagitty_filter_prior_edges <- function(
  edges,
  confidence_threshold,
  exposure = NULL,
  outcome = NULL
) {
  confidence_threshold <- suppressWarnings(as.numeric(confidence_threshold[[1L]]))
  if (!is.finite(confidence_threshold) || confidence_threshold < 0) {
    stop("`confidence_threshold` must be finite and non-negative.", call. = FALSE)
  }
  confidence <- suppressWarnings(as.numeric(edges$confidence_score))
  keep <- is.finite(confidence) & confidence >= confidence_threshold

  # The queried exposure -> outcome edge defines the estimand and may have been
  # manually accepted upstream. Keep it even when its database score is below
  # the confounder-network threshold; it is removed from the back-door graph
  # during validation and therefore cannot create a spurious back-door path.
  if (!is.null(exposure) && !is.null(outcome)) {
    keep_query <- toupper(trimws(as.character(edges$tf))) ==
      toupper(trimws(as.character(exposure[[1L]]))) &
      toupper(trimws(as.character(edges$target))) ==
      toupper(trimws(as.character(outcome[[1L]])))
    keep_query[is.na(keep_query)] <- FALSE
    keep <- keep | keep_query
  }
  filtered <- edges[keep, , drop = FALSE]
  attr(filtered, "confidence_threshold") <- confidence_threshold
  attr(filtered, "edges_before_threshold") <- nrow(edges)
  attr(filtered, "edges_after_threshold") <- nrow(filtered)
  filtered
}

find_adjustment_dagitty_sanitize_node <- function(x) {
  y <- gsub("[^A-Za-z0-9_]", "_", x)
  y <- gsub("_+", "_", y)
  y <- sub("^_+", "", y)
  ifelse(grepl("^[A-Za-z]", y), y, paste0("N_", y))
}

find_adjustment_dagitty_make_node_map <- function(nodes) {
  safe <- find_adjustment_dagitty_sanitize_node(nodes)
  out <- character(length(nodes))
  seen <- integer(0)
  names(seen) <- character(0)

  for (i in seq_along(nodes)) {
    base <- safe[[i]]
    if (!base %in% names(seen)) {
      seen[[base]] <- 1L
      out[[i]] <- base
    } else {
      seen[[base]] <- seen[[base]] + 1L
      out[[i]] <- sprintf("%s_%d", base, seen[[base]])
    }
  }

  data.frame(node = nodes, dagitty_id = out, stringsAsFactors = FALSE)
}

find_adjustment_dagitty_direct_distance_map <- function(graph, target_node) {
  if (!target_node %in% igraph::V(graph)$name) {
    stop(sprintf("Target node '%s' not found in graph.", target_node))
  }

  rev_graph <- igraph::reverse_edges(graph)
  bfs_res <- igraph::bfs(
    rev_graph,
    root = target_node,
    mode = "out",
    unreachable = FALSE,
    dist = TRUE,
    order = FALSE
  )

  dist_vals <- as.numeric(bfs_res$dist)
  dist_vals[dist_vals < 0] <- Inf
  names(dist_vals) <- igraph::V(graph)$name
  dist_vals
}

find_adjustment_dagitty_build_local_dag <- function(edges, tf_query, gene_query, graph = NULL, reverse_graph = NULL) {
  g <- if (is.null(graph)) {
    igraph::graph_from_data_frame(edges[, c("tf", "target"), drop = FALSE], directed = TRUE)
  } else {
    graph
  }
  rev_g <- if (is.null(reverse_graph)) {
    igraph::reverse_edges(g)
  } else {
    reverse_graph
  }

  anc_gene <- igraph::subcomponent(rev_g, gene_query, mode = "out")
  anc_tf <- igraph::subcomponent(rev_g, tf_query, mode = "out")

  # For backdoor adjustment we only need the ancestor graph of exposure/outcome.
  relevant_nodes <- unique(c(igraph::V(g)$name[anc_gene], igraph::V(g)$name[anc_tf], tf_query, gene_query))
  sub_g <- igraph::induced_subgraph(g, vids = relevant_nodes)
  sub_edges <- igraph::as_data_frame(sub_g, what = "edges")
  sub_edges <- merge(sub_edges, edges, by.x = c("from", "to"), by.y = c("tf", "target"), all.x = TRUE, sort = FALSE)
  colnames(sub_edges)[colnames(sub_edges) == "from"] <- "tf"
  colnames(sub_edges)[colnames(sub_edges) == "to"] <- "target"

  nodes <- sort(unique(c(sub_edges$tf, sub_edges$target, tf_query, gene_query)))
  dist_map_gene <- find_adjustment_dagitty_direct_distance_map(g, gene_query)
  dist_map_tf <- find_adjustment_dagitty_direct_distance_map(g, tf_query)
  d_to_gene <- as.numeric(dist_map_gene[nodes])
  d_to_tf <- as.numeric(dist_map_tf[nodes])

  node_meta <- data.frame(
    node = nodes,
    d_to_gene = d_to_gene,
    d_to_tf = d_to_tf,
    stringsAsFactors = FALSE
  )

  node_meta$class_rank <- 2L
  node_meta$class_rank[is.finite(node_meta$d_to_tf) & is.finite(node_meta$d_to_gene) & node_meta$node != tf_query & node_meta$node != gene_query] <- 0L
  node_meta$class_rank[node_meta$node == tf_query] <- 1L
  node_meta$class_rank[node_meta$node == gene_query] <- 3L

  node_meta$order_score <- with(
    node_meta,
    class_rank * 1e6 +
      ifelse(is.finite(d_to_tf), d_to_tf, 999) * 1e4 +
      ifelse(is.finite(d_to_gene), d_to_gene, 999)
  )

  order_map <- setNames(rank(node_meta$order_score, ties.method = "first"), node_meta$node)
  sub_edges$from_rank <- order_map[sub_edges$tf]
  sub_edges$to_rank <- order_map[sub_edges$target]
  dag_edges <- sub_edges[sub_edges$from_rank < sub_edges$to_rank, , drop = FALSE]
  dag_edges <- dag_edges[!duplicated(dag_edges[, c("tf", "target")]), , drop = FALSE]

  list(
    dag_edges = dag_edges,
    node_meta = node_meta
  )
}

find_adjustment_dagitty_choose_recommended_set <- function(adj_sets_df, dag_edges_out, tf_query, gene_query) {
  if (nrow(adj_sets_df) == 0) {
    return(list(
      recommended_set_id = NA_integer_,
      recommended_variables = character(0),
      recommended_table = data.frame(
        set_id = integer(0),
        variable = character(0),
        variable_confidence = integer(0),
        stringsAsFactors = FALSE
      ),
      set_scores = data.frame(
        set_id = integer(0),
        set_size = integer(0),
        total_confidence = integer(0),
        stringsAsFactors = FALSE
      )
    ))
  }

  vars <- sort(unique(adj_sets_df$variable))
  variable_conf <- vapply(vars, function(v) {
    vals <- dag_edges_out$confidence_score[dag_edges_out$tf == v | dag_edges_out$target == v]
    vals <- vals[is.finite(vals)]
    if (length(vals) == 0) {
      0L
    } else {
      as.integer(max(vals))
    }
  }, integer(1))
  var_conf_df <- data.frame(
    variable = vars,
    variable_confidence = as.integer(variable_conf),
    stringsAsFactors = FALSE
  )

  scored_df <- merge(adj_sets_df, var_conf_df, by = "variable", all.x = TRUE, sort = FALSE)
  scored_df$variable_confidence[is.na(scored_df$variable_confidence)] <- 0L

  split_sets <- split(scored_df, scored_df$set_id)
  set_scores <- do.call(
    rbind,
    lapply(split_sets, function(chunk) {
      data.frame(
        set_id = chunk$set_id[[1]],
        set_size = nrow(chunk),
        total_confidence = sum(chunk$variable_confidence),
        stringsAsFactors = FALSE
      )
    })
  )
  set_scores <- set_scores[order(set_scores$set_size, -set_scores$total_confidence, set_scores$set_id), , drop = FALSE]
  best_id <- set_scores$set_id[[1]]
  recommended_table <- scored_df[scored_df$set_id == best_id, c("set_id", "variable", "variable_confidence"), drop = FALSE]
  recommended_table <- recommended_table[order(-recommended_table$variable_confidence, recommended_table$variable), , drop = FALSE]
  rownames(recommended_table) <- NULL

  list(
    recommended_set_id = best_id,
    recommended_variables = recommended_table$variable,
    recommended_table = recommended_table,
    set_scores = set_scores
  )
}

find_adjustment_dagitty_prepare_target_tree <- function(graph, target_node) {
  if (!target_node %in% igraph::V(graph)$name) {
    stop(sprintf("Target node '%s' not found in graph.", target_node))
  }

  bfs_res <- igraph::bfs(
    graph,
    root = target_node,
    mode = "all",
    unreachable = TRUE,
    dist = TRUE,
    parent = TRUE,
    order = FALSE
  )

  vertex_names <- igraph::V(graph)$name
  dist_vals <- bfs_res$dist
  names(dist_vals) <- vertex_names

  father_idx <- bfs_res$parent
  father_names <- rep(NA_character_, length(father_idx))
  valid_father <- !is.na(father_idx) & father_idx > 0
  father_names[valid_father] <- vertex_names[father_idx[valid_father]]
  names(father_names) <- vertex_names

  list(
    graph = graph,
    target = target_node,
    distances = dist_vals,
    fathers = father_names
  )
}

find_adjustment_dagitty_best_path_metric <- function(tree_info, from_node) {
  graph <- tree_info$graph
  target_node <- tree_info$target

  if (!from_node %in% names(tree_info$distances)) {
    return(list(distance = Inf, avg_confidence = 0))
  }

  if (identical(from_node, target_node)) {
    return(list(distance = 0, avg_confidence = 0))
  }

  dist_val <- tree_info$distances[[from_node]]
  if (!is.finite(dist_val)) {
    return(list(distance = Inf, avg_confidence = 0))
  }

  total_conf <- 0
  step_count <- 0
  current <- from_node

  while (!identical(current, target_node)) {
    parent <- tree_info$fathers[[current]]
    if (is.na(parent) || !nzchar(parent)) {
      return(list(distance = Inf, avg_confidence = 0))
    }

    edge_id <- igraph::get_edge_ids(graph, vp = c(current, parent), directed = FALSE)
    if (length(edge_id) == 0 || edge_id[[1]] <= 0) {
      return(list(distance = Inf, avg_confidence = 0))
    }

    conf_val <- igraph::E(graph)$confidence_score[[edge_id[[1]]]]
    if (is.finite(conf_val)) {
      total_conf <- total_conf + conf_val
    }

    step_count <- step_count + 1
    current <- parent
  }

  avg_conf <- if (step_count == 0) 0 else total_conf / step_count
  list(distance = step_count, avg_confidence = avg_conf)
}

find_adjustment_dagitty_build_distance_graph <- function(edges) {
  edge_dt <- data.table::as.data.table(edges[, c("tf", "target", "confidence_score"), drop = FALSE])
  edge_dt[, node1 := ifelse(tf <= target, tf, target)]
  edge_dt[, node2 := ifelse(tf <= target, target, tf)]
  collapsed <- edge_dt[
    ,
    .(confidence_score = as.integer(max(confidence_score, na.rm = TRUE))),
    by = .(node1, node2)
  ]
  collapsed[!is.finite(confidence_score), confidence_score := 0L]
  graph <- igraph::graph_from_data_frame(
    collapsed[, .(from = node1, to = node2, confidence_score)],
    directed = FALSE
  )
  igraph::E(graph)$confidence_score <- collapsed$confidence_score
  graph
}

find_adjustment_dagitty_direct_effect_to_gene <- function(edges, gene_query, regulators) {
  regulator_vec <- unique(as.character(regulators))
  out <- rep("unknown", length(regulator_vec))
  names(out) <- regulator_vec

  edge_sub <- edges[edges$target == gene_query & edges$tf %in% regulator_vec, c("tf", "effect"), drop = FALSE]
  if (nrow(edge_sub) == 0) {
    return(out)
  }

  split_effects <- split(trimws(as.character(edge_sub$effect)), edge_sub$tf)
  for (reg in names(split_effects)) {
    vals <- split_effects[[reg]]
    vals <- vals[!is.na(vals) & vals != ""]
    vals <- vals[vals %in% c("activation", "repression", "unknown", "mixed")]
    vals_known <- vals[vals %in% c("activation", "repression")]

    out[[reg]] <- if (length(vals_known) == 0) {
      "unknown"
    } else if (all(vals_known == "activation")) {
      "activation"
    } else if (all(vals_known == "repression")) {
      "repression"
    } else {
      "unknown"
    }
  }

  out
}

find_adjustment_dagitty_effect_to_direction <- function(effect) {
  effect <- tolower(trimws(as.character(effect)))
  out <- integer(length(effect))
  out[effect %in% c("activation", "activate", "positive", "up")] <- 1L
  out[effect %in% c("repression", "repress", "negative", "down")] <- -1L
  out
}

find_adjustment_dagitty_direct_edge_metrics <- function(
  edges,
  gene_query,
  regulators
) {
  regulators <- unique(trimws(as.character(regulators)))
  out <- data.frame(
    tf = regulators,
    distance_to_gene_A = rep(Inf, length(regulators)),
    effect_on_gene_A = rep("unknown", length(regulators)),
    direction_to_gene_A = integer(length(regulators)),
    avg_confidence_to_gene_A = numeric(length(regulators)),
    stringsAsFactors = FALSE
  )
  if (length(regulators) == 0L) {
    return(out)
  }

  direct <- edges[
    edges$target == gene_query & edges$tf %in% regulators,
    c("tf", "effect", "confidence_score"),
    drop = FALSE
  ]
  if (nrow(direct) == 0L) {
    return(out)
  }

  grouped <- split(seq_len(nrow(direct)), direct$tf)
  for (regulator in names(grouped)) {
    block <- direct[grouped[[regulator]], , drop = FALSE]
    directions <- unique(find_adjustment_dagitty_effect_to_direction(block$effect))
    known_directions <- unique(directions[directions != 0L])
    direction <- if (length(known_directions) == 1L) {
      known_directions[[1]]
    } else {
      0L
    }
    effect <- if (direction == 1L) {
      "activation"
    } else if (direction == -1L) {
      "repression"
    } else if (length(known_directions) > 1L) {
      "mixed"
    } else {
      "unknown"
    }
    confidence <- suppressWarnings(as.numeric(block$confidence_score))
    confidence <- confidence[is.finite(confidence)]
    row_index <- match(regulator, out$tf)
    out$distance_to_gene_A[[row_index]] <- 1
    out$effect_on_gene_A[[row_index]] <- effect
    out$direction_to_gene_A[[row_index]] <- as.integer(direction)
    out$avg_confidence_to_gene_A[[row_index]] <- if (length(confidence)) {
      max(confidence)
    } else {
      0
    }
  }
  out
}

find_adjustment_dagitty_build_recommended_metrics <- function(
  edges,
  tf_query,
  gene_query,
  recommended_variables,
  beta = 2,
  direct_edges = edges
) {
  beta <- suppressWarnings(as.numeric(beta[[1]]))
  if (!is.finite(beta) || beta < 0) {
    stop("`beta` must be a non-negative number.")
  }

  graph <- find_adjustment_dagitty_build_distance_graph(edges)
  tf_tree <- find_adjustment_dagitty_prepare_target_tree(graph, tf_query)

  vars <- unique(c(recommended_variables, tf_query, gene_query))
  backdoor_vars <- setdiff(vars, tf_query)
  backdoor_paths <- tf_adjustment_score_paths(
    dag_edges = edges,
    adjustment_set = backdoor_vars,
    outcome = gene_query,
    exposure = tf_query
  )$details
  target_path <- tf_adjustment_score_paths(
    dag_edges = edges,
    adjustment_set = tf_query,
    outcome = gene_query
  )$details
  directed_paths <- rbind(backdoor_paths, target_path)
  directed_paths <- directed_paths[
    match(vars, directed_paths$variable), , drop = FALSE
  ]
  fallback_vars <- directed_paths$variable[
    directed_paths$variable %in% recommended_variables & !is.finite(directed_paths$path_length)
  ]
  tf_paths <- if (length(fallback_vars)) tf_adjustment_score_paths(
    dag_edges = edges, adjustment_set = fallback_vars, outcome = tf_query
  )$details else NULL
  out <- lapply(vars, function(v) {
    path_index <- match(v, directed_paths$variable)
    path_distance <- directed_paths$path_length[[path_index]]
    path_confidence <- directed_paths$path_mean_confidence[[path_index]]
    path_direction <- directed_paths$path_direction[[path_index]]
    path_effect <- directed_paths$path_effect[[path_index]]
    to_tf <- find_adjustment_dagitty_best_path_metric(tf_tree, v)
    # For a separator with no directed Y-path, use its directed path to T,
    # never an undirected shortcut through the specially retained query edge.
    if (!is.finite(path_distance) && v %in% recommended_variables) {
      tf_index <- match(v, tf_paths$variable)
      to_tf <- list(distance = tf_paths$path_length[[tf_index]],
                    avg_confidence = tf_paths$path_mean_confidence[[tf_index]])
    }
    avg_to_gene <- if (is.finite(path_distance) && path_distance > 0) {
      path_confidence
    } else {
      0
    }
    avg_to_tf <- if (is.finite(to_tf$distance)) to_tf$avg_confidence else 0
    overall_avg <- if (is.finite(path_distance) && path_distance > 0 && is.finite(to_tf$distance)) {
      (beta * avg_to_gene + avg_to_tf) / (1 + beta)
    } else if (is.finite(path_distance) && path_distance > 0) {
      avg_to_gene
    } else if (is.finite(to_tf$distance)) {
      avg_to_tf
    } else {
      0
    }
    data.frame(
      tf = v,
      distance_to_gene_A = path_distance,
      effect_on_gene_A = path_effect,
      direction_to_gene_A = path_direction,
      avg_confidence_to_gene_A = avg_to_gene,
      selected_path_to_gene_A = directed_paths$selected_path[[path_index]],
      distance_to_tf_B = if (is.finite(to_tf$distance)) to_tf$distance else Inf,
      avg_confidence_to_tf_B = avg_to_tf,
      overall_avg_confidence = overall_avg,
      stringsAsFactors = FALSE
    )
  })
  out <- do.call(rbind, out)
  out <- out[order(out$tf), , drop = FALSE]
  rownames(out) <- NULL
  out
}

find_adjustment_dagitty_recommended_metrics <- function(
  tf,
  gene,
  network = NULL,
  edge_file = "tf_union_output/tf_gene_merged_human_weighted_clean.csv",
  recommended_set_file = NULL,
  outdir = "adjustment_output",
  beta = 2,
  overall_confidence_threshold = 3,
  build_if_missing = TRUE,
  search_starts = 1L,
  search_cores = 4L,
  search_seed = 123L,
  max_adjustment_sets = NULL,
  write_full_outputs = FALSE,
  write_files = TRUE
) {
  tf_query <- trimws(as.character(tf))
  gene_query <- trimws(as.character(gene))
  if (!nzchar(tf_query) || !nzchar(gene_query)) {
    stop("`tf` and `gene` must be non-empty strings.")
  }
  overall_confidence_threshold <- suppressWarnings(as.numeric(overall_confidence_threshold[[1]]))
  if (!is.finite(overall_confidence_threshold)) {
    stop("`overall_confidence_threshold` must be a finite number.")
  }

  local_dag_edge_file <- file.path(outdir, sprintf("dagitty_local_dag_edges_%s_%s.csv", tf_query, gene_query))
  bundle <- find_adjustment_dagitty_extract_bundle(network = network, edge_file = edge_file)
  use_default_edge <- is.null(edge_file) || identical(edge_file, "tf_union_output/tf_gene_merged_human_weighted_clean.csv")
  resolved_edge_file <- if (!is.null(bundle)) {
    if (!is.null(bundle$metadata$source_csv) && nzchar(bundle$metadata$source_csv)) {
      bundle$metadata$source_csv
    } else {
      edge_file
    }
  } else if (use_default_edge && file.exists(local_dag_edge_file)) {
    normalizePath(local_dag_edge_file, winslash = "/", mustWork = TRUE)
  } else {
    find_adjustment_dagitty_extract_edge_file(network = network, edge_file = edge_file)
  }
  source_edges <- if (!is.null(bundle)) bundle$edges else find_adjustment_dagitty_read_edges(resolved_edge_file)
  edges <- find_adjustment_dagitty_filter_prior_edges(
    edges = source_edges,
    confidence_threshold = overall_confidence_threshold,
    exposure = tf_query,
    outcome = gene_query
  )

  if (is.null(recommended_set_file) || !nzchar(recommended_set_file)) {
    recommended_set_file <- file.path(outdir, sprintf("dagitty_recommended_adjustment_set_%s_%s.csv", tf_query, gene_query))
  }
  if (!file.exists(recommended_set_file)) {
    if (!isTRUE(build_if_missing)) {
      stop(sprintf("Recommended adjustment set file not found: %s", recommended_set_file))
    }

    generated <- find_adjustment_dagitty_run(
      tf = tf_query,
      gene = gene_query,
      network = network,
      edge_file = edge_file,
      outdir = outdir,
      beta = beta,
      overall_confidence_threshold = overall_confidence_threshold,
      search_starts = search_starts,
      search_cores = search_cores,
      search_seed = search_seed,
      max_adjustment_sets = max_adjustment_sets,
      write_full_outputs = write_full_outputs,
      write_files = write_files
    )

    return(list(
      tf = tf_query,
      gene = gene_query,
      metrics = generated$recommended_adjustment_metrics,
      filtered_metrics = generated$filtered_recommended_adjustment_metrics,
      overall_confidence_threshold = overall_confidence_threshold,
      files = generated$files,
      generated_recommended_set = TRUE
    ))
  }

  rec_df <- utils::read.csv(recommended_set_file, stringsAsFactors = FALSE, check.names = FALSE)
  if (!("variable" %in% colnames(rec_df))) {
    stop("Recommended adjustment set file must contain a 'variable' column.")
  }

  recommended_variables <- unique(trimws(as.character(rec_df$variable)))
  recommended_variables <- recommended_variables[!is.na(recommended_variables) & recommended_variables != ""]
  metrics <- find_adjustment_dagitty_build_recommended_metrics(
    edges = edges,
    tf_query = tf_query,
    gene_query = gene_query,
    recommended_variables = recommended_variables,
    beta = beta,
    direct_edges = source_edges
  )
  keep_forced <- metrics$tf %in% c(tf_query, gene_query)
  keep_adjustment_path <- metrics$tf %in% recommended_variables
  filtered_metrics <- metrics[keep_forced | keep_adjustment_path, , drop = FALSE]
  filtered_metrics$meets_prior_edge_threshold <-
    keep_forced[keep_forced | keep_adjustment_path] |
    ifelse(is.finite(filtered_metrics$distance_to_gene_A) &
             filtered_metrics$distance_to_gene_A > 0,
           filtered_metrics$avg_confidence_to_gene_A,
           filtered_metrics$avg_confidence_to_tf_B) >= overall_confidence_threshold
  if (any(!filtered_metrics$meets_prior_edge_threshold)) {
    stop(
      "Saved adjustment set is incompatible with the requested prior-edge confidence threshold.",
      call. = FALSE
    )
  }
  rownames(filtered_metrics) <- NULL

  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }
  metrics_file <- file.path(outdir, sprintf("dagitty_recommended_adjustment_metrics_%s_%s.csv", tf_query, gene_query))
  filtered_metrics_file <- file.path(
    outdir,
    sprintf("dagitty_recommended_adjustment_metrics_filtered_%s_%s.csv", tf_query, gene_query)
  )
  if (write_files) {
    utils::write.csv(filtered_metrics, file = filtered_metrics_file, row.names = FALSE)
    if (isTRUE(write_full_outputs)) {
      utils::write.csv(metrics, file = metrics_file, row.names = FALSE)
    }
  }

  list(
    tf = tf_query,
    gene = gene_query,
    metrics = metrics,
    filtered_metrics = filtered_metrics,
    overall_confidence_threshold = overall_confidence_threshold,
    files = list(
      recommended_set = normalizePath(recommended_set_file, winslash = "/", mustWork = TRUE),
      recommended_adjustment_metrics = normalizePath(metrics_file, winslash = "/", mustWork = FALSE),
      recommended_adjustment_metrics_filtered = normalizePath(filtered_metrics_file, winslash = "/", mustWork = FALSE),
      edge_file_used = resolved_edge_file
    )
  )
}

find_adjustment_dagitty_run <- function(
  tf,
  gene,
  network = NULL,
  edge_file = "tf_union_output/tf_gene_merged_human_weighted_clean.csv",
  outdir = "adjustment_output",
  beta = 2,
  overall_confidence_threshold = 3,
  search_starts = 1L,
  search_cores = 4L,
  search_seed = 123L,
  max_iterations = 1000L,
  max_pair_checks = 200000L,
  max_ancestor_candidates = 500L,
  max_pairs_per_ancestor = 20L,
  max_adjustment_sets = NULL,
  include_candidate_confounders = TRUE,
  write_full_outputs = FALSE,
  write_files = TRUE
) {
  find_adjustment_dagitty_install_if_missing("igraph")
  tf_query <- trimws(as.character(tf))
  gene_query <- trimws(as.character(gene))
  if (!nzchar(tf_query) || !nzchar(gene_query)) {
    stop("`tf` and `gene` must be non-empty strings.")
  }
  search_starts <- suppressWarnings(as.integer(search_starts[[1L]]))
  search_cores <- suppressWarnings(as.integer(search_cores[[1L]]))
  search_seed <- suppressWarnings(as.integer(search_seed[[1L]]))
  if (anyNA(c(search_starts, search_cores, search_seed)) ||
      search_starts < 1L || search_cores < 1L) {
    stop(
      "`search_starts` and `search_cores` must be positive integers; `search_seed` must be an integer.",
      call. = FALSE
    )
  }
  if (!is.null(max_adjustment_sets)) {
    warning(
      "`max_adjustment_sets` is ignored; minimum-cut search returns one minimum-cardinality set.",
      call. = FALSE
    )
  }

  bundle <- find_adjustment_dagitty_extract_bundle(network = network, edge_file = edge_file)
  edge_file <- if (is.null(bundle)) {
    find_adjustment_dagitty_extract_edge_file(network = network, edge_file = edge_file)
  } else if (!is.null(bundle$metadata$source_csv) && nzchar(bundle$metadata$source_csv)) {
    bundle$metadata$source_csv
  } else {
    ""
  }
  outdir <- normalizePath(outdir, winslash = "/", mustWork = FALSE)
  if (!dir.exists(outdir)) {
    dir.create(outdir, recursive = TRUE, showWarnings = FALSE)
  }

  source_edges <- if (!is.null(bundle)) bundle$edges else find_adjustment_dagitty_read_edges(edge_file)
  all_nodes <- sort(unique(c(source_edges$tf, source_edges$target)))

  if (!tf_query %in% all_nodes) {
    stop(sprintf("TF '%s' not found in edge file.", tf_query))
  }
  if (!gene_query %in% all_nodes) {
    stop(sprintf("Gene '%s' not found in edge file.", gene_query))
  }
  edges <- find_adjustment_dagitty_filter_prior_edges(
    edges = source_edges,
    confidence_threshold = overall_confidence_threshold,
    exposure = tf_query,
    outcome = gene_query
  )
  edges_before_threshold <- attr(edges, "edges_before_threshold")
  edges_after_threshold <- attr(edges, "edges_after_threshold")
  trusted_nodes <- unique(c(edges$tf, edges$target))
  if (!(tf_query %in% trusted_nodes) || !(gene_query %in% trusted_nodes)) {
    stop(
      sprintf(
        "TF `%s` and gene `%s` are not both present after applying confidence threshold %s.",
        tf_query, gene_query, format(overall_confidence_threshold, trim = TRUE)
      ),
      call. = FALSE
    )
  }

  local_res <- find_adjustment_dagitty_build_local_dag(
    edges,
    tf_query,
    gene_query,
    graph = NULL,
    reverse_graph = NULL
  )
  dag_edges <- local_res$dag_edges
  node_meta <- local_res$node_meta

  local_nodes <- sort(unique(c(dag_edges$tf, dag_edges$target, tf_query, gene_query)))
  node_map <- find_adjustment_dagitty_make_node_map(local_nodes)
  id_map <- setNames(node_map$dagitty_id, node_map$node)

  dagitty_lines <- c("dag {")
  for (id in node_map$dagitty_id) {
    tags <- character(0)
    if (id == id_map[[tf_query]]) {
      tags <- c(tags, "exposure")
    }
    if (id == id_map[[gene_query]]) {
      tags <- c(tags, "outcome")
    }
    if (length(tags) > 0) {
      dagitty_lines <- c(dagitty_lines, sprintf("  %s [%s]", id, paste(tags, collapse = ",")))
    }
  }
  if (nrow(dag_edges) > 0) {
    edge_strings <- sprintf("  %s -> %s", id_map[dag_edges$tf], id_map[dag_edges$target])
    dagitty_lines <- c(dagitty_lines, edge_strings)
  }
  dagitty_lines <- c(dagitty_lines, "}")
  dagitty_string <- paste(dagitty_lines, collapse = "\n")

  dag_edges_out <- dag_edges[, c("tf", "target", "effect", "confidence_score", "supporting_databases"), drop = FALSE]
  dag_edges_out$direction <- find_adjustment_dagitty_effect_to_direction(
    dag_edges_out$effect
  )
  recursive_search <- tf_adjustment_search(
    dag_edges = dag_edges_out,
    exposure = tf_query,
    outcome = gene_query,
    n_starts = search_starts,
    cores = search_cores,
    seed = search_seed,
    max_iterations = max_iterations,
    max_pair_checks = max_pair_checks,
    max_ancestor_candidates = max_ancestor_candidates,
    max_pairs_per_ancestor = max_pairs_per_ancestor
  )
  adj_sets_df <- do.call(rbind, lapply(recursive_search$all_starts, function(item) {
    data.frame(
      set_id = rep(item$start_id, length(item$final_set)),
      variable = item$final_set,
      stringsAsFactors = FALSE
    )
  }))
  candidate_df <- if (isTRUE(include_candidate_confounders)) {
    data.frame(
      variable = recursive_search$initial_set,
      stringsAsFactors = FALSE
    )
  } else {
    data.frame(variable = character(0), stringsAsFactors = FALSE)
  }
  winner <- recursive_search$winner
  winner_paths <- recursive_search$winner_score$details
  winner_paths <- winner_paths[
    match(winner$final_set, winner_paths$variable), , drop = FALSE
  ]
  recommended_table <- data.frame(
    set_id = rep(winner$start_id, length(winner$final_set)),
    variable = winner$final_set,
    variable_confidence = winner_paths$path_mean_confidence,
    path_length = winner_paths$path_length,
    path_confidence_sum = winner_paths$path_confidence_sum,
    selected_path = winner_paths$selected_path,
    stringsAsFactors = FALSE
  )
  set_scores <- data.frame(
    set_id = recursive_search$multistart_summary$start_id,
    set_size = recursive_search$multistart_summary$final_adjustment_size,
    total_confidence = recursive_search$multistart_summary$total_path_confidence,
    path_edge_mean_confidence = recursive_search$multistart_summary$path_edge_mean_confidence,
    final_valid = recursive_search$multistart_summary$final_valid,
    stringsAsFactors = FALSE
  )
  recommended_res <- list(
    recommended_set_id = winner$start_id,
    recommended_variables = winner$final_set,
    recommended_table = recommended_table,
    set_scores = set_scores
  )
  node_meta_out <- merge(node_meta, node_map, by = "node", all.x = TRUE, sort = FALSE)
  recommended_metrics <- find_adjustment_dagitty_build_recommended_metrics(
    edges = dag_edges_out,
    tf_query = tf_query,
    gene_query = gene_query,
    recommended_variables = recommended_res$recommended_variables,
    beta = beta,
    direct_edges = edges
  )
  keep_forced <- recommended_metrics$tf %in% c(tf_query, gene_query)
  # A separator can be necessary even without a directed path to Y (e.g. a
  # collider-control parent). Never discard members of the validated set.
  keep_adjustment_path <- recommended_metrics$tf %in% recommended_res$recommended_variables
  filtered_recommended_metrics <- recommended_metrics[
    keep_forced | keep_adjustment_path, , drop = FALSE
  ]
  filtered_recommended_metrics$meets_prior_edge_threshold <-
    keep_forced[keep_forced | keep_adjustment_path] |
    ifelse(is.finite(filtered_recommended_metrics$distance_to_gene_A) &
             filtered_recommended_metrics$distance_to_gene_A > 0,
           filtered_recommended_metrics$avg_confidence_to_gene_A,
           filtered_recommended_metrics$avg_confidence_to_tf_B) >= overall_confidence_threshold
  if (any(!filtered_recommended_metrics$meets_prior_edge_threshold)) {
    stop(
      "Internal error: an adjustment path violates the prior-edge confidence threshold.",
      call. = FALSE
    )
  }
  if (!all(winner$final_set %in% filtered_recommended_metrics$tf)) {
    stop("Internal error: adjustment members were lost from regression metrics.", call. = FALSE)
  }
  rownames(filtered_recommended_metrics) <- NULL
  node_meta_out$role <- "other"
  node_meta_out$role[node_meta_out$node == tf_query] <- "exposure"
  node_meta_out$role[node_meta_out$node == gene_query] <- "outcome"
  node_meta_out$role[node_meta_out$node %in% candidate_df$variable] <- "candidate_confounder"
  node_meta_out$role[node_meta_out$node %in% recommended_res$recommended_variables] <- "recommended_adjustment"
  node_meta_out$role[node_meta_out$node == tf_query] <- "exposure"
  node_meta_out$role[node_meta_out$node == gene_query] <- "outcome"

  file_map <- list(
    dag_edges = file.path(outdir, sprintf("dagitty_local_dag_edges_%s_%s.csv", tf_query, gene_query)),
    dag_nodes = file.path(outdir, sprintf("dagitty_local_dag_nodes_%s_%s.csv", tf_query, gene_query)),
    candidate_confounders = file.path(outdir, sprintf("dagitty_candidate_confounders_%s_%s.csv", tf_query, gene_query)),
    minimal_adjustment_sets = file.path(outdir, sprintf("dagitty_minimal_adjustment_sets_%s_%s.csv", tf_query, gene_query)),
    recommended_adjustment_set = file.path(outdir, sprintf("dagitty_recommended_adjustment_set_%s_%s.csv", tf_query, gene_query)),
    recommended_adjustment_metrics = file.path(outdir, sprintf("dagitty_recommended_adjustment_metrics_%s_%s.csv", tf_query, gene_query)),
    recommended_adjustment_metrics_filtered = file.path(
      outdir,
      sprintf("dagitty_recommended_adjustment_metrics_filtered_%s_%s.csv", tf_query, gene_query)
    ),
    adjustment_set_scores = file.path(outdir, sprintf("dagitty_adjustment_set_scores_%s_%s.csv", tf_query, gene_query)),
    randomized_search_summary = file.path(outdir, sprintf("recursive_adjustment_multistart_%s_%s.csv", tf_query, gene_query)),
    recursive_operations = file.path(outdir, sprintf("recursive_adjustment_operations_%s_%s.csv", tf_query, gene_query)),
    dagitty_graph = file.path(outdir, sprintf("dagitty_graph_%s_%s.txt", tf_query, gene_query)),
    summary = file.path(outdir, sprintf("dagitty_summary_%s_%s.txt", tf_query, gene_query))
  )

  if (write_files) {
    utils::write.csv(filtered_recommended_metrics, file = file_map$recommended_adjustment_metrics_filtered, row.names = FALSE)
    if (isTRUE(write_full_outputs)) {
      utils::write.csv(dag_edges_out, file = file_map$dag_edges, row.names = FALSE)
      utils::write.csv(node_meta_out, file = file_map$dag_nodes, row.names = FALSE)
      utils::write.csv(candidate_df, file = file_map$candidate_confounders, row.names = FALSE)
      utils::write.csv(adj_sets_df, file = file_map$minimal_adjustment_sets, row.names = FALSE)
      utils::write.csv(recommended_res$recommended_table, file = file_map$recommended_adjustment_set, row.names = FALSE)
      utils::write.csv(recommended_metrics, file = file_map$recommended_adjustment_metrics, row.names = FALSE)
      utils::write.csv(recommended_res$set_scores, file = file_map$adjustment_set_scores, row.names = FALSE)
      utils::write.csv(recursive_search$multistart_summary, file = file_map$randomized_search_summary, row.names = FALSE)
      utils::write.csv(winner$operations, file = file_map$recursive_operations, row.names = FALSE)
      writeLines(dagitty_string, con = file_map$dagitty_graph)
    }
  }

  summary_lines <- c(
    sprintf("TF: %s", tf_query),
    sprintf("Gene: %s", gene_query),
    sprintf("Edge file: %s", edge_file),
    sprintf("Prior edges before confidence threshold: %d", edges_before_threshold),
    sprintf("Prior edges at confidence >= %s: %d", format(overall_confidence_threshold, trim = TRUE), edges_after_threshold),
    sprintf("Local DAG nodes: %d", length(local_nodes)),
    sprintf("Local DAG edges: %d", nrow(dag_edges_out)),
    sprintf("Candidate confounders: %d", nrow(candidate_df)),
    sprintf("Adjustment algorithm: %s", recursive_search$algorithm),
    sprintf("Minimum cardinality certified: %s", recursive_search$optimality_certified),
    "Confidence is reported, not optimized among equal-size minimum cuts.",
    sprintf("Solver runs: %d", recursive_search$n_starts),
    sprintf("Search cores: %d", recursive_search$cores),
    sprintf("Winning start: %d", winner$start_id),
    sprintf("Winning adjustment size: %d", length(winner$final_set)),
    sprintf("Winning path-edge mean confidence: %.6f", recursive_search$winner_score$path_edge_mean_confidence),
    sprintf("Search elapsed seconds: %.3f", recursive_search$search_seconds),
    sprintf("Beta: %s", format(beta, trim = TRUE)),
    sprintf("Overall confidence threshold: %s", format(overall_confidence_threshold, trim = TRUE)),
    sprintf(
      "Recommended adjustment set: %s",
      if (length(recommended_res$recommended_variables) == 0) "none" else paste(recommended_res$recommended_variables, collapse = ";")
    ),
    sprintf("DAG edge file: %s", file_map$dag_edges),
    sprintf("Adjustment set file: %s", file_map$minimal_adjustment_sets),
    sprintf("Recommended set file: %s", file_map$recommended_adjustment_set),
    sprintf("Recommended metrics file: %s", file_map$recommended_adjustment_metrics),
    sprintf("Filtered recommended metrics file: %s", file_map$recommended_adjustment_metrics_filtered)
  )

  if (write_files && isTRUE(write_full_outputs)) {
    writeLines(summary_lines, con = file_map$summary)
  }

  result <- list(
    tf = tf_query,
    gene = gene_query,
    network = network,
    edge_file = edge_file,
    prior_edge_confidence_threshold = overall_confidence_threshold,
    prior_edge_counts = c(
      before = edges_before_threshold,
      after = edges_after_threshold
    ),
    dag_edges = dag_edges_out,
    dag_nodes = node_meta_out,
    candidate_confounders = candidate_df,
    minimal_adjustment_sets = adj_sets_df,
    recommended_adjustment_set = recommended_res$recommended_table,
    recommended_adjustment_metrics = recommended_metrics,
    filtered_recommended_adjustment_metrics = filtered_recommended_metrics,
    adjustment_set_scores = recommended_res$set_scores,
    randomized_search = recursive_search,
    adjustment_search = recursive_search,
    dagitty_graph = dagitty_string,
    files = file_map,
    summary = summary_lines
  )
  class(result) <- c("find_adjustment_dagitty_result", class(result))
  result
}

find_adjustment_dagitty <- function(
  tf,
  gene,
  network = NULL,
  edge_file = "tf_union_output/tf_gene_merged_human_weighted_clean.csv",
  outdir = "adjustment_output",
  beta = 2,
  overall_confidence_threshold = 3,
  search_starts = 1L,
  search_cores = 4L,
  search_seed = 123L,
  max_adjustment_sets = NULL,
  include_candidate_confounders = TRUE,
  write_full_outputs = FALSE,
  write_files = TRUE
) {
  find_adjustment_dagitty_run(
    tf = tf,
    gene = gene,
    network = network,
    edge_file = edge_file,
    outdir = outdir,
    beta = beta,
    overall_confidence_threshold = overall_confidence_threshold,
    search_starts = search_starts,
    search_cores = search_cores,
    search_seed = search_seed,
    max_adjustment_sets = max_adjustment_sets,
    include_candidate_confounders = include_candidate_confounders,
    write_full_outputs = write_full_outputs,
    write_files = write_files
  )
}

print.find_adjustment_dagitty_result <- function(x, ...) {
  cat(paste(x$summary, collapse = "\n"), "\n")
  invisible(x)
}

find_adjustment_dagitty_main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  required_args <- c("--tf", "--gene")
  missing_required <- required_args[!vapply(required_args, function(x) any(grepl(paste0("^", x, "="), args)), logical(1))]
  if (length(missing_required) > 0) {
    stop(sprintf("Missing required argument(s): %s", paste(missing_required, collapse = ", ")))
  }

  result <- find_adjustment_dagitty(
    tf = find_adjustment_dagitty_get_arg(args, "--tf"),
    gene = find_adjustment_dagitty_get_arg(args, "--gene"),
    edge_file = find_adjustment_dagitty_get_arg(args, "--edge_file", "tf_union_output/tf_gene_merged_human_weighted_clean.csv"),
    outdir = find_adjustment_dagitty_get_arg(args, "--outdir", "adjustment_output"),
    beta = as.numeric(find_adjustment_dagitty_get_arg(args, "--beta", "2")),
    overall_confidence_threshold = as.numeric(find_adjustment_dagitty_get_arg(args, "--overall_confidence_threshold", "3")),
    search_starts = as.integer(find_adjustment_dagitty_get_arg(args, "--search_starts", "1")),
    search_cores = as.integer(find_adjustment_dagitty_get_arg(args, "--search_cores", "4")),
    search_seed = as.integer(find_adjustment_dagitty_get_arg(args, "--search_seed", "123")),
    write_full_outputs = tolower(find_adjustment_dagitty_get_arg(args, "--write_full_outputs", "false")) == "true",
    write_files = TRUE
  )

  print(result)
}

if (sys.nframe() == 0) {
  find_adjustment_dagitty_main()
}
